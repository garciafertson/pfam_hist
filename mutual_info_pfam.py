'''the script takes in the csv in proteins2fasta fortmat with pfam list
   reads line by line and process 2x2 contingency table for pair of PFAMs,
   only most frequent PFAMs in BGC set, top 95th pecentile.
   *Segundo filtro, los pfams se encuentran el al menos en al menos 
   el 5% del total de BGC, ~freq/totalBGCnimber > 0.05
'''
import argparse
import numpy as np
from collections import defaultdict
import csv
import math
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
from scipy.cluster import hierarchy
import pandas as pd
from numpy import linalg as LA

def getQuasiDiag(link):
    # Sort clustered items by distance
    link = link.astype(int)
    sortIx = pd.Series([link[-1, 0], link[-1, 1]])
    numItems = link[-1, 3]  # number of original items
    while sortIx.max() >= numItems:
        sortIx.index = range(0, sortIx.shape[0] * 2, 2)  # make space
        df0 = sortIx[sortIx >= numItems]  # find clusters
        i = df0.index
        j = df0.values - numItems
        sortIx[i] = link[j, 0]  # item 1
        df0 = pd.Series(link[j, 1], index=i + 1)
        sortIx = sortIx.append(df0)  # item 2
        sortIx = sortIx.sort_index()  # re-sort
        sortIx.index = range(sortIx.shape[0])  # re-index
    return sortIx.tolist()

class bgc_pfam(object):
    def __init__(self,contig_id):
        self.contig_id=contig_id
        self.pfam=[]
        self.locustag=[]
        self.g_pos_start=[]
        self.g_pos_end=[]
        self.p_pos_start=[]
        self.p_pos_end=[]
        self.bp_len=0
        self.pfam_set=set()
        self.posdict={}

    def add_pfam_from_line(self, line):
        self.contig_id=line[0]
        self.pfam.append(line[6])
        self.locustag.append(line[1])
        self.g_pos_start.append(int(line[3]))
        self.g_pos_end.append(int(line[4]))
        self.p_pos_start.append(int(line[7]))
        self.p_pos_end.append(int(line[8]))

    def calculate_bgc_len(self):
        self.bp_len=max(self.g_pos_end)-min(self.g_pos_start)
        return(self.bp_len)

    def calculate_pfam_set(self):
        self.pfam_set=set(self.pfam)

    def pfam_count(self, pfam_list):
        dim_pfam=len(pfam_list)
        contingency_add=np.zeros((dim_pfam,dim_pfam,2,2))
        for i, pfami in enumerate(pfam_list):
            for j in range(i,dim_pfam):
                pfamj=pfam_list[j]
                if pfami not in self.pfam_set and pfamj not in self.pfam_set:
                    contingency_add[i][j][0][0]+=1
                    contingency_add[j][i][0][0]+=1
                elif pfami not in self.pfam_set and pfamj in self.pfam_set:
                    contingency_add[i][j][0][1]+=1
                    contingency_add[j][i][1][0]+=1
                elif pfami in self.pfam_set and pfamj not in self.pfam_set:
                    contingency_add[i][j][1][0]+=1
                    contingency_add[j][i][0][1]+=1
                #for now do not consider if the pfam 
                #correspond to the same sequence
                else:
                    contingency_add[i][j][1][1]+=1
                    contingency_add[j][i][1][1]+=1
        return contingency_add

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", dest="input", required=True,
                      help="input .csv with pfam list", metavar="FILE")
parser.add_argument("-p", "--pfam", dest="pfam", required=True,
                      help="most freq pfam list", metavar="FILE")
options = parser.parse_args()

pfam_list=[]
freq_list=[]
bgc_number=0

#fill in list of most frequent pfams in bgc set
with open(options.pfam) as f:
    for line in f:
        line=line.strip()
        pfam,freq=line.split("\t")
        pfam_list.append(pfam)
        freq_list.append(freq)
#print(pfam_list)
#contig_id  0,  locus_tag    1,    protein_id  2,
#gene_start 3,  gene_end     4,    gene_strand 5,
#pfam_id    6,  domain_start 7,    domain_end  8,
#evalue     9,  bitscore    10,

#contingency tensor, presence/absence pairs of most freq pfams
pfam_cont=np.zeros([len(pfam_list),len(pfam_list),2,2])
#mutual information matrix, pairs of most freq pfams
pfam_MI=np.zeros([len(pfam_list), len(pfam_list)])
pfam_Fet=np.zeros([len(pfam_list), len(pfam_list)])
#list of bp length for each BGC in complete set
bgc_lengths=[]

prev_id=''
#read file with PFAM list from BGC set and fill contingency tensor
with open(options.input) as f:
    reader=csv.reader(f, delimiter=",")
    next(reader)
    for line in reader:
        bgc_number+=1
        if prev_id == line[0]:
            new_contig_id=False
        else:
            new_contig_id=True
        if new_contig_id:
            try:
                bgc
            except:
                bgc=bgc_pfam(line[0])
                bgc.add_pfam_from_line(line)
                continue
            #print(bgc.contig_id)
            bp_len=bgc.calculate_bgc_len()
            bgc_lengths.append(bp_len)
            bgc.calculate_pfam_set()
            #print(bgc.pfam_set)
            add_pfam_freq=bgc.pfam_count(pfam_list)
            pfam_cont += add_pfam_freq
            #print(pfam_cont[1][2])
            bgc= bgc_pfam(line[0])
            bgc.add_pfam_from_line
        else:
            bgc.add_pfam_from_line(line)
        prev_id=line[0]


#### calculate Mutual Infromation Matrix from contingency tensor
for i in range(len(pfam_list)):
    for j in range(i, len(pfam_list)):
        total=float(pfam_cont[i][j].sum())
        Fy=pfam_cont[i][j].sum(axis=0)
        Fx=pfam_cont[i][j].sum(axis=1)
        MI=0
        for x in (0,1):
            for y in (0,1):
                pxy=pfam_cont[i][j][x][y]/total
                px=Fx[x]/total
                py=Fy[y]/total
                ratio=pxy/(px*py)
                if ratio>0: MI += pxy * math.log(ratio,2)
                else: continue
        pfam_MI[i][j]=MI
        pfam_MI[j][i]=MI
#        _,pvalue=stats.fisher_exact(pfam_cont[i][j])
#        pfam_Fet[i][j]=pvalue
#        pfam_Fet[j][i]=pvalue

threshold=0.1
distance=(pfam_MI)**.5
linkage=hierarchy.linkage(distance, method='complete')
sortIx=getQuasiDiag(linkage)

w,v=LA.eig(distance)
print(len(w))
print (w.shape)
print (v.shape)
#w,v=LA.eig(pfam_MI)
#lmd=np.asmatrix(np.diag(w))
#ev=np.asmatrix(v)
cl=np.zeros([len(w),len(w)])
print(cl.shape)
print(w[0:20])

#iw,v=LA.eig(pfam_MI)
for i in range(len(w)):
    for j in range(len(w)):
        tmp=v[i,0:20]*w[0:20]
        cl[i,j]=np.dot(tmp, v[j,0:20])

print(cl)
cl=np.asarray(cl)
MI_sort=np.zeros((len(sortIx),len(sortIx)))
for i,a in enumerate(sortIx):
    for j,b in enumerate(sortIx):
        MI_sort[i][j]=pfam_MI[a][b]

#clusters=hierarchy.fcluster(linkage, threshold, criterion="distance")
plt.subplot(121)
plt.imshow(cl,cmap='hot', interpolation='nearest')
plt.subplot(122)
hierarchy.dendrogram(linkage, color_threshold=threshold)
plt.show()

#plt.show()
#print(clusters)
#print (pfam_MI)

