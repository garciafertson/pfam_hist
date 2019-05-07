'''the script takes in the csv in proteins2fasta fortmat with pfam list
   reads line by line and process 2x2 contingency table for pairs of most frequent BGC
   (uper 70 percentil) only work with the top BGC composing the 30% of the total pfam
   recieve a file or the list with those top (most frequent) pfam
   en al menos el 10% del total de BGC
'''
import argparse
import numpy as np
from collections import defaultdict
class bgc_pfam(object):
    def __init__(self,contig_id):
        self.contig_id=contig_id
        self.pfam=[]
        self.locustag=[]
        self.g_pos_start=[]
#        self.g_pos_end=[]
        self.p_pos_start=[]
        self.p_pos_end=[]
        self.bp_len=0
        self.pfam_set=set()
        self.posdict={}

    def add_pfam_from_str(self, string):
        line=string.split(",")
        self.contig_id=line[0]
        self.pfam.append(line[6])
        self.g_pos_start.append(line[3])
        self.g_pos_end.append(line[4])
        self.p_pos_start.append(line[7])
        self.p_pos_end.append(line[8])

    def calulcate_bgc_len(self):
        self.bp_len=min(self.g_pos_start)-max(self.g_pos_end)

    def calulate_pfam_set(self):
        self.pfam_set=set(self.pfam)

    def add_bgc_cont_matrix(self, glob_pfam_set):
        pfam_list=list(glob_pfam_set)
        dim_pfam=len(pfam_list)
        contignency_add=np.zeros((dim_pfam,dim_pfam,2,2))
        for i, pfami in enumerate(pfam_list):
            for j in range(i+1,dim_pfam):
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
                else:
                    for k,p1 in enumerate(self.pfam):
                        if pfami = p1:
                            posg=

                        if()# if in different gene add one paired pressence
                        if() # else same gene, if in different dommain add one paired presence
                             # else same gene,same domain, do nothing, a
                


