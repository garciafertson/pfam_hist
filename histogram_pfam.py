import sys, argparse
from collections import Counter
import csv
import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns
import pandas as pd
#import Image
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", dest="input", required=True,
                      help="input .csv with pfam list", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", required=True,
                      help="Output for tsv count table", metavar="FILE")
    options = parser.parse_args()

    with open(options.input) as domtbl:
        reader=csv.reader(domtbl, delimiter=',')
        counter=Counter(line[6] for line in reader if line[6])
    _,values=zip(*counter.items())

    #####
    #plot histogram, bin by frequencies, *y log values, 
    kwargs=dict(density=True, stacked=True, alpha=0.5)
#    geom=np.random.geometric(p=1/2.718, size=5000)
#    pois=np.random.poisson(lam=1, size=5000)
    plt.hist(values,bins=max(values), **kwargs)
#    plt.hist(geom,bins=max(geom), **kwargs, color="r")
#    plt.hist(pois,bins=max(pois), **kwargs, color='g' )
    v=np.array(values)
    per=np.percentile(v,95)
    plt.axis([0, 100, 0, 0.5])
    plt.title("Pfam histogram")
    plt.xlabel("times found")
    plt.ylabel("density")
    plt.vlines(per,ymin=0,ymax=1,color='r', 
            lw=0.3, label="95th percentile")
   # plt.show()
    figure_name=options.output+"_pfam_histogram.png"
    plt.savefig(figure_name)
    plt.close()

    #######
    #plot top 20 most frequent pfam
    labels,values=zip(*counter.most_common(20))
    xs= [i+0.1 for i, _ in enumerate(labels)]
    plt.bar(xs,values)
    plt.xticks([i for i,_ in enumerate(labels)], labels,
            size=8, rotation=45)
    plt.ylabel("frequency")
    plt.title("20 Most common PFAMs")
    figure_name=options.output+"_top20_Pfams.png"
    plt.savefig(figure_name)

    #save tsv with 95th pecentile most frequent pfams
    labelsper= [label for label in counter.keys() if counter[label]>per ]
    tablename=options.output+".csv"
    with open(tablename, 'w') as out:
        for label in labelsper:
            out.write("%s\t%s\n" %(label, counter[label]))


