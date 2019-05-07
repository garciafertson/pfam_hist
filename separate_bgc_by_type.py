'''the script separates bgc.gbk files by type 
   reads in a file with the list of genebank filenames with absolute
   path, recovers the product type in cluster FEATURE and create a new file for
   each type found with the filenames'''
#/usr/local/bin/python3.7
import sys
from collections import defaultdict
from Bio import SeqIO

filename = sys.argv[1]
type_dict = defaultdict(list)

with open(filename, 'r') as f:
    for line in f:
        name=line.strip()
        with open(name, 'r') as gbk:
            for record in SeqIO.parse(gbk, "gb"):
                for feature in record.features:
                    if feature.type=="cluster":
                        assert len(feature.qualifiers['product'])>0
                        bgc_type=feature.qualifiers['product'][0]
                        type_dict[bgc_type].append(name)
                        break

for bgc_type in type_dict.keys():
    outfile=bgc_type+".txt"
    with open (outfile, "w") as out:
        for name in type_dict[bgc_type]:
            out.write("%s\n" %name )


