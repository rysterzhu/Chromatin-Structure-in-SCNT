#--coding:utf-8
'''
Created on 2017年12月8日
@author: Zhu Qianshu
'''
import sys,getopt
import numpy as np
import gzip 

resolution = 10000
size_file = "/home/qszhu/ann/hic-pro/chrom_mm10.sizes"

################usage########################################################################
usage = "Usage: " + sys.argv[0]
usage += """
Usage:python 
    <required>:
        -i interaction raw file
        -b compartment PC1 bed file
        -s chromosome size file
        -o output file
        -r resolution.
        -c short distance cut off.

    [Options]:

"""
try:
    opts,args = getopt.getopt(sys.argv[1:], "hi:o:s:r:", 
                              ["help"])
except getopt.GetoptError:
    sys.exit(usage)

for opt in opts:
    if opt[0] == '-h': sys.exit(usage)
    elif opt[0] == '-i': input_file = opt[1]
    elif opt[0] == '-o': output_file = opt[1]
    elif opt[0] == '-s': size_file = opt[1]
    elif opt[0] == "-r": resolution = int(opt[1])
###########################################################################################   
    
chrSize = {} # chr: [start_step,end_step]
chrProbability = {} # # chr1:{0:...,1:...}
start = 1
for line in open(size_file,"rU"):
    fs = line.rstrip().rsplit("\t")
    chrSize[fs[0]] = [start,start+int(np.ceil(float(fs[1])/resolution))-1]
    chrProbability[fs[0]] = {}
    start += int(np.ceil(float(fs[1])/resolution))       
print chrSize

with open(input_file,"r") as f:
    for line in f:
        fs = line.rstrip().rsplit("\t")
        f1,f2,f3 = int(fs[0]),int(fs[1]),int(fs[2])
        
        for c in chrSize.keys():
            if f1 >= chrSize[c][0] and f1 <= chrSize[c][1]: c1 = c
            if f2 >= chrSize[c][0] and f2 <= chrSize[c][1]: c2 = c
        if c1 != c2: continue
        
        if chrProbability[c1].has_key(f2-f1):
            chrProbability[c1][f2-f1] += f3
        else:
            chrProbability[c1][f2-f1] = f3
            
            
outfile = open(output_file,"w")
for c in chrProbability.keys():
    for l in sorted(chrProbability[c].keys()):
        outfile.write(c + "\t" + str(l) + "\t" + str(chrProbability[c][l]) + "\n")





