#--coding:utf-8
'''
modify: 2018年10月18日 add PC1 cutoff
Created on 2017年12月8日
@author: Zhu Qianshu
'''
import sys,getopt
import numpy as np
import gzip 

resolution = 25000
size_file = "/home/qszhu/ann/hic-pro/chrom_mm10.sizes"
shortCutoff = 2e6
pc1Cutoff = 0.0

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
        -c short distance cut off. [default 2M]
        -p pc1 value cut off for compartment [default 0]
        

    [Options]:

"""
try:
    opts,args = getopt.getopt(sys.argv[1:], "hi:b:o:s:r:c:p:", 
                              ["help"])
except getopt.GetoptError:
    sys.exit(usage)

for opt in opts:
    if opt[0] == '-h': sys.exit(usage)
    elif opt[0] == '-i': input_file = opt[1]
    elif opt[0] == '-b': bg_file = opt[1]
    elif opt[0] == '-o': output_file = opt[1]
    elif opt[0] == '-s': size_file = opt[1]
    elif opt[0] == "-r": resolution = int(opt[1])
    elif opt[0] == "-c": shortCutoff = float(opt[1])
    elif opt[0] == "-p": pc1Cutoff = float(opt[1])
###########################################################################################   
    
chrSize = {} # chr: [start_step,end_step]
chrDiff = {} # chr: [diff_compartment_bins,diff_compartment_interactions,same_compartment_bins,same_compartment_interactions]
start = 1
for line in open(size_file,"rU"):
    fs = line.rstrip().rsplit("\t")
    chrSize[fs[0]] = [start,start+int(np.ceil(float(fs[1])/resolution))-1]
    chrDiff[fs[0]] = [0,0,0,0]
    start += int(np.ceil(float(fs[1])/resolution))       
print chrSize

pc1 = {}
firstL = True
with open(bg_file, "r") as f:
    for line in f:
        if firstL:
            firstL = False
            continue
        fs = line.rstrip().rsplit("\t")
        s1 = chrSize[fs[0]][0] + int(np.floor(float(int(fs[1]))/resolution))
        p = float(fs[3])
        if np.abs(p) >= pc1Cutoff:
            pc1[s1] = p
print "read PC1 done."

with open(input_file,"r") as f:
    for line in f:
        fs = line.rstrip().rsplit("\t")
        f1,f2,f3 = int(fs[0]),int(fs[1]),int(fs[2])
        
        if np.abs(f2 - f1) < shortCutoff/resolution: continue
        
        for c in chrSize.keys():
            if f1 >= chrSize[c][0] and f1 <= chrSize[c][1]: c1 = c
            if f2 >= chrSize[c][0] and f2 <= chrSize[c][1]: c2 = c
        if c1 != c2: continue
        
        try:
            temp = pc1[f1] * pc1[f2] #has no key error
        except:
            continue
        
        if temp < 0: #如果这个interaction对应的两个位置的PC1值的乘积大于0，说明它们在相同类型的compartment内；若小于0，则在不同类型的compartment内
            chrDiff[c1][0] += 1
            chrDiff[c1][1] += f3
        elif temp > 0:
            chrDiff[c1][2] += 1
            chrDiff[c1][3] += f3
        else:
            print "warning PC1==0:", fs, pc1[f1], pc1[f2]

outfile = open(output_file,"w")
outfile.write("Chr\tDiff_bins\tDiff_interactions\tSame_bins\tSame_interactions\n")
for c in chrDiff.keys():
    outfile.write(c + "\t" + str(chrDiff[c][0]) + "\t" + str(chrDiff[c][1]) + "\t"
                   + str(chrDiff[c][2]) + "\t" + str(chrDiff[c][3]) +"\n" )
    
outfile.close()           
                
    
    
    