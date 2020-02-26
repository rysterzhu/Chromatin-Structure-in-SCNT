#--coding:utf-8
'''
Created on 2018年6月5日

@author: tj1
'''
import sys,getopt
import numpy as np

resolution = 10000
if __name__ == '__main__':
    ################usage###################
    usage = "Usage: " + sys.argv[0]
    usage += """
    Usage:python 
        <required>:
            -i genetrack index file.[chrom    index    forward    reverse]
            -I input index file
            -o output file
            -c chrom.
    
        [Options]:
    
    """
    try:
        opts,args = getopt.getopt(sys.argv[1:], "hi:o:s:c:", 
                                  ["help"])
    except getopt.GetoptError:
        sys.exit(usage)
    
    for opt in opts:
        if opt[0] == '-h': sys.exit(usage)
        elif opt[0] == '-i': input_file = opt[1]
        elif opt[0] == '-o': output_file = opt[1]
        elif opt[0] == '-s': size_file = opt[1]
        elif opt[0] == "-c": chrN = opt[1]

    start = 1
    for line in open(size_file,"rU"):
        fs = line.rstrip().rsplit("\t")
        chrom = fs[0]
        length = int(np.ceil(float(fs[1])/resolution))
        if chrom == chrN:
            break
        start += length
    
    lines = open(input_file,"rU").readlines()
    out = open(output_file,'w')
    
    print start,"-",len(lines),"==",length
    if len(lines) != length: 
        print "error"
    else:
        for i in range(length):
            fs = lines[i].rsplit("\t")
            for j in range(i,length):
                if float(fs[j]) != 0.0:
                    out.write(str(start+i)+'\t'+str(start+j)+"\t"+fs[j]+"\n")
    out.close()
    
    
    
    
    
    
    
    