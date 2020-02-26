#--coding:utf-8
'''
Created on 2017年12月8日
@author: Zhu Qianshu
'''
import sys,getopt
import numpy as np
import gzip 

resolution = 100000
size_file = "/home/qszhu/ann/hic-pro/chrom_mm10.sizes"

if __name__ == '__main__':
    ################usage###################
    usage = "Usage: " + sys.argv[0]
    usage += """
    Usage:python 
        <required>:
            -i input file.#chr1    3050000    chr1    3150000    85    1.659249e-05    6.992300e-03
            -s size file
            -o output file
            -r resolution.
    
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
        
        

    chrSize = {} # chr: [start_point,start_step]
    start = 1
    for line in open(size_file,"rU"):
        fs = line.rstrip().rsplit("\t")
        chrSize[fs[0]] = start
#         print fs[0],start
        start += int(np.ceil(float(fs[1])/resolution))
        
    print chrSize
    if input_file.endswith(".gz"):
        infile = gzip.GzipFile(input_file, "r").readlines()
    else:
        infile = open(input_file, "r").readlines()
    outfile = open(output_file,"w")
    percentage_finish = 0
    for i in range(1,len(infile)):
        if (i %  (len(infile) / 10)== 0):
            print 'finish ', percentage_finish, '%'
            percentage_finish += 10
        fs = infile[i].rstrip().rsplit("\t")
        s1 = chrSize[fs[0]] + int(np.floor(float(int(fs[1]))/resolution))
        s2 = chrSize[fs[2]] + int(np.floor(float(int(fs[3]))/resolution))
#         if fs[0] == "chr14" and fs[1]=="3050000" and fs[3]=="3150000":
#             print chrSize[fs[0]],fs,s1,s2
        outfile.write(str(s1) + "\t" + str(s2) + "\t" + fs[4] + "\t" + fs[5] + "\t" + fs[6] + "\n")
    
    outfile.close()
    
    
    
    
   
        
