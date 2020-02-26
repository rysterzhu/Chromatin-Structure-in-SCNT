#--coding:utf-8
'''
Created on 2018年12月2日
@author: Zhu Qianshu
2018年12月21日： add loop count
'''
import sys,getopt
import numpy as np
import gzip 

resolution = 40000
digest_file = "/home/qszhu/ann/hic-pro/MboI_resfrag_mm10.bed"
chrom = "chr1"

if __name__ == '__main__':
    ################usage###################
    usage = "Usage: " + sys.argv[0]
    usage += """
    Usage:python 
        <required>:
            -i input file.#chr1    3050000    chr1    3150000    85    1.659249e-05    6.992300e-03
            -f fithic file. #chr1    fragmentMid1    chr2    fragmentMid2    contactCount    p-value    q-valu
            -d digest file
            -o output file
            -r resolution.
            -c chromosome.
    
        [Options]:
    
    """
    try:
        opts,args = getopt.getopt(sys.argv[1:], "hi:o:f:d:r:c:", 
                                  ["help"])
    except getopt.GetoptError:
        sys.exit(usage)
    
    for opt in opts:
        if opt[0] == '-h': sys.exit(usage)
        elif opt[0] == '-i': input_file = opt[1]
        elif opt[0] == '-o': output_file = opt[1]
        elif opt[0] == '-f': fithic_file = opt[1]
        elif opt[0] == '-d': digest_file = opt[1]
        elif opt[0] == "-r": resolution = int(opt[1])
        elif opt[0] == '-c': chrom = opt[1]
    
    digest = {} #name:(start,end)
    for line in open(digest_file,"rU"):
        if not line.startswith(chrom):
            continue
        fs = line.rstrip().rsplit("\t")
        digest[fs[3]] = (int(fs[1]),int(fs[2]))
    
    if fithic_file.endswith(".gz"):
        infile = gzip.GzipFile(fithic_file, "r").readlines()
    else:
        infile = open(fithic_file, "rU").readlines()
    loops = {}    
    count=0
    for i in range(1,len(infile)):
        fs = infile[i].rstrip().rsplit("\t")
        if not (fs[0]==chrom and fs[2]==chrom):
            continue
            
        start = int(fs[1])
        end = int(fs[3])
        count = int(fs[4])
        
        if start>end: start,end = end,start
        k=str(start/resolution) + "-" + str(end/resolution)
        loops[k] = count
    
    outfile = open(output_file,"w")
    outfile.write("chr1\tstart1\tend1\tstrand1\tchr2\tstart2\tend2\tstrand2\tid\tloop\tlcount\n")
    if input_file.endswith(".gz"):
        infile = gzip.GzipFile(input_file, "r").readlines()
    else:
        infile = open(input_file, "rU").readlines()
    percentage_finish = 0
    for i in range(1,len(infile)):
        if (i %  (len(infile) / 10)== 0):
            #print 'write contacts: ', percentage_finish, '%'
            percentage_finish += 10
        fs = infile[i].rstrip().rsplit("\t")
        if not (fs[1]==chrom and fs[4]==chrom):
            continue
        
        id = fs[0]
        start = int(fs[2])
        strandS = fs[3]
        end = int(fs[5])
        strandE = fs[6]
        resS = fs[8]
        resE = fs[9]
        
        if start>end: start,end = end,start
        k=str(start/resolution) + "-" + str(end/resolution) #will be int to string
        if loops.has_key(k):    #the contacts in loops
            if strandS == "+":     #if the contact strand is +, then the fragment is position to the digest site of downstream
                start2 = digest[resS][1] 
            else:
                start2 = start #if the contact strand is -, then the fragment is the digest site of upstream to the position
                start = digest[resS][0] 
            if strandE == "+":
                end2 = digest[resE][1]
            else:
                end2 = end
                end = digest[resE][0]
            if start2 < start or end2 < end:
                print "error: " + str(start) + "\t" + str(start2) + "\t" + "\t" + str(end) + "\t" + str(end2) + "\t" + id + "\t" + k
            
            outfile.write(chrom + "\t" + str(start) + "\t" + str(start2) + "\t" + strandS + "\t" +
                          chrom + "\t" + str(end) + "\t" + str(end2) + "\t" + strandS + "\t" +
                          id + "\t" + k + "\t" + str(loops[k]) + "\n")
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        