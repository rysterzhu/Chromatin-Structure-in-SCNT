#--coding:utf-8
'''
Created on 2018年12月2日
@author: Zhu Qianshu
2018年12月21日： add loop count
'''
import sys,getopt
import numpy as np
import gzip 

digest_file = "/home/qszhu/ann/hic-pro/MboI_resfrag_mm10.bed"
chrom = "chr1"
fdr=500
if __name__ == '__main__':
    ################usage###################
    usage = "Usage: " + sys.argv[0]
    usage += """
    Usage:python 
        <required>:
            -i input file.#chr1    3050000    chr1    3150000    85    1.659249e-05    6.992300e-03
            -d digest file
            -o output file #chr1    start1    end1    strand1    chr2    start2    end2    strand2    id    loop    lcount
            -r output file for density
            -f fdr for degest length
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
        elif opt[0] == '-r': output_file2 = opt[1]
        elif opt[0] == '-d': digest_file = opt[1]
        elif opt[0] == '-c': chrom = opt[1]
        elif opt[0] == "-f": fdr = int(opt[1])
    
    digest = {} #name:(start,end)
    for line in open(digest_file,"rU"):
        if not line.startswith(chrom):
            continue
        fs = line.rstrip().rsplit("\t")
        digest[fs[3]] = (int(fs[1]),int(fs[2]))
    
    
    outfile = open(output_file,"w")
    outfile.write("chr1\tstart1\tend1\tstrand1\tchr2\tstart2\tend2\tstrand2\tid\tloop\tlcount\n")
    outfile2 = open(output_file2,"w")
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
        if strandS == "+":     #if the contact strand is +, then the fragment is position to the digest site of downstream
            start2 = digest[resS][1] 
            if fdr > 0 and start2-start > fdr: start=start2-fdr
        else:
            start2 = start #if the contact strand is -, then the fragment is the digest site of upstream to the position
            start = digest[resS][0] 
            if fdr > 0 and start2-start > fdr: start2=start+fdr
        if strandE == "+":
            end2 = digest[resE][1]
            if fdr > 0 and end2-end > fdr: end=end2-fdr
        else:
            end2 = end
            end = digest[resE][0]
            if fdr > 0 and end2-end > fdr: end2=end+fdr
         
        outfile.write(chrom + "\t" + str(start) + "\t" + str(start2) + "\t" + strandS + "\t" +
                      chrom + "\t" + str(end) + "\t" + str(end2) + "\t" + strandS + "\t" +
                      id + "\n")
        outfile2.write(str(int((end2+end-start2-start)/2))+"\n")
    outfile.close()
    outfile2.close()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        