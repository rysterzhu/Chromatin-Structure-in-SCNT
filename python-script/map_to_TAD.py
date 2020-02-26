#--coding:utf-8
'''
Created on 2018年7月26日

@author: tj1
'''
import numpy as np
from skimage.transform import resize
import sys,getopt
import cPickle

bed_file = "40000_mm10.bed"
input_file = "test.matrix"
tad_file = "morula.filt.tad"
output_file = "RTI.test.bg"
resize_tad_length = 1.2e6
extend_length = 4e5
resolution = 40000
order = 3

def parse_bed(bed_file):
    step = {} #{chr:(first_sit,end_site),...}
    isF = True
    first = 0
    chrom = None
    for line in open(bed_file,"rU"):
        fs = line.rstrip().rsplit("\t")
        if 0 == int(fs[1]):
            if not isF:
                step[chrom] = (first,int(fs[3]) - 1)
            isF = False    
            chrom = fs[0]
            first = int(fs[3])
    
    step[chrom] = (first,int(fs[3]) - 1)
    return step
    
def parse_tad(tad_file):
    tad = {} #{chrom:[(start,end),...],...}
    for line in open(tad_file,"rU"):
        fs = line.rstrip().rsplit("\t")
        try:
            chrom = fs[0]
            start = int(fs[1])
            end = int(fs[2])
        except:
            continue
        tad.setdefault(chrom,[]).append((start,end))
        
    return tad
        
def parse_matrix(input_file,step):
    mat = {} # {chrom:{(start, end):value,...}, end > start, start:x end:y
    for line in open(input_file,"rU"):
        fs = line.rstrip().rsplit("\t")
        try:
            start = int(fs[0])
            end = int(fs[1])
            value = float(fs[2])
        except:
            continue
        for c in step.keys():
            if step[c][0] <= start <= step[c][1]:
                if step[c][0] <= end <= step[c][1]:
                    chrom = c
                else:
                    chrom = None          # filt the diff chromomosome interactions
                break
        if chrom != None:    
            mat.setdefault(chrom, {})[(start,end)] = value

    return mat

def extract_square_around_tad(mat,tad,step,resolution,resize_tad_length,extend_length,order = order):
    resize_box_length = int((resize_tad_length + 2 * extend_length)/resolution)
    output_box = np.zeros((resize_box_length,resize_box_length))
    isF=False
    n=0
    for chrom in tad.keys():
        for (tadS,tadE) in tad[chrom]:
            a = int(tadS/resolution) + step[chrom][0] # tad start 
            b = int(tadE/resolution) + step[chrom][0] # tad end
            if a > b: a,b = b,a
            d = b - a + 1 # tad length
            e = int(np.round(d*float(extend_length/resize_tad_length))) # extend length
            a = a - e
            b = b + e
            d = b - a + 1
            tadM = np.zeros((d,d))
            for i in range(d):
                for j in range(d):
                    try:
                        if i > j:
                            tadM[i,j] = mat[chrom][(a+j,a+i)]
                        else:
                            tadM[i,j] = mat[chrom][(a+i,a+j)]
                    except:
                        pass
            resize_tadM = resize(tadM, (resize_box_length,resize_box_length), order=order) # default order=3
            if isF:
                print chrom,tadS,tadE,a,b,d
                print tadM
                print resize_tadM
                isF = False
            output_box = output_box + resize_tadM
            n += 1
    return output_box/n


if __name__ == '__main__':
    ################usage###################
    usage = "Usage: " + sys.argv[0]
    usage += """
    Usage:python TAD_analysis.consolidation_of_tad
        <required>:
            -i HiC-Pro iced matirx.# [site site interaction]
            -I input tad bed file. # bed file
            -b annotation bed of iced matrix; # chr    start    end site
            -o output file  # chr    start    end    in_score    in_num    out_score    out_num    ratio
        [Options]:
            -r resolution of matrix. [default = 40k]
            -s resize tad length. [default 1.2 million ]
            -e extend tad length. [default 0.4 million ]
    """
    try:
        opts,args = getopt.getopt(sys.argv[1:], "hi:o:I:s:b:r:e:", 
                                  ["help"])
    except getopt.GetoptError:
        sys.exit(usage)
    
    for opt in opts:
        if opt[0] == '-h': sys.exit(usage)
        elif opt[0] == '-i': input_file = opt[1]
        elif opt[0] == '-o': output_file = opt[1]
        elif opt[0] == '-I': tad_file = opt[1]
        elif opt[0] == '-b': bed_file = opt[1]
        elif opt[0] == "-s": resize_tad_length = int(opt[1])
        elif opt[0] == "-e": extend_length = int(opt[1])
        elif opt[0] == "-r": resolution = int(opt[1])
        
    tad = parse_tad(tad_file)
    step = parse_bed(bed_file)
    mat = parse_matrix(input_file, step)
#     cPickle.dump(mat, open("test.dump","w"), protocol=0)
#     mat = cPickle.load(open("test.dump",'r'))
    output_box = extract_square_around_tad(mat, tad, step, resolution, resize_tad_length, extend_length, order)
#     np.save(open("output_box.np","w"),output_box)
    [x,y] = output_box.shape
    o = open(output_file,"w")
    for i in range(x):
        for j in range(y):
            o.write(str(output_box[i,j]) + "\t")
        o.write("\n")
    o.close()
    








