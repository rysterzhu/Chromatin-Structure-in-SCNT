#--coding:utf-8
'''
Created on 2018年3月27日
@author: Zhu Qianshu
'''
import sys,getopt
import cPickle
import numpy as np


bed_file = "40000_mm10.bed"
input_file = "test.matrix"
tad_file = "morula2.filt.tad"
output_file = "RTI.test.bg"
cutoff = -1
resolution = 40000

def calc_consolidation(a, b, contacts,cutoff): #contacts: {(start, end):value,...}
    if a > b: a,b = b,a
    D=b-a  #D is distance of a to b, start from 0
    intraTAD,outTAD = [0.0,0.0],[0.0,0.0]
    for x in range(a-int((D+1)/2),b+1):
        for y in range(a,b+int((D+1)/2)+2):
            if x+cutoff <= y and y <= x+D and y+x <= 2*b+1 and y+x >= 2*a-1 and x > 0 and y > 0 :   #cutoff must be interge, would be 0
                if x>=a and y<=b:
                    intraTAD[1] +=1
                    try:
                        intraTAD[0] += contacts[(x,y)]
                    except:
                        pass  #intraTAD[0] +=0
                else:
                    outTAD[1] +=1
                    try:
                        outTAD[0] += contacts[(x,y)]
                    except:
                        pass
    
    try:
        score = float((intraTAD[0]/intraTAD[1])/(outTAD[0]/outTAD[1]))
    except:
        score = 0.0
    return intraTAD[0],intraTAD[1], outTAD[0], outTAD[1],score

def calc_RTI(a, b, contacts):
    if a > b: a,b = b,a
    D=b-a
    inTAD = []  #inTAD[0] are all of first line interactions in TAD: I1_TAD 
    outTAD = [] #outTAD[0] are all of first line interactions out of TAD: I1_out_TAD 
    RTIs = []
    #print a,b,D
    for i in range(0,D+1):
        inTAD.append([])
        outTAD.append([])
    for x in range(a-int((D+1)/2),b+1):
        for y in range(a,b+int((D+1)/2)+2):
            if x <= y and y <= x+D and y+x <= 2*b+1 and y+x >= 2*a-1 and x > 0 and y > 0 :  # y is bigger than or equal to x, y is less than x+D
                if x>=a and y<=b:
                    try:
                        inTAD[y-x].append(contacts[(x,y)])
                        #print x,y,contacts[(x,y)],"in"
                    except:
                        inTAD[y-x].append(0.0) #if contacts has no this interaction, that means it's 0
                else:
                    try:
                        outTAD[y-x].append(contacts[(x,y)])
                        #print x,y,contacts[(x,y)],"out"
                    except:
                        outTAD[y-x].append(0.0)

    for i in range(0,D+1):  # from 0 to b-a
        #print "inTAD",inTAD[i]
        #print "outTAD",outTAD[i]
        in_median = 0.0 if len(inTAD[i]) == 0 else np.median(inTAD[i])
        out_median = 0.0 if len(outTAD[i]) == 0 else np.median(outTAD[i])
        
        if in_median != 0.0 and out_median != 0.0:         #RTI = median( median(I(n)_TAD)/median(I(n)_outTAD) )
            RTIs.append(in_median/out_median)
        elif in_median != 0.0 and out_median == 0.0:
            RTIs.append(float("inf"))
        elif in_median == 0.0 and out_median != 0.0:
            RTIs.append(0.0)
        elif in_median == 0.0 and out_median == 0.0:
            pass
    #print "RTIs",RTIs
    return 0.0 if len(RTIs) == 0 else np.median(RTIs)

def calc_condense(a, b, contacts,resolution):
    #mat should be raw contacts
    if a > b: a,b = b,a
    ll=[] #length list
    for i in range(a,b+1):
        for j in range(i,b+1):
            try:
                contact = int(contacts[(i,j)])
            except:
                contact = 0
            ll += [int((j-i+0.5) * resolution)]*contact
    return (0.0,0.0) if len(ll) == 0 else (np.mean(ll),np.median(ll))
            
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

def main_consolidation(tad, step, mat, resolution,cutoff,output_file):
    out = open(output_file,'w')
    out.write("chr\tstart\tend\tin_score\tin_num\tout_score\tout_num\tratio\n")
    for chrom in tad.keys():
        for (tadS,tadE) in tad[chrom]:
            tadA = int(tadS/resolution) + step[chrom][0]
            tadB = int(tadE/resolution) + step[chrom][0]
            a,b,c,d,score = calc_consolidation(tadA, tadB, mat[chrom],cutoff/resolution)
            out.write("%s\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n"%(chrom,tadS,tadE,a,b,c,d,score))
    
def main_RTI(tad, step, mat, resolution,output_file):
    out = open(output_file,'w')
    out.write("chr\tstart\tend\tRTI\n")
    for chrom in tad.keys():
        for (tadS,tadE) in tad[chrom]:
            tadA = int(tadS/resolution) + step[chrom][0]
            tadB = int(tadE/resolution) + step[chrom][0]
            score = calc_RTI(tadA, tadB, mat[chrom])
            out.write("%s\t%d\t%d\t%.3f\n"%(chrom,tadS,tadE,score))

def main_condense(tad, step, mat, resolution,output_file):
    out = open(output_file,'w')
    out.write("chr\tstart\tend\tmean\tmedian\n")
    for chrom in tad.keys():
        for (tadS,tadE) in tad[chrom]:
            tadA = int(tadS/resolution) + step[chrom][0]
            tadB = int(tadE/resolution) + step[chrom][0]
            mean,median = calc_condense(tadA, tadB, mat[chrom],resolution)
            out.write("%s\t%d\t%d\t%.2f\t%.2f\n"%(chrom,tadS,tadE,mean,median))


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
            -c cutoff of short distance region(for consolidation) or -1 (for Relative TAD intensity). [default = 400k]
    """
    try:
        opts,args = getopt.getopt(sys.argv[1:], "hi:o:I:c:b:r:", 
                                  ["help"])
    except getopt.GetoptError:
        sys.exit(usage)
    
    for opt in opts:
        if opt[0] == '-h': sys.exit(usage)
        elif opt[0] == '-i': input_file = opt[1]
        elif opt[0] == '-o': output_file = opt[1]
        elif opt[0] == '-I': tad_file = opt[1]
        elif opt[0] == '-b': bed_file = opt[1]
        elif opt[0] == "-c": 
            try:
                cutoff = int(opt[1])
            except:
                cutoff = opt[1]
        elif opt[0] == "-r": resolution = int(opt[1])
        
    tad = parse_tad(tad_file)
    step = parse_bed(bed_file)
    mat = parse_matrix(input_file, step)
#     cPickle.dump(mat, open("test.dump","w"), protocol=0)
#     mat = cPickle.load(open("test.dump",'r'))
    if cutoff == "condense":
        main_condense(tad, step, mat, resolution, output_file)
    elif cutoff == "RTI" or cutoff < 0 :
        main_RTI(tad, step, mat, resolution, output_file)
    else:
        main_consolidation(tad, step, mat, resolution, cutoff, output_file)


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        