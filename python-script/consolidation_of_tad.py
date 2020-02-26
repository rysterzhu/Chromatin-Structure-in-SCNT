#--coding:utf-8
'''
Created on 2018年3月27日
@author: Zhu Qianshu
'''
import sys,getopt
import cPickle


bed_file = "40000_mm10.bed"
input_file = "morula_40000_iced.matrix"
tad_file = "morula.filt.tad"
output_file = "morula_40000.tadscore.40k.bg"
cutoff = 400000
resolution = 40000

def calc_consolidation(a, b, contacts,cutoff): #contacts: {(start, end):value,...}
    if a > b: a,b = b,a
    D=b-a+1
    intraTAD,outTAD = [0.0,0.0],[0.0,0.0]
    for x in range(a-int(D/2),b+1):
        for y in range(a,b+int(D/2)+2):
            if x+cutoff <= y and y <= x+b-a and y+x <= 2*b and y+x >= 2*a:
                if x>=a and y<=b:
                    intraTAD[1] +=1
                    try:
                        intraTAD[0] += contacts[(x,y)]
                    except:
                        pass
                else:
                    outTAD[1] +=1
                    try:
                        outTAD[0] += contacts[(x,y)]
                    except:
                        pass
    
    try:
        score = (intraTAD[0]/intraTAD[1])/(outTAD[0]/outTAD[1])
    except:
        score = 0
    return intraTAD[0],intraTAD[1], outTAD[0], outTAD[1],score

def parse_bed(bed_file):
    step = {} #{chr:(first_sit,end_site),...}
    isF = True
    first = 0
    for line in open(bed_file,"rU"):
        fs = line.rstrip().rsplit("\t")
        if 0 == int(fs[1]):
            if not isF:
                step[chr] = (first,int(fs[3]) - 1)
            isF = False    
            chr = fs[0]
            first = int(fs[3])
    
    step[chr] = (first,int(fs[3]) - 1)
    return step
    
def parse_tad(tad_file):
    tad = {} #{chr:[(start,end),...],...}
    for line in open(tad_file,"rU"):
        fs = line.rstrip().rsplit("\t")
        try:
            chr = fs[0]
            start = int(fs[1])
            end = int(fs[2])
        except:
            continue
        tad.setdefault(chr,[]).append((start,end))
        
    return tad
        
def parse_matrix(input_file,step):
    mat = {} # {chr:{(start, end):value,...},
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
                    chr = c
                else:
                    chr = None
                break
        if chr != None:    
            mat.setdefault(chr, {})[(start,end)] = value

    return mat

def main(tad, step, mat, resolution,cutoff,output_file):
    out = open(output_file,'w')
    out.write("chr\tstart\tend\tin_score\tin_num\tout_score\tout_num\tratio\n")
    for chr in tad.keys():
        for (tadS,tadE) in tad[chr]:
            tadA = int(tadS/resolution) + step[chr][0]
            tadB = int(tadE/resolution) + step[chr][0]
            a,b,c,d,score = calc_consolidation(tadA, tadB, mat[chr],cutoff/resolution)
            out.write("%s\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n"%(chr,tadS,tadE,a,b,c,d,score))
    



if __name__ == '__main__':
    ################usage###################
    usage = "Usage: " + sys.argv[0]
    usage += """
    Usage:python TAD_analysis.consolidation_of_tad
        <required>:
            -i HiC-Pro iced matirx.# [site site interaction]
            -I input tad bed file. # bed file
            -b annotation bed of iced matrix; # bed site
            -o output file  # chr	start	end	in_score	in_num	out_score	out_num	ratio
            -c cutoff of short distance region. [default = 400k]
            -r resolution of matrix. [default = 40k]
    
        [Options]:
    
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
        elif opt[0] == "-c": cutoff = int(opt[1])
        elif opt[0] == "-r": resolution = int(opt[1])
        
    tad = parse_tad(tad_file)
    step = parse_bed(bed_file)
    mat = parse_matrix(input_file, step)
#     cPickle.dump(mat, open("test.dump","w"), protocol=0)
#     mat = cPickle.load(open("test.dump",'r'))
#     tad = parse_tad(tad_file)
    
    main(tad, step, mat, resolution, cutoff, output_file)


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        