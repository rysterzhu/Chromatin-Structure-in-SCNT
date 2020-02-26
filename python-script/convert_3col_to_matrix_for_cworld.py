#--coding:utf-8
'''
Created on 2017年12月8日
@author: Zhu Qianshu
'''
import sys,getopt

OFS="\t"
ORS="\n"


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
        opts,args = getopt.getopt(sys.argv[1:], "hi:o:I:c:", 
                                  ["help"])
    except getopt.GetoptError:
        sys.exit(usage)
    
    for opt in opts:
        if opt[0] == '-h': sys.exit(usage)
        elif opt[0] == '-i': input_file = opt[1]
        elif opt[0] == '-o': output_file = opt[1]
        elif opt[0] == '-I': bed_file = opt[1]
        elif opt[0] == "-c": chr = opt[1]

###test
# input_file = "nt2h.matrix"
# bed_file = "40000_mm10.bed"
# output_file = "nt2h.mat"

###        

regions = {}
firstFlag = True
for line in open(bed_file,"rU"):
    fs = line.rstrip().rsplit("\t")
    if not fs[0] == chr:
        continue
    k = int(fs[3])
    if firstFlag:
        MIN = k
        MAX = k
        firstFlag = False
    if MIN > k: MIN = k
    if MAX < k: MAX = k
    regions[k] = fs[3] + "|mm10|" + fs[0] + ":" + str(int(fs[1])+1) + "-" + fs[2]
print MIN,MAX
melt = {}        
for line in open(input_file,"rU"):
    fs = line.rstrip().rsplit("\t")
    #if int(fs[0]) < MIN or int(fs[1]) < MIN or int(fs[0]) > MAX or int(fs[1]) > MAX: continue
    if MIN <= int(fs[0]) <= MAX and MIN <= int(fs[1]) <= MAX:
        melt[fs[0] + "-" + fs[1]] = fs[2]




out = open(output_file,"w")
out.write(OFS)
for i in range(MIN,MAX):
    out.write(regions[i] + OFS)
out.write(regions[MAX] + ORS)

for i in range(MIN,MAX+1):
    out.write(regions[i] + OFS)
    for j in range(MIN,MAX):
        if i < j:
            k = str(i) + "-" + str(j)
        else:
            k = str(j) + "-" + str(i)
        if melt.has_key(k):
            out.write(melt[k] + OFS)
        else:
            out.write("nan" + OFS)
    k = str(i) + "-" + str(MAX)
    if melt.has_key(k):
        out.write(melt[k] + ORS)
    else:
        out.write("nan" + ORS)
        
out.close()




















