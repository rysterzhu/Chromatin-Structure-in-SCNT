#--coding:utf-8
'''
Created on 2017年12月3日
@author: Zhu Qianshu
'''
import sys,getopt


if __name__ == '__main__':
    ################usage###################
    usage = "Usage: " + sys.argv[0]
    usage += """
    Usage:python 
        <required>:
            -i input matrix file
            -I input bed file
    
        [Options]:
            --outdir_cword: output of each chrom for insulation score with cword
            --outdir_origin: output of each chrom as original format
            --out_intra: output file for intrachromosome 
            --out_inter: output file for interchromosome
            
    
    """
    try:
        opts,args = getopt.getopt(sys.argv[1:], "hi:o:I:c:", 
                                  ["help","outdir_cword=","outdir_origin=","out_intra=","out_inter="])
    except getopt.GetoptError:
        sys.exit(usage)
    
    for opt in opts:
        if opt[0] == '-h': sys.exit(usage)
        elif opt[0] == '-i': input_file = opt[1]
        elif opt[0] == '-I': bed_file = opt[1]
    
    regions = {}
    chroms = []    
    for line in open(bed_file,"rU"):
        fs = line.rstrip().rsplit("\t")
        regions[fs[3]] = fs[0:3]#fs[0] + ":" + fs[1] + "-" + fs[2]
        chroms.append(fs[0])    
    
    for opt in opts:
        if opt[0] == '--outdir_cword': 
            out_cword = {}
            for i in chroms: out_cword[i] = open(opt[1] + "/" + i + ".matrix","w")
        elif opt[0] == '--outdir_origin': 
            out_origin = {}
            for i in chroms: out_origin[i] = open(opt[1] + "/" + i + ".matrix","w")
        elif opt[0] == '--out_intra': out_intra = open(opt[1],"w")
        elif opt[0] == '--out_inter': out_inter = open(opt[1],"w")
    


    
    for line in open(input_file,"rU"):
        fs = line.rstrip().rsplit("\t")
        if regions[fs[0]][0] == regions[fs[1]][0]:
            if locals().has_key("out_intra"):
                out_intra.write(line)
            if locals().has_key("out_origin"):
                out_origin[regions[fs[0]][0]].write(line)
            if locals().has_key("out_cword"):
                a = regions[fs[0]]
                b = regions[fs[1]]
                out_cword[regions[fs[0]][0]].write("\t".join((a[0] + ":" + a[1] + "-" + a[2], 
                                                              b[0] + ":" + b[1] + "-" + b[2],fs[2])) + "\n")
        else:
            if locals().has_key("out_inter"):
                out_inter.write(line)
                
    
    if locals().has_key("out_intra"):
        out_intra.close()
    if locals().has_key("out_origin"):
        for i in chroms: out_origin[i].close()
    if locals().has_key("out_cword"):
        for i in chroms: out_cword[i].close()
    if locals().has_key("out_inter"):
        out_inter.close()                

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        