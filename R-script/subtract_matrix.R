library(magrittr)
library(preprocessCore)
library(data.table)

args=commandArgs(T)
res=args[4]
ddir=args[1] #"~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M"
wdir=args[2] #"~/workspace/8.NT-HiC/4.HiCPlotter/d.20180830_subtract/6.NT_whole_genome"
if(args[3] == "NT"){
  ss=c("cc","6h","12h","e2cell","l2cell","8cell","icm")
}else if(args[3] == "NF"){
  ss=c("PN5_zygote","early_2cell","late_2cell","8cell","ICM")
}else{
  ss=strsplit(args[3],":")[[1]]
}

paste0(ss,collapse = ",") %>% message
cdata1 = fread(paste0(ddir,"/",ss[1],"_",res,"_iced.matrix"), header=F,col.names=c("x","y",ss[1]),key="x,y")
for(i in 2:length(ss)){
	cdata2 = fread(paste0(ddir,"/",ss[i],"_",res,"_iced.matrix"), header=F,col.names=c("x","y",ss[i]),key="x,y")
	cdata = merge(cdata1,cdata2,all=T)
	cdata[is.na(cdata)] = 0.0
	temp <- cdata[,-2:-1] %>% as.matrix %>% 
	 # normalize.quantiles %>% 
	  as.data.table

	temp2 = cdata[,.(x,y,z=temp[,2]-temp[,1])]
	colnames(temp2) %>% paste0(collapse = ",") %>% message
	fwrite(temp2,paste0(wdir,"/",ss[i],"-",ss[i-1],".",res,".matrix"),sep="\t",col.names=F)
	cdata1 = cdata2
}
