library(magrittr)
library(preprocessCore)
library(data.table)
library(pheatmap)

ss=c("cc","6h","12h","e2cell","l2cell","8cell","icm")

wdir="~/workspace/8.NT-HiC/4.HiCPlotter/d.20180830_subtract/9.TAD_change"
ddir="~/workspace/8.NT-HiC/4.HiCPlotter/d.20180830_subtract/8.NT_intra/0.data"
config=fread("~/ann/hic-pro/40000_mm10.bed",col.names = c("chr","start","end","pos"),key = "pos")

for(s1 in 5:6){
	s2=s1+1
	tad=fread(paste0(wdir,"/",ss[s2],"-",ss[s1],".tad"),col.names = c("chr","start","end"))
	cdata = fread(paste0(ddir,"/",ss[s2],"-",ss[s1],".matrix"), header=F,col.names=c("x","y","z"),key="x,y")

	temp1 = merge(tad,config[,-3],by=c("chr","start"))
	temp2 = merge(tad,config[,-2],by=c("chr","end"))
	temp3 = merge(temp1,temp2,by=c("chr","start","end"))
	message("read done.")

	fractionTAD = matrix(NA, nrow = dim(temp3)[1], ncol = 3)
	colnames(fractionTAD) = c("increase","decrease","total")
	for(i in 1:dim(temp3)[1]){
	  a=temp3[i,pos.x]
	  b=temp3[i,pos.y]
	  z = cdata[x>=a&y<=b,z]
	  fractionTAD[i,] = c(length(which(z>0)),length(which(z<0)),length(z))
	  message(i," ",appendLF = F)
	}
	temp4 = cbind(temp3,fractionTAD)

	write.table(temp4,paste0(wdir,"/fraction_TAD_",ss[s2],"-",ss[s1],".tab"),sep = "\t",col.names = T,row.names = F,quote = F)
	message(ss[s2],"-",ss[s1]," TAD done.")
}