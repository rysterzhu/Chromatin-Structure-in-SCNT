library(magrittr)
library(preprocessCore)
library(data.table)
library(pheatmap)

args=commandArgs(T)
ddir=args[1] #"~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M"
wdir=args[2] #"~/workspace/8.NT-HiC/4.HiCPlotter/d.20180830_subtract/6.NT_whole_genome"

if(args[3] == "NT"){
  ss=c("cc","6h","12h","e2cell","l2cell","8cell","icm")
}else if(args[3] == "NF"){
  ss=c("PN5_zygote","early_2cell","late_2cell","8cell","ICM")
}else{
  ss=strsplit(args[3],":")[[1]]
}
temp=strsplit(args[4],":")[[1]]
config=fread("~/ann/hic-pro/40000_mm10.bed",col.names = c("chr","start","end","pos"),key = "pos")
a=config[chr==temp[1]&start==as.numeric(temp[2]),pos]
b=config[chr==temp[1]&start==as.numeric(temp[3]),pos]
maX = 20
height = 50

message("a= ",a)
message("b= ",b)

fox = data.table(x=rep(a:b,each=b-a+1),y=rep(a:b,times=b-a+1),key="x,y")
for(i in 2:length(ss)){
  temp2 = fread(paste0(ddir,"/",ss[i],"-",ss[i-1],".matrix"), header=F,col.names=c("x","y","z"),key="x,y")
  temp2 = temp2[x >= a & y <= b,]
  message(ss[i],"-",ss[i-1],appendLF=F)
  paste0(" min = ",min(temp2[,3]),", max = ",max(temp2[,3])) %>% message
  temp2 = merge(fox,temp2,all=T)
  temp2[is.na(temp2)] = 0.0
  temp2[z > maX,z := maX]
  temp2[z < -maX,z := -maX]
  temp2[x>y | x<y-height,z:=NA]
  
  temp3 = as.matrix(dcast(temp2,formula = x~y,value.var = "z",fill=0))[,-1]

  pheatmap(temp3,cluster_cols = F ,cluster_rows = F,border_color = NA, silent = F,
    color= colorRampPalette(c("deepskyblue","black","yellow"))(101),legend = F,
    #color= colorRampPalette(c("white","red"))(2000),
    #breaks = unique(c(seq(min(test5),0,length.out = 1001),seq(0,max(test5),length.out = 1000))),
    show_rownames = F , show_colnames = F,cellwidth = 5,cellheight = 5,
    filename = paste0(wdir,"/",ss[i],"-",ss[i-1],".pdf"))
  pheatmap(temp3,cluster_cols = F ,cluster_rows = F,border_color = NA, silent = F,
    color= colorRampPalette(c("deepskyblue","black","yellow"))(101),legend = F,
    #color= colorRampPalette(c("white","red"))(2000),
    #breaks = unique(c(seq(min(test5),0,length.out = 1001),seq(0,max(test5),length.out = 1000))),
    show_rownames = F , show_colnames = F,cellwidth = 5,cellheight = 5,
    filename = paste0(wdir,"/",ss[i],"-",ss[i-1],".png"))
}