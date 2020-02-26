wdir=~/workspace/8.NT-HiC/p.loops/b.FIND/1.e2cell-sample
ddir=~/workspace/8.NT-HiC/5.maps/4.except/1.raw
ddir2=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/3.raw

awk 'NR==FNR&&$1=="chr7"&&$2==10000000{a=$4} NR==FNR&&$1=="chr7"&&$2==14000000{b=$4} 
NR>FNR&&$1>=a&&$2<=b{$1=$1-a+1;$2=$2-a+1;print}' $ConfigHP/100000_mm10.bed $ddir/e2cell_100000.matrix > $wdir/NT.matrix &
awk 'NR==FNR&&$1=="chr7"&&$2==10000000{a=$4} NR==FNR&&$1=="chr7"&&$2==14000000{b=$4} 
NR>FNR&&$1>=a&&$2<=b{$1=$1-a+1;$2=$2-a+1;print}' $ConfigHP/100000_mm10.bed $ddir2/early_2cell_100000.matrix > $wdir/NF.matrix &
awk 'NR==FNR&&$1=="chr7"&&$2==10000000{a=$4} NR==FNR&&$1=="chr7"&&$2==14000000{b=$4} 
NR>FNR&&$1>=a&&$2<=b{$1=$1-a+1;$2=$2-a+1;print}' $ConfigHP/100000_mm10.bed $ddir/kdm4d-e2cell_100000.matrix > $wdir/kdm4d.matrix &

##########replicates
wdir=~/workspace/8.NT-HiC/p.loops/b.FIND/2.e2cell-replicates
ddir=~/workspace/8.NT-HiC/5.maps/2.replicate/1.raw
ddir2=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/7.replicate/3.raw
mkd $wdir

awk 'NR==FNR&&$1=="chr7"&&$2==10000000{a=$4} NR==FNR&&$1=="chr7"&&$2==14000000{b=$4} 
NR>FNR&&$1>a&&$2<=b{$1=$1-a;$2=$2-a;print;print $2,$1,$3}' $ConfigHP/100000_mm10.bed $ddir2/early_2cell_rep4_100000.matrix > $wdir/NF-rep2.matrix &
awk 'NR==FNR&&$1=="chr7"&&$2==10000000{a=$4} NR==FNR&&$1=="chr7"&&$2==14000000{b=$4} 
NR>FNR&&$1>a&&$2<=b{$1=$1-a;$2=$2-a;print;print $2,$1,$3}' $ConfigHP/100000_mm10.bed $ddir2/early_2cell_rep5_100000.matrix > $wdir/NF-rep3.matrix &
awk 'NR==FNR&&$1=="chr7"&&$2==10000000{a=$4} NR==FNR&&$1=="chr7"&&$2==14000000{b=$4} 
NR>FNR&&$1>a&&$2<=b{$1=$1-a;$2=$2-a;print;print $2,$1,$3}' $ConfigHP/100000_mm10.bed $ddir2/early_2cell_rep6_100000.matrix > $wdir/NF-rep1.matrix &

awk 'NR==FNR&&$1=="chr7"&&$2==10000000{a=$4} NR==FNR&&$1=="chr7"&&$2==14000000{b=$4} 
NR>FNR&&$1>a&&$2<=b{$1=$1-a;$2=$2-a;print;print $2,$1,$3}' $ConfigHP/100000_mm10.bed $ddir/e2cell_rep3_100000.matrix > $wdir/NT-rep1.matrix &
awk 'NR==FNR&&$1=="chr7"&&$2==10000000{a=$4} NR==FNR&&$1=="chr7"&&$2==14000000{b=$4} 
NR>FNR&&$1>a&&$2<=b{$1=$1-a;$2=$2-a;print;print $2,$1,$3}' $ConfigHP/100000_mm10.bed $ddir/e2cell_rep4_100000.matrix > $wdir/NT-rep2.matrix &

awk 'NR==FNR&&$1=="chr7"&&$2==10000000{a=$4} NR==FNR&&$1=="chr7"&&$2==14000000{b=$4} 
NR>FNR&&$1>a&&$2<=b{$1=$1-a;$2=$2-a;print;print $2,$1,$3}' $ConfigHP/100000_mm10.bed $ddir/kdm4d-e2cell_rep1_100000.matrix > $wdir/kdm4d-rep1.matrix &
awk 'NR==FNR&&$1=="chr7"&&$2==10000000{a=$4} NR==FNR&&$1=="chr7"&&$2==14000000{b=$4} 
NR>FNR&&$1>a&&$2<=b{$1=$1-a;$2=$2-a;print;print $2,$1,$3}' $ConfigHP/100000_mm10.bed $ddir/kdm4d-e2cell_rep2_100000.matrix > $wdir/kdm4d-rep2.matrix &


##################################可以直接用HiC-Pro的数据，注意看vignette!
#注意0-100000是第一个bin，位置是1；因此11100000-11200000是第112个bin，位置是112
wdir=~/workspace/8.NT-HiC/p.loops/b.FIND/3.e2cell-chr7-iced
ddir=~/workspace/8.NT-HiC/5.maps/2.replicate/2.iced
ddir2=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/7.replicate/4.iced

#test chr7
a=$(grep chr7 /home/share/HiC-Pro-annotations/100000_mm10.bed | head -1 | cut -f 4)
b=$(grep chr7 /home/share/HiC-Pro-annotations/100000_mm10.bed | tail -n 1 | cut -f 4)
grep chr7 /home/share/HiC-Pro-annotations/100000_mm10.bed > $wdir/100000_chr7.bed
for i in $ddir/e2cell_rep[34]*_100000_iced.matrix $ddir2/early_2cell*_100000_iced.matrix; do o=$(basename $i _100000_iced.matrix)
awk -v a=$a -v b=$b '$1>=a&&$2<=b{print}' $i > $wdir/$o.chr7.matrix &
done


###test 40k in chr 7
wdir=~/workspace/8.NT-HiC/p.loops/b.FIND/4.e2cell-chr7-iced-40k
ddir=~/workspace/8.NT-HiC/5.maps/2.replicate/2.iced
ddir2=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/7.replicate/4.iced
mkdir $wdir
#test chr7
a=$(grep chr7 /home/share/HiC-Pro-annotations/40000_mm10.bed | head -1 | cut -f 4)
b=$(grep chr7 /home/share/HiC-Pro-annotations/40000_mm10.bed | tail -n 1 | cut -f 4)
grep chr7 /home/share/HiC-Pro-annotations/40000_mm10.bed > $wdir/40000_chr7.bed
for i in $ddir/e2cell_rep[34]*_40000_iced.matrix $ddir2/early_2cell*_40000_iced.matrix; do o=$(basename $i _40000_iced.matrix)
awk -v a=$a -v b=$b '$1>=a&&$2<=b{print}' $i > $wdir/$o.chr7.matrix &
done

wdir=~/workspace/8.NT-HiC/p.loops/b.FIND/7.kdm4d-NT-chr7-iced
ddir=~/workspace/8.NT-HiC/5.maps/2.replicate/2.iced
ddir2=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/7.replicate/4.iced
a=$(grep chr7 /home/share/HiC-Pro-annotations/100000_mm10.bed | head -1 | cut -f 4)
b=$(grep chr7 /home/share/HiC-Pro-annotations/100000_mm10.bed | tail -n 1 | cut -f 4)
grep chr7 /home/share/HiC-Pro-annotations/100000_mm10.bed > $wdir/100000_chr7.bed
for i in $ddir/e2cell_rep[34]*_100000_iced.matrix $ddir/kdm4d*_100000_iced.matrix; do o=$(basename $i _100000_iced.matrix)
awk -v a=$a -v b=$b '$1>=a&&$2<=b{print}' $i > $wdir/$o.chr7.matrix &
done


######################################################all chromosome
wdir=~/workspace/8.NT-HiC/p.loops/b.FIND/5.e2cell-allchr-iced
ddir=~/workspace/8.NT-HiC/5.maps/2.replicate/2.iced
ddir2=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/7.replicate/4.iced

grep -v "chrM|chrY" /home/share/HiC-Pro-annotations/100000_mm10.bed > $wdir/100000_mm10.bed
for i in $ddir/e2cell_rep[34]*_100000_iced.matrix $ddir2/early_2cell*_100000_iced.matrix; do o=$(basename $i _100000_iced.matrix)
awk 'NR==FNR&&$2==0{begin[$1]=$4} NR==FNR{end[$1]=$4} 
NR>FNR{for(c in begin){if($1>=begin[c]&&$2<=end[c]){print;next}}}
' 100000_mm10.bed $i > $wdir/$o.intra.matrix &
done #又慢又容易出错

for chr in {1..19} X; do chr="chr"$chr
awk -v chr=$chr '$1==chr{print}' /home/share/HiC-Pro-annotations/100000_mm10.bed > $wdir/0.data/100000_$chr.bed
a=$(head -1 $wdir/0.data/100000_$chr.bed| cut -f 4)
b=$(tail -1 $wdir/0.data/100000_$chr.bed| cut -f 4)
for i in $ddir/e2cell_rep[34]*_100000_iced.matrix $ddir2/early_2cell*_100000_iced.matrix; do o=$(basename $i _100000_iced.matrix)
awk -v a=$a -v b=$b '$1>=a&&$2<=b{print}' $i > $wdir/0.data/$o.$chr.matrix &
done;done

for chr in {1..19} X; do chr="chr"$chr
nohup Rscript $wdir/5.e2cell-allchr-iced.FIND.R $chr > $wdir/logs/$chr.R.log &
done

#####
for chr in {1..19} X; do chr="chr"$chr
awk -v chr=$chr 'NR==FNR{a[$1"\t"$2]=$3} NR>FNR{print chr,$0,a[$1"\t"$2]}' $wdir/4.foldchange/NT-NF.$chr.fc.txt $wdir/2.qvalues/NT-NF.quantreg.$chr.qvalue.txt >> NT-NF.allchr.qvalue-FC.tab 
done
awk '$1=="chr7"&&$2==112&&$3==131' NT-NF.allchr.qvalue-FC.tab #Zscan4d, col4 qvalue, all <0.01; col5 FC,NF/NT;
chr7	112	131	2.55743647326161e-68	4.30747816744012

awk '$5<0.5{a+=1} $5>2{b+=1} END{print a,b}' NT-NF.allchr.qvalue-FC.tab
448136	303819
######Se-P analysis
k27ac=/mnt/guru/home4/qszhu/workspace/9.NT-ChIP/2.public/a.Renbin/5.homer/1.bed/2cell_H3K27ac.super5.bed
wdir=~/workspace/8.NT-HiC/p.loops/b.FIND/5.e2cell-allchr-iced/6.SE-P
promoter=~/ann/mm10_promoter_1000-1000.bed

awk '$5>=0.4&&$5<=2.5{next} $3-$2<=2||$3-$2>30||$4>1e-5{next}
$5>1{f="NF";a+=1} $5<1{f="NT";b+=1}
{print $1,$2*100000-100000,$2*100000,f"-"NR > "left.bin.bed"} {print $1,$3*100000-100000,$3*100000,f"-"NR > "right.bin.bed"} 
END{print a,b}' ../NT-NF.allchr.qvalue-FC.tab


#左到右，右到左；分别得到两边与enhancer和promoter的交集，enhancer要以强度排序
intersectBed -a left.bin.bed -b $k27ac -r -e -f 0.5 -wo | cut -f 1,2,3,4,8,10 | sort -k6g,6 > left.enhancer  
intersectBed -a right.bin.bed -b $promoter -r -e -f 0.5 -wo | cut -f 1,2,3,4,8,9 > right.promoter
awk 'NR==FNR{left[$4]=$0} NR>FNR{if($4 in left){print $0,left[$4]}}' left.enhancer right.promoter > left2right.ep

intersectBed -a right.bin.bed -b $k27ac -r -e -f 0.5 -wo | cut -f 1,2,3,4,8,10 | sort -k6g,6 > right.enhancer  
intersectBed -a left.bin.bed -b $promoter -r -e -f 0.5 -wo | cut -f 1,2,3,4,8,9 > left.promoter
awk 'NR==FNR{left[$4]=$0} NR>FNR{if($4 in left){print $0,left[$4]}}' right.enhancer left.promoter > right2left.ep

cat left2right.ep right2left.ep >all.ep
awk '$4~"NF"{nf[$6]} $4~"NT"{nt[$6]} END{for(i in nf){if(!(i in nt))print i,"NF"};for(i in nt){if(!(i in nf))print i,"NT"}}' all.ep > all.gene.tab

cut -f 2 all.gene.tab | sort | uniq -c


##################################################################################################################chr15 500k
wdir=~/workspace/8.NT-HiC/p.loops/b.FIND/6.chr15-iced-500k
ddir=~/workspace/8.NT-HiC/5.maps/2.replicate/2.iced
ddir2=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/7.replicate/4.iced
mkd $wdir/0.data $wdir/1.RData $wdir/2.qvalues $wdir/3.pdfs $wdir/logs
chr=chr15
awk -v chr=$chr '$1==chr{print}' /home/share/HiC-Pro-annotations/500000_mm10.bed > $wdir/0.data/500000_$chr.bed
a=$(head -1 $wdir/0.data/500000_$chr.bed| cut -f 4)
b=$(tail -1 $wdir/0.data/500000_$chr.bed| cut -f 4)

for i in $ddir/e2cell_rep[34]*_500000_iced.matrix $ddir/icm*_500000_iced.matrix $ddir/l2cell*_500000_iced.matrix $ddir/8cell*_500000_iced.matrix ; do o=$(basename $i _500000_iced.matrix)
awk -v a=$a -v b=$b '$1>=a&&$2<=b{print}' $i > $wdir/0.data/NT-$o.matrix &
done

for i in $ddir2/early*_500000_iced.matrix $ddir2/ICM*_500000_iced.matrix $ddir2/late*_500000_iced.matrix $ddir2/8cell*_500000_iced.matrix ; do o=$(basename $i _500000_iced.matrix)
awk -v a=$a -v b=$b '$1>=a&&$2<=b{print}' $i > $wdir/0.data/NF-$o.matrix &
done




