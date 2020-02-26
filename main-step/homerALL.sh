#######20180414 all_except_norm_depth_20180412:
wdir=~/workspace/8.NT-HiC/h.homer_ALL/1.all_except_norm_depth_20180412
ddir=~/workspace/8.NT-HiC/d.HiC-Pro_analysis/1.merge_samples_except_20180411/2.allValidPairs

bash ~/codes/hic_homer_threads.sh -w $wdir -d $ddir -r1 25000 -r2 100000 -t 4

for i in *bedGraph; do sed -i 's|/home/qszhu/workspace/8.NT-HiC/h.homer_ALL/1.all_except_norm_depth_20180412/a.PC1/||' $i;done
#########20180414 xiewei homer all  norm depth:
wdir=~/workspace/8.NT-HiC/h.homer_ALL/1.all_except_norm_depth_20180412
ddir=~/workspace/8.NT-HiC/d.HiC-Pro_analysis/1.merge_samples_except_20180411/2.allValidPairs

bash ~/codes/hic_homer_threads.sh -w $wdir -d $ddir -r1 25000 -r2 100000 -t 4


#################compartment change analysis
awk '{$1=$1;print}' a2b.bed | sortBed -i - | mergeBed -i - -d 1| intersectBed -a ~/ann/mm10.RefSeq.bed -b - -r -e -f 1 -wa -u |cut -f 4 | sort | uniq > a2b.genes
awk '{$1=$1;print}' b2a.bed | sortBed -i - | mergeBed -i - -d 1| intersectBed -a ~/ann/mm10.RefSeq.bed -b - -r -e -f 1 -wa -u |cut -f 4 | sort | uniq > b2a.genes
wc -l *genes
  752 a2b.genes
  791 b2a.genes
  

#######################################2018年5月13日 all_except 6h_rep1 6h_rep2 morula_rep2
wdir=~/workspace/8.NT-HiC/h.homer_ALL/2.all_except_20180513
ddir=~/workspace/8.NT-HiC/2.merge_sample/4.except/1.allValidPairs
bash ~/codes/hic_homer_threads.sh -w $wdir -d $ddir -r1 25000 -r2 100000 -t 4

for i in *bedGraph; do sed -i 's|/home/qszhu/workspace/8.NT-HiC/h.homer_ALL/2.all_except_20180513/2.PC1/||' $i;done

#2018年5月16日 更换2h 6h的数据
wdir=~/workspace/8.NT-HiC/h.homer_ALL/2.all_except_20180513/3.temp
ln -s $ddir/2h_allValidPairs $wdir
ln -s $ddir/6h_allValidPairs $wdir
bash ~/codes/hic_homer_threads.sh -w $wdir -d $wdir -r1 25000 -r2 100000 -t 4


######2018年6月09日 两次添加的数据较多，因此全部重跑
wdir=~/workspace/8.NT-HiC/h.homer_ALL/3.except_100k_20180609
ddir=~/workspace/8.NT-HiC/2.merge_sample/4.except/1.allValidPairs
bash ~/codes/hic_homer_threads.sh -w $wdir -d $ddir -r1 100000 -r2 400000 -t 8
for i in *bedGraph; do sed -i 's|/home/qszhu/workspace/8.NT-HiC/h.homer_ALL/3.except_100k_20180609/2.PC1/||' $i;done
#分辨率低看起来不对劲
wdir=~/workspace/8.NT-HiC/h.homer_ALL/4.except_25k_20180609
ddir=~/workspace/8.NT-HiC/2.merge_sample/4.except/1.allValidPairs
bash ~/codes/hic_homer_threads.sh -w $wdir -d $ddir -r1 25000 -r2 100000 -t 4
for i in *bedGraph; do sed -i 's|/home/qszhu/workspace/8.NT-HiC/h.homer_ALL/4.except_25k_20180609/2.PC1/||' $i;done

######2018年6月09日 每个重复单独跑
wdir=~/workspace/8.NT-HiC/h.homer_ALL/5.rep_25k_20180610
ddir=~/workspace/8.NT-HiC/2.merge_sample/2.replicate/1.allValidPairs
bash ~/codes/hic_homer_threads.sh -w $wdir -d $ddir -r1 25000 -r2 100000 -t 4 &


##########2018年6月19日 定义compartment的紊乱程度
ddir=~/workspace/8.NT-HiC/h.homer_ALL/4.except_25k_20180609/2.PC1
wdir=~/workspace/8.NT-HiC/h.homer_ALL/4.except_25k_20180609/3.length_compartment
mkdir -p $wdir
for i in $ddir/*bedGraph; do o=$wdir/$(basename $i .PC1.bedGraph).txt
awk 'NR>1{if($4<0){$4="B"}else{$4="A"};print}' $i | \
awk 'NR==1{a=$1;b=$4;n=1} NR>1{if($1==a&&$4==b){n+=1}else{print a,b,n;a=$1;b=$4;n=1} }'  > $o &
done 

for i in ${keys[@]}; do 
echo $i `wc -l $i.txt`
done

#reps
ddir=~/workspace/8.NT-HiC/h.homer_ALL/5.rep_25k_20180610/2.PC1
wdir=~/workspace/8.NT-HiC/h.homer_ALL/5.rep_25k_20180610/3.length_compartment

#############2018年6月20日 compartment间的interaction比例
ddir=~/workspace/8.NT-HiC/5.maps/4.except/1.raw
bdir=~/workspace/8.NT-HiC/h.homer_ALL/4.except_25k_20180609/2.PC1
wdir=~/workspace/8.NT-HiC/h.homer_ALL/4.except_25k_20180609/4.interaction_between_compartment

python ~/codes/interactions_between_compartments.py -i $ddir/morula_25000.matrix -b $bdir/morula.PC1.bedGraph -o $wdir/morula.txt -r 25000


#####################################################################
#####2018年6月29日 更换te 1h icm
ddir=~/workspace/8.NT-HiC/2.merge_sample/4.except/1.allValidPairs
wdir=~/workspace/8.NT-HiC/h.homer_ALL/4.except_25k_20180609/0.add_temp
mkdir $wdir
for i in te 1h icm; do
ln -s $ddir/${i}_allValidPairs $wdir
done
bash ~/codes/hic_homer_threads.sh -w $wdir -d $wdir -r1 25000 -r2 100000 -t 4

wdir=~/workspace/8.NT-HiC/h.homer_ALL/6.rep_25k_20180629
ddir=~/workspace/8.NT-HiC/2.merge_sample/2.replicate/1.allValidPairs
bash ~/codes/hic_homer_threads.sh -w $wdir -d $ddir -r1 25000 -r2 100000 -t 4 &

#add 100k 400k 
wdir=~/workspace/8.NT-HiC/h.homer_ALL/6.rep_25k_20180629
ddir=~/workspace/8.NT-HiC/2.merge_sample/2.replicate/1.allValidPairs
mkdir -p $wdir
bash ~/codes/hic_homer_threads2.sh -w $wdir -d $ddir -r1 100000 -r2 400000 -t 4 -o 3.PC1_100-400k &


########################################################################
#Homer 更新新版本
wdir=~/workspace/8.NT-HiC/h.homer_ALL/7.except_newHomer_20180705
ddir=~/workspace/8.NT-HiC/2.merge_sample/4.except/1.allValidPairs
bash ~/codes/hic_makeTagDirectory.sh -w $wdir -d $ddir -t 8 &

wdir=~/workspace/8.NT-HiC/h.homer_ALL/7.except_newHomer_20180705
odir=2.PC1_50-200k #4.PC1_50-200k #3.TAD
mkdir -p $wdir/$odir/logs; mkfifo $wdir/$odir/logs/temp.fifo; exec 1000<>$wdir/$odir/logs/temp.fifo
for((i=0;i<8;i++));do echo >&1000; done
for i in $wdir/1.tags/*tag; do read -u 1000 
{
k=$(basename $i .tag)
nohup runHiCpca.pl $wdir/$odir/$k $i/ -res 50000 -window 200000 -genome mm10 -cpu 4 > $wdir/$odir/logs/$k.log
#nohup findTADsAndLoops.pl find $i/ -o $wdir/$odir/$k -res 40000 -window 200000 -minDist 200000 -minDist 5000000 -cpu 4 > $wdir/$odir/logs/$k.log
echo >&1000 
echo $i" done"
}& done &


odir=4.PC1_50-200k ##5.TAD  #    #2.PC1_50-100k
mkdir -p $wdir/$odir/logs
for i in $wdir/1.tags/*tag; do 
k=$(basename $i .tag)
nohup runHiCpca.pl $wdir/$odir/$k $i/ -res 50000 -window 200000 -genome mm10 -cpu 4 > $wdir/$odir/logs/$k.log &
#nohup findTADsAndLoops.pl find $i/ -o $wdir/$odir/$k -cpu 8 -genome mm10  > $wdir/$odir/logs/$k.log &
done

for i in *bedGraph; do 
bedGraphToBigWig $i ~/ann/mm10.chromSizes bw/${i/bedGraph/bw} &
done

#compartment length
ddir=~/workspace/8.NT-HiC/h.homer_ALL/7.except_newHomer_20180705/4.PC1_50-200k
wdir=~/workspace/8.NT-HiC/h.homer_ALL/7.except_newHomer_20180705/3.length_compartment
mkdir -p $wdir
for i in $ddir/*bedGraph; do o=$wdir/$(basename $i .PC1.bedGraph).txt
awk 'NR>1{if($4<0){$4="B"}else{$4="A"};print}' $i | \
awk 'NR==1{a=$1;b=$4;n=1} NR>1{if($1==a&&$4==b){n+=1}else{print a,b,n;a=$1;b=$4;n=1} }'  > $o &
done 

for i in ${keys[@]}; do 
echo $i `wc -l $i.txt`
done

######2018年7月29日 icm seq2
rm -r 1.tags/icm.tag/ 0.input/* 4.PC1_50-200k/icm.PC1.* 
wdir=~/workspace/8.NT-HiC/h.homer_ALL/7.except_newHomer_20180705
ln -s ~/workspace/8.NT-HiC/2.merge_sample/4.except/1.allValidPairs/icm_allValidPairs $wdir/0.input/
bash ~/codes/hic_makeTagDirectory.sh -w $wdir -d $wdir/0.input -t 8 &
nohup runHiCpca.pl $wdir/4.PC1_50-200k/icm $wdir/1.tags/icm.tag -res 50000 -window 200000 -genome mm10 -cpu 4 > $wdir/4.PC1_50-200k/logs/icm.log &

###更改05h 1h 2h 
rm -r 1.tags/[012]* 0.input/* 4.PC1_50-200k/[012]* 
wdir=~/workspace/8.NT-HiC/h.homer_ALL/7.except_newHomer_20180705
ln -s ~/workspace/8.NT-HiC/2.merge_sample/4.except/1.allValidPairs/[012]*_allValidPairs $wdir/0.input/
bash ~/codes/hic_makeTagDirectory.sh -w $wdir -d $wdir/0.input -t 8 &
for i in $wdir/0.input/*; do k=$(basename $i _allValidPairs)
nohup runHiCpca.pl $wdir/4.PC1_50-200k/$k $wdir/1.tags/$k.tag -res 50000 -window 200000 -genome mm10 -cpu 4 > $wdir/4.PC1_50-200k/logs/$k.log &
done

##更改12h 增加12h_rep3
wdir=~/workspace/8.NT-HiC/h.homer_ALL/7.except_newHomer_20180705
rm -r $wdir/1.tags/12* $wdir/0.input/* $wdir/4.PC1_50-200k/12* 
ln -s ~/workspace/8.NT-HiC/2.merge_sample/4.except/1.allValidPairs/12*_allValidPairs $wdir/0.input/
bash ~/codes/hic_makeTagDirectory.sh -w $wdir -d $wdir/0.input -t 8 &
for i in $wdir/0.input/*; do k=$(basename $i _allValidPairs)
nohup runHiCpca.pl $wdir/4.PC1_50-200k/$k $wdir/1.tags/$k.tag -res 50000 -window 200000 -genome mm10 -cpu 4 > $wdir/4.PC1_50-200k/logs/$k.log &
done
for i in *bedGraph; do sed -i 's|/home/qszhu/workspace/8.NT-HiC/h.homer_ALL/7.except_newHomer_20180705/4.PC1_50-200k/||' $i;done

#100-400kb compartment
for i in $wdir/1.tags/*; do k=$(basename $i .tag)
nohup runHiCpca.pl $wdir/5.PC1_100-400k/$k $wdir/1.tags/$k.tag -res 100000 -window 400000 -genome mm10 -cpu 8 > $wdir/5.PC1_100-400k/$k.log &
done

#40-200kb compartment
wdir=~/workspace/8.NT-HiC/h.homer_ALL/7.except_newHomer_20180705
for i in $wdir/1.tags/*; do k=$(basename $i .tag)
nohup runHiCpca.pl $wdir/9.PC1_40-200k/$k $wdir/1.tags/$k.tag -res 40000 -window 200000 -genome mm10 -cpu 8 > $wdir/9.PC1_40-200k/$k.log &
done

#300-600kb compartment
wdir=~/workspace/8.NT-HiC/h.homer_ALL/7.except_newHomer_20180705
for i in $wdir/1.tags/*; do k=$(basename $i .tag)
nohup runHiCpca.pl $wdir/2.PC1_300-600k/$k $wdir/1.tags/$k.tag -res 300000 -window 600000 -genome mm10 -cpu 16 > $wdir/2.PC1_300-600k/$k.log &
done

for i in $wdir/9.PC1_40-200k/*.PC1.txt;do k=$(basename $i .PC1.txt)
(findHiCCompartments.pl $i -thresh 0.25 |awk '{print $2,$3,$4,"A"}' > $wdir/6.findHiCCompartments/$k.txt
findHiCCompartments.pl $i -thresh 0.25 -opp |awk '{print $2,$3,$4,"B"}' >> $wdir/6.findHiCCompartments/$k.txt) &
done

######################################
#analyzHiC
wdir=~/workspace/8.NT-HiC/h.homer_ALL/7.except_newHomer_20180705
#make resolution
bash ~/codes/hic_homer_threads.sh -w $wdir -d ~/workspace/8.NT-HiC/2.merge_sample/4.except/1.allValidPairs/ -t 8 &
#observer/expected
#creat bg
res=600000
for i in $wdir/1.tags/*tag; do o=$(basename $i .tag)
nohup analyzeHiC $wdir/1.tags/$o -res $res -bgonly -cpu 16 >> $wdir/logs/${o/.tag/.log} &
done

for i in $wdir/1.tags/*tag; do o=$(basename $i .tag)
#analyzeHiC $i -pos chr19:35,000,000-45,000,000 -res 40000 -distNorm -cpu 6 > $wdir/7.distNorm/$o.chr19.35-45M.txt &
#analyzeHiC $i -chr chr19 -res 100000 -distNorm -cpu 6 > $wdir/7.distNorm/$o.chr19.all.txt &
Rscript ~/R/8.NTHiC/plot_pheatmap_homer_nofilt_0.R $wdir/7.distNorm/$o.chr19.35-45M.txt &
Rscript ~/R/8.NTHiC/plot_pheatmap_homer_nofilt_0.R $wdir/7.distNorm/$o.chr19.all.txt &
done

for i in $wdir/1.tags/*tag; do o=$wdir/8.corr/$(basename $i .tag)
(analyzeHiC $i -chr chr19 -res 100000 -window 400000 -corr -distNorm -cpu 16 -o $o.chr19.all.txt 
Rscript ~/R/8.NTHiC/plot_pheatmap_homer_nofilt_0.R $o.chr19.all.txt 
)&
done

for i in $wdir/1.tags/*tag; do o=$wdir/b.corr-300-600/$(basename $i .tag)
(analyzeHiC $i -chr chr19 -res 300000 -window 600000 -corr -distNorm -cpu 16 -o $o.chr19.all.txt 
Rscript ~/R/8.NTHiC/plot_pheatmap_homer_nofilt_0.R $o.chr19.all.txt 
analyzeHiC $i -chr chr12 -start 0 -end 121000000 -res 300000 -window 600000 -corr -distNorm -cpu 16 -o $o.chr12.txt 
Rscript ~/R/8.NTHiC/plot_pheatmap_homer_nofilt_0.R $o.chr12.txt 
)&
done


for i in $wdir/1.tags/*tag; do o=$wdir/a.dlr/$(basename $i .tag)
analyzeHiC $i -chr chr12 -start 0 -end 121000000 -res 300000 -window 600000 -corr -distNorm -cpu 16 -o $o.all.txt  &
done

########################################
#########reps
########################################
wdir=~/workspace/8.NT-HiC/h.homer_ALL/8.rep_newHomer_20180705
ddir=~/workspace/8.NT-HiC/2.merge_sample/2.replicate/1.allValidPairs
bash ~/codes/hic_makeTagDirectory.sh -w $wdir -d $ddir -t 8 &

wdir=~/workspace/8.NT-HiC/h.homer_ALL/8.rep_newHomer_20180705
odir=4.PC1_50-200k #3.TAD    #2.PC1_50-100k
mkdir -p $wdir/$odir/logs; mkfifo $wdir/$odir/logs/temp.fifo; exec 1000<>$wdir/$odir/logs/temp.fifo
for((i=0;i<8;i++));do echo >&1000; done
for i in $wdir/1.tags/*tag; do read -u 1000 
{
k=$(basename $i .tag)
nohup runHiCpca.pl $wdir/$odir/$k $i/ -res 50000 -window 200000 -genome mm10 -cpu 4 > $wdir/$odir/logs/$k.log
#nohup findTADsAndLoops.pl find $i/ -o $wdir/$odir/$k -res 40000 -window 200000 -minDist 200000 -minDist 5000000 -cpu 4 > $wdir/$odir/logs/$k.log
echo >&1000 
echo $i" done"
}& done &
for i in *bedGraph; do sed -i 's|/home/qszhu/workspace/8.NT-HiC/e.ES/3.homer/2.PC1_50-200k/xw_rep1 ||' $i;done

##增加12h_rep3
wdir=~/workspace/8.NT-HiC/h.homer_ALL/8.rep_newHomer_20180705
rm $wdir/0.input/*
ln -s ~/workspace/8.NT-HiC/2.merge_sample/2.replicate/1.allValidPairs/12h_rep3*_allValidPairs $wdir/0.input/
bash ~/codes/hic_makeTagDirectory.sh -w $wdir -d $wdir/0.input -t 8 &
for i in $wdir/0.input/*; do k=$(basename $i _allValidPairs)
nohup runHiCpca.pl $wdir/4.PC1_50-200k/$k $wdir/1.tags/$k.tag -res 50000 -window 200000 -genome mm10 -cpu 4 > $wdir/4.PC1_50-200k/logs/$k.log &
done
for i in *bedGraph; do sed -i 's|/home/qszhu/workspace/8.NT-HiC/h.homer_ALL/8.rep_newHomer_20180705/2.PC1_50-100k||' $i;done


wdir=~/workspace/8.NT-HiC/h.homer_ALL/8.rep_newHomer_20180705
for res in 100000 400000 40000 200000; do
for i in $wdir/1.tags/*tag; do o=$(basename $i .tag)
nohup analyzeHiC $wdir/1.tags/$o.tag -res $res -bgonly -cpu 16 >> $wdir/logs/${o/.tag/.log}
done &
done

#40-200kb compartment
wdir=~/workspace/8.NT-HiC/h.homer_ALL/8.rep_newHomer_20180705
mkdir -p $wdir/9.PC1_40-200k
for i in $wdir/1.tags/*; do k=$(basename $i .tag)
nohup runHiCpca.pl $wdir/9.PC1_40-200k/$k $wdir/1.tags/$k.tag -res 40000 -window 200000 -genome mm10 -cpu 8 > $wdir/9.PC1_40-200k/$k.log &
done

#100-400kb compartment
mkdir -p $wdir/5.PC1_100-400k
for i in $wdir/1.tags/*; do k=$(basename $i .tag)
nohup runHiCpca.pl $wdir/5.PC1_100-400k/$k $wdir/1.tags/$k.tag -res 100000 -window 400000 -genome mm10 -cpu 2 > $wdir/5.PC1_100-400k/$k.log &
done

################################
#correlation of correlation matrix
###############################
wdir=~/workspace/8.NT-HiC/h.homer_ALL/7.except_newHomer_20180705
for i in $wdir/1.tags/*tag; do o=$wdir/c.correlation_of_correlation/$(basename $i .tag)
for j in $(seq 1 19) X; do 
nohup analyzeHiC $i -res 300000 -window 600000 -chr chr$j -corr -distNorm -cpu 16 -o $o.chr$j.txt > $o.chr$j.log
done &
done

wdir=~/workspace/8.NT-HiC/h.homer_ALL/7.except_newHomer_20180705/c.correlation_of_correlation
for i in ${NT[@]}; do i=${i%%:*}
awk 'FNR==1{a=4} FNR>1{for(i=a;i<=NF;i++){print $i};a++}' $i.chr*txt >> $i.temp
echo $i
done
#check
paste `for i in ${NT[@]}; do i=${i%%:*};echo -n $i".temp ";done` > all.tab

#做与NF-ICM的相关性
wdir=~/workspace/8.NT-HiC/h.homer_ALL/7.except_newHomer_20180705/c.correlation_of_correlation
i=~/workspace/8.NT-HiC/3.public_data/f.XW_analysis/7.newHomer/1.tags/ICM.tag;
o=$wdir/$(basename $i .tag)
for j in $(seq 1 19) X; do 
nohup analyzeHiC $i -res 300000 -window 600000 -chr chr$j -corr -distNorm -cpu 6 -o $o.chr$j.txt > $o.chr$j.log &
done
wait
awk 'FNR==1{a=4} FNR>1{for(i=a;i<=NF;i++){print $i};a++}' ICM.chr*txt >> ICM.temp
paste `for i in ${NT[@]} ICM; do i=${i%%:*};echo -n $i".temp ";done` > all-ICM.tab