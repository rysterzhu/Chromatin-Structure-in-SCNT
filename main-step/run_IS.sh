wdir=~/workspace/8.NT-HiC/f.IS_ALL/1.all_except_norm_depth_20180412
ddir=~/workspace/8.NT-HiC/d.HiC-Pro_analysis/1.merge_samples_except_20180411/4.norm_to_depth 
bash ~/codes/hic_insulation_score.sh -w $wdir -d $ddir
keys=("cc" "05h" "1h" "2h" "6h" "e2cell" "l2cell" "4cell" "8cell" "morula")
for i in ${keys[@]}; do 
echo -en $i"\t"
cat $i.bed | wc -l
done
#for insulation score plot
for i in ${keys[@]}; do 
awk -v i=${i} '{print i,$5}' $i.bed >> all.ins
done


####################
wdir=~/workspace/8.NT-HiC/f.IS_ALL/2.norm_test
ddir=~/workspace/8.NT-HiC/d.HiC-Pro_analysis/1.merge_samples_except_20180411/7.norm_test
bash ~/codes/hic_insulation_score.sh -w $wdir -d $ddir
keys=("cc" "05h" "1h" "2h" "6h" "e2cell" "l2cell" "4cell" "8cell" "morula")
for i in ${keys[@]}; do 
echo -en $i"\t"
cat $i.bed | wc -l
done
#for insulation score plot
for i in ${keys[@]}; do 
awk -v i=${i} '{print i,$5}' $i.bed >> all.ins
done



#########################################################################plot 
wdir=~/workspace/8.NT-HiC/f.IS_ALL/1.all_except_norm_depth_20180412
ddir=~/workspace/8.NT-HiC/f.IS_ALL/1.all_except_norm_depth_20180412/cat_insulation


for k in ${keys[@]}; do 
#echo "track type=bedGraph name='$k--is120000--nt0--ids80000--ss40000--imiqrMean.insulation'" > $wdir/$k.bedGraph
(bedSort $wdir/cat_insulation/$k.bedGraph $wdir/cat_insulation/$k.bedGraph
bedGraphToBigWig $wdir/cat_insulation/$k.bedGraph ~/ann/mm10.chromSizes $wdir/cat_insulation/$k.bw )&
done


(computeMatrix reference-point -S $ddir/cc.bw $ddir/05h.bw $ddir/1h.bw \
-R $wdir/cat_boundary/filted/cc.filted.bed -p 16 -a 1400000 -b 1400000 --referencePoint center -bs 40000 --skipZeros -o $wdir/plotProfile/1.gz
plotProfile -m $wdir/plotProfile/1.gz -out $wdir/plotProfile/1.pdf --perGroup --regionsLabel "" --refPointLabel boundary \
-T "Insulation score around CC TAD boundary" --legendLocation upper-right --samplesLabel "Donor Cell" "Transfer 0.5h" "Transfer 1h" \
--colors '#2b8cbe' '#a6bddb' '#ece7f2') &


(computeMatrix reference-point -S $ddir/1h.bw $ddir/2h.bw $ddir/6h.bw $ddir/e2cell.bw $ddir/l2cell.bw $ddir/4cell.bw $ddir/8cell.bw $ddir/morula.bw \
-R $wdir/cat_boundary/filted/morula.filted.bed -p 16 -a 1400000 -b 1400000 --referencePoint center -bs 40000 --skipZeros -o $wdir/plotProfile/2.gz
plotProfile -m $wdir/plotProfile/2.gz -out $wdir/plotProfile/2.pdf --perGroup --regionsLabel "" --refPointLabel boundary \
-T "Insulation score around morula TAD boundary" --legendLocation upper-right \
--samplesLabel "Activate" "Activated 1h" "Activated 6h" "Early 2 Cell" "Late 2 Cell" "4 Cell" "8 Cell" "Morula" \
--colors '#fff7fb' '#ece7f2' '#d0d1e6' '#a6bddb' '#74a9cf' '#3690c0' '#0570b0' '#034e7b') &



awk '$3-$2>0{a=($3-$2)/2;$2=$2-a;$3=$3+a;print $0}' ~/workspace/8.NT-HiC/g.DI_ALL/1.all_except_norm_depth_20180412/unoverlap_tads/morula.tad > morula.2tad.bed
awk '$3-$2>0{a=($3-$2)/2;$2=$2-a;$3=$3+a;print $0}' ~/workspace/8.NT-HiC/g.DI_ALL/1.all_except_norm_depth_20180412/unoverlap_tads/cc.tad > cc.2tad.bed

wdir=~/workspace/8.NT-HiC/f.IS_ALL/1.all_except_norm_depth_20180412/plotProfile
ddir=~/workspace/8.NT-HiC/f.IS_ALL/1.all_except_norm_depth_20180412/cat_insulation
key=is_tad1
color=('#2b8cbe' '#a6bddb' '#ece7f2')
(computeMatrix scale-regions -S $ddir/cc.bw $ddir/05h.bw $ddir/1h.bw \
-R $wdir/cc.2tad.bed -p 16 -m 2000000 -a 1000000 -b 1000000 -bs 40000 --skipZeros -o $wdir/$key.gz
plotProfile -m $wdir/$key.gz -out $wdir/$key.pdf --perGroup --regionsLabel "" --startLabel "start" --endLabel "end" -T "Insulation Score around CC TAD(DI)" \
--legendLocation upper-right --colors ${color[@]} --samplesLabel "Donor Cell" "Transfer 0.5h" "Transfer 1h" ) &

key=is_tad2
color=('#fff7fb' '#ece7f2' '#d0d1e6' '#a6bddb' '#74a9cf' '#3690c0' '#0570b0' '#034e7b')
(computeMatrix scale-regions -S $ddir/1h.bw $ddir/2h.bw $ddir/6h.bw $ddir/e2cell.bw $ddir/l2cell.bw $ddir/4cell.bw $ddir/8cell.bw $ddir/morula.bw \
-R $wdir/morula.2tad.bed -p 16 -m 2000000 -a 1000000 -b 1000000 -bs 40000 --skipZeros -o $wdir/$key.gz
plotProfile -m $wdir/$key.gz -out $wdir/$key.pdf --perGroup --regionsLabel "" --startLabel "start" --endLabel "end" -T "Insulation Score around morula TAD(DI)" \
--legendLocation upper-right --colors ${color[@]} --samplesLabel "Activate" "Activated 1h" "Activated 6h" "Early 2 Cell" "Late 2 Cell" "4 Cell" "8 Cell" "Morula" ) &

##########################################################
##########2018年5月14日 + 3h 12h + 16日 + 2h_rep3 -6h_rep4
##########################################################
wdir=~/workspace/8.NT-HiC/f.IS_ALL/3.all_except_norm_mean_20180514
ddir=~/workspace/8.NT-HiC/5.maps/4.except/6.norm_to_mean
bash ~/codes/hic_insulation_score.sh -w $wdir -d $ddir &

for i in ${keys[@]}; do 
echo -en $i"\t"
cat $i.bed | wc -l
done
#for insulation score plot
for i in ${keys[@]}; do 
awk -v i=${i} '{print i,$5}' $i.bed >> all.ins
done

echo -e "sample\tnumber" > counts.tab
for i in *bed;do k=${i%%.*}; c=`cat $i|wc -l `; echo -e "$k\t$c" >> counts.tab; done

######plotProfile
for k in ${keys[@]}; do 
#echo "track type=bedGraph name='$k--is120000--nt0--ids80000--ss40000--imiqrMean.insulation'" > $wdir/$k.bedGraph
(bedSort $wdir/cat_insulation/$k.bedGraph $wdir/cat_insulation/$k.bedGraph
bedGraphToBigWig $wdir/cat_insulation/$k.bedGraph ~/ann/mm10.chromSizes $wdir/cat_insulation/$k.bw )&
done

wdir=~/workspace/8.NT-HiC/f.IS_ALL/3.all_except_norm_mean_20180514
ddir=~/workspace/8.NT-HiC/f.IS_ALL/3.all_except_norm_mean_20180514/cat_insulation

(computeMatrix reference-point -S $ddir/cc.bw $ddir/05h.bw $ddir/1h.bw \
-R $wdir/cat_boundary/filted/cc.filted.bed -p 16 -a 1400000 -b 1400000 --referencePoint center -bs 40000 --skipZeros -o $wdir/plotProfile/1.gz
plotProfile -m $wdir/plotProfile/1.gz -out $wdir/plotProfile/1.pdf --perGroup --regionsLabel "" --refPointLabel boundary \
-T "Insulation score around CC TAD boundary" --legendLocation upper-right --samplesLabel "Donor Cell" "Transfer 0.5h" "Transfer 1h" \
--colors '#2b8cbe' '#a6bddb' '#ece7f2') &

(computeMatrix reference-point -S $ddir/1h.bw $ddir/2h.bw $ddir/3h.bw $ddir/6h.bw $ddir/12h.bw $ddir/e2cell.bw $ddir/l2cell.bw $ddir/4cell.bw $ddir/8cell.bw $ddir/morula.bw \
-R $wdir/cat_boundary/filted/morula.filted.bed -p 16 -a 1400000 -b 1400000 --referencePoint center -bs 40000 --skipZeros -o $wdir/plotProfile/3.gz
plotProfile -m $wdir/plotProfile/3.gz -out $wdir/plotProfile/3.pdf --perGroup --regionsLabel "" --refPointLabel boundary \
-T "Insulation score around morula TAD boundary" --legendLocation upper-right \
--samplesLabel "Activate 0h" "Activated 1h" "Activated 3h" "Activated 6h" "Activated 12h" "Early 2 Cell" "Late 2 Cell" "4 Cell" "8 Cell" "Morula" \
--colors '#a6cee3' '#1f78b4' '#b2df8a' '#33a02c' '#fb9a99' '#e31a1c' '#fdbf6f' '#ff7f00' '#cab2d6' '#6a3d9a') &

(computeMatrix reference-point -S $ddir/cc.bw $ddir/2h.bw $ddir/3h.bw $ddir/6h.bw \
-R $wdir/cat_boundary/filted/cc.filted.bed -p 16 -a 1400000 -b 1400000 --referencePoint center -bs 40000 --skipZeros -o $wdir/plotProfile/4.gz
plotProfile -m $wdir/plotProfile/4.gz -out $wdir/plotProfile/4.pdf --perGroup --regionsLabel "" --refPointLabel boundary \
-T "Insulation score around cc TAD boundary" --legendLocation upper-right \
--samplesLabel "Cumulus Cell" "Activated 1h" "Activated 3h" "Activated 6h" \
--colors '#e41a1c' '#377eb8' '#4daf4a' '#984ea3') &

taddir=~/workspace/8.NT-HiC/g.DI_ALL/6.except_norm_mean_20180521/2.unoverlap_tads
(computeMatrix scale-regions -S $ddir/cc.bw $ddir/2h.bw $ddir/3h.bw $ddir/6h.bw \
-R $taddir/cc.tad -p 16 -p 16 -m 1200000 -a 600000 -b 600000 -bs 40000 --skipZeros -o $wdir/plotProfile/5.gz
plotProfile -m $wdir/plotProfile/5.gz -out $wdir/plotProfile/5.pdf --perGroup --regionsLabel "" --startLabel "start" --endLabel "end" \
-T "Insulation score around cc TAD" --legendLocation upper-right \
--samplesLabel "Cumulus Cell" "Activated 1h" "Activated 3h" "Activated 6h" \
--colors '#e41a1c' '#377eb8' '#4daf4a' '#984ea3') &

(computeMatrix reference-point -S $ddir/12h.bw $ddir/e2cell.bw $ddir/l2cell.bw $ddir/4cell.bw $ddir/8cell.bw $ddir/morula.bw \
-R $wdir/cat_boundary/filted/morula.filted.bed -p 16 -a 1400000 -b 1400000 --referencePoint center -bs 40000 --skipZeros -o $wdir/plotProfile/6.gz
plotProfile -m $wdir/plotProfile/6.gz -out $wdir/plotProfile/6.pdf --perGroup --regionsLabel "" --refPointLabel boundary \
-T "Insulation score around morula TAD boundary" --legendLocation upper-right \
--samplesLabel "Activated 12h" "Early 2 Cell" "Late 2 Cell" "4 Cell" "8 Cell" "Morula" \
--colors '#fb9a99' '#e31a1c' '#fdbf6f' '#ff7f00' '#cab2d6' '#6a3d9a') &

taddir=~/workspace/8.NT-HiC/g.DI_ALL/6.except_norm_mean_20180521/2.unoverlap_tads
(computeMatrix scale-regions -S $ddir/12h.bw $ddir/e2cell.bw $ddir/l2cell.bw $ddir/4cell.bw $ddir/8cell.bw $ddir/morula.bw \
-R $taddir/morula.tad -p 16 -p 16 -m 1200000 -a 600000 -b 600000 -bs 40000 --skipZeros -o $wdir/plotProfile/7.gz
plotProfile -m $wdir/plotProfile/7.gz -out $wdir/plotProfile/7.pdf --perGroup --regionsLabel "" --startLabel "start" --endLabel "end" \
-T "Insulation score around morula TAD" --legendLocation upper-right \
--samplesLabel "Activated 12h" "Early 2 Cell" "Late 2 Cell" "4 Cell" "8 Cell" "Morula" \
--colors '#fb9a99' '#e31a1c' '#fdbf6f' '#ff7f00' '#cab2d6' '#6a3d9a') &

(computeMatrix reference-point -S $ddir/cc.bw $ddir/6h.bw $ddir/12h.bw \
-R $wdir/cat_boundary/filted/cc.filted.bed -p 16 -a 1400000 -b 1400000 --referencePoint center -bs 40000 --skipZeros -o $wdir/plotProfile/8.gz
plotProfile -m $wdir/plotProfile/8.gz -out $wdir/plotProfile/8.pdf --perGroup --regionsLabel "" --refPointLabel boundary \
-T "Insulation score around cc TAD boundary" --legendLocation upper-right \
--samplesLabel "Cumulus Cell" "Activated 6h" "Activated 12h" \
--colors '#e41a1c' '#377eb8' '#4daf4a') &

taddir=~/workspace/8.NT-HiC/g.DI_ALL/6.except_norm_mean_20180521/2.unoverlap_tads
(computeMatrix scale-regions -S $ddir/cc.bw $ddir/6h.bw $ddir/12h.bw \
-R $taddir/cc.tad -p 16 -p 16 -m 1200000 -a 600000 -b 600000 -bs 40000 --skipZeros -o $wdir/plotProfile/9.gz
plotProfile -m $wdir/plotProfile/9.gz -out $wdir/plotProfile/9.pdf --perGroup --regionsLabel "" --startLabel "start" --endLabel "end" \
-T "Insulation score around cc TAD" --legendLocation upper-right \
--samplesLabel "Cumulus Cell" "Activated 6h" "Activated 12h" \
--colors '#e41a1c' '#377eb8' '#4daf4a') &


##########################################################
##########2018年5月31日 + 12h_rep1_seq2 1h_rep3_seq1 2h_rep3_seq2 icm_rep1_seq1 te_rep1_seq1
########2018年6月29日 + te_rep1_seq2 1h_rep3_seq2 icm_rep2_seq1 te_rep3_seq1
#####################2018年7月4日 12h_rep2加测  05h_rep3  修改去除的rep
##########################################################
wdir=~/workspace/8.NT-HiC/f.IS_ALL/4.except_norm_mean_raw_20180531
ddir=~/workspace/8.NT-HiC/5.maps/4.except/5.norm_to_mean
bash ~/codes/hic_insulation_score.sh -w $wdir -d $ddir -r 40000 &

for i in ${keys[@]}; do 
echo -en $i"\t"
cat $i.bed | wc -l
done
#for insulation score plot
for i in ${keys[@]}; do 
awk -v i=${i} '{print i,$5}' $i.bed >> all.ins
done

echo -e "sample\tnumber" > counts.tab
for i in *bed;do k=${i%%.*}; c=`cat $i|wc -l `; echo -e "$k\t$c" >> counts.tab; done

######plotProfile
for k in ${keys[@]}; do 
#echo "track type=bedGraph name='$k--is120000--nt0--ids80000--ss40000--imiqrMean.insulation'" > $wdir/$k.bedGraph
(bedSort $wdir/cat_insulation/$k.bedGraph $wdir/cat_insulation/$k.bedGraph
bedGraphToBigWig $wdir/cat_insulation/$k.bedGraph ~/ann/mm10.chromSizes $wdir/cat_insulation/$k.bw )&
done

ddir=~/workspace/8.NT-HiC/f.IS_ALL/4.except_norm_mean_raw_20180531/cat_insulation
taddir=~/workspace/8.NT-HiC/g.DI_ALL/10.except_norm_1e7_20180727/2.unoverlap_tads
(computeMatrix scale-regions -S $ddir/cc.bw $ddir/6h.bw $ddir/12h.bw \
-R $taddir/icm.tad -p 32 -m 1200000 -a 600000 -b 600000 -bs 40000 --skipZeros -o $wdir/plotProfile/9.gz
plotProfile -m $wdir/plotProfile/9.gz -out $wdir/plotProfile/9.pdf --perGroup --regionsLabel "" --startLabel "start" --endLabel "end" \
-T "Insulation score around icm TAD" --legendLocation upper-right \
--samplesLabel "Cumulus Cell" "Activated 6h" "Activated 12h" \
--colors '#e41a1c' '#377eb8' '#4daf4a') &
####################################################################################


wdir=~/workspace/8.NT-HiC/f.IS_ALL/3.all_except_norm_mean_20180514
ddir=~/workspace/8.NT-HiC/f.IS_ALL/3.all_except_norm_mean_20180514/cat_insulation

#
wdir=~/workspace/8.NT-HiC/f.IS_ALL/6.except_norm_mean_raw_20180713
ddir=~/workspace/8.NT-HiC/5.maps/4.except/5.norm_to_mean
bash ~/codes/hic_insulation_score.sh -w $wdir -d $ddir -r 40000 &

######plotProfile
for k in ${keys[@]}; do 
#echo "track type=bedGraph name='$k--is120000--nt0--ids80000--ss40000--imiqrMean.insulation'" > $wdir/$k.bedGraph
(bedSort $wdir/cat_insulation/$k.bedGraph $wdir/cat_insulation/$k.bedGraph
bedGraphToBigWig $wdir/cat_insulation/$k.bedGraph ~/ann/mm10.chromSizes $wdir/cat_insulation/$k.bw )&
done

ddir=~/workspace/8.NT-HiC/f.IS_ALL/6.except_norm_mean_raw_20180713/cat_insulation
taddir=~/workspace/8.NT-HiC/g.DI_ALL/10.except_norm_1e7_20180727/2.unoverlap_tads
(computeMatrix scale-regions -S $ddir/cc.bw $ddir/6h.bw $ddir/12h.bw \
-R $taddir/icm.tad -p 16 -m 1200000 -a 600000 -b 600000 -bs 40000 --skipZeros -o $wdir/plotProfile/9.gz
plotProfile -m $wdir/plotProfile/9.gz -out $wdir/plotProfile/9.pdf --perGroup --regionsLabel "" --startLabel "start" --endLabel "end" \
-T "Insulation score around icm TAD" --legendLocation upper-right \
--samplesLabel "Cumulus Cell" "Activated 6h" "Activated 12h" \
--colors '#e41a1c' '#377eb8' '#4daf4a') &


#######################################2018年8月2日 norm to 100M; norm 没有变化
wdir=~/workspace/8.NT-HiC/f.IS_ALL/8.except_norm_mean_raw_20180802
ddir=~/workspace/8.NT-HiC/5.maps/4.except/5.norm_to_mean
bash ~/codes/hic_insulation_score.sh -w $wdir -d $ddir -r 40000 &

######plotProfile
for k in ${keys[@]}; do 
#echo "track type=bedGraph name='$k--is120000--nt0--ids80000--ss40000--imiqrMean.insulation'" > $wdir/$k.bedGraph
(bedSort $wdir/cat_insulation/$k.bedGraph $wdir/cat_insulation/$k.bedGraph
bedGraphToBigWig $wdir/cat_insulation/$k.bedGraph ~/ann/mm10.chromSizes $wdir/cat_insulation/$k.bw )&
done

ddir=~/workspace/8.NT-HiC/f.IS_ALL/8.except_norm_mean_raw_20180802/cat_insulation
taddir=~/workspace/8.NT-HiC/g.DI_ALL/10.except_norm_1e7_20180727/2.unoverlap_tads
(computeMatrix scale-regions -S $ddir/cc.bw $ddir/6h.bw $ddir/12h.bw \
-R $taddir/icm.tad -p 16 -m 1200000 -a 600000 -b 600000 -bs 40000 --skipZeros -o $wdir/plotProfile/9.gz
plotProfile -m $wdir/plotProfile/9.gz -out $wdir/plotProfile/9.pdf --perGroup --regionsLabel "" --startLabel "start" --endLabel "end" \
-T "Insulation score around icm TAD" --legendLocation upper-right \
--samplesLabel "Cumulus Cell" "Activated 6h" "Activated 12h" \
--colors '#e41a1c' '#377eb8' '#4daf4a') &


##########################################################
########分reps
##########################################################
wdir=~/workspace/8.NT-HiC/f.IS_ALL/5.allRep_norm_mean_20180630
ddir=~/workspace/8.NT-HiC/5.maps/2.replicate/5.norm_to_mean
bash ~/codes/hic_insulation_score.sh -w $wdir -d $ddir -r 100000

######plotProfile
for i in $wdir/cat_insulation/*bedGraph; do k=$(basename $i .bedGraph)
#echo "track type=bedGraph name='$k--is120000--nt0--ids80000--ss40000--imiqrMean.insulation'" > $wdir/$k.bedGraph
(bedSort $i $wdir/plotProfile/0.bw/$k.bg
bedGraphToBigWig $wdir/plotProfile/0.bw/$k.bg ~/ann/mm10.chromSizes $wdir/plotProfile/0.bw/$k.bw )&
done
#部分没有的rep，报错

wdir=~/workspace/8.NT-HiC/f.IS_ALL/5.allRep_norm_mean_20180630
ddir=~/workspace/8.NT-HiC/f.IS_ALL/5.allRep_norm_mean_20180630/plotProfile/0.bw/
taddir=~/workspace/8.NT-HiC/g.DI_ALL/10.except_norm_1e7_20180727/2.unoverlap_tads
(computeMatrix scale-regions -S $ddir/cc_rep4.bw $ddir/6h_rep3.bw $ddir/6h_rep5.bw $ddir/12h_rep1.bw $ddir/12h_rep2.bw \
-R $taddir/icm.tad -p 16 -m 1200000 -a 600000 -b 600000 -bs 40000 --skipZeros -o $wdir/plotProfile/9.gz
plotProfile -m $wdir/plotProfile/9.gz -out $wdir/plotProfile/9.pdf --perGroup --regionsLabel "" --startLabel "start" --endLabel "end" \
-T "Insulation score around icm TAD" --legendLocation upper-right \
--samplesLabel "Cumulus Cell" "Activated 6h1" "Activated 6h2" "Activated 12h1" "Activated 12h2" \
--colors '#fb9a99' '#e31a1c' '#fdbf6f' '#ff7f00' '#cab2d6') &

##########################################################
########2018年8月1日 norm_to_100M 与 NF对应
##########################################################
wdir=~/workspace/8.NT-HiC/f.IS_ALL/7.except_norm_1e8_20180801
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
bash ~/codes/hic_insulation_score.sh -w $wdir -d $ddir -r 40000 &

######plotProfile
for k in ${keys[@]}; do 
#echo "track type=bedGraph name='$k--is120000--nt0--ids80000--ss40000--imiqrMean.insulation'" > $wdir/$k.bedGraph
(bedSort $wdir/cat_insulation/$k.bedGraph $wdir/cat_insulation/$k.bedGraph
bedGraphToBigWig $wdir/cat_insulation/$k.bedGraph ~/ann/mm10.chromSizes $wdir/cat_insulation/$k.bw )&
done

ddir=~/workspace/8.NT-HiC/f.IS_ALL/7.except_norm_1e8_20180801/cat_insulation
taddir=~/workspace/8.NT-HiC/g.DI_ALL/10.except_norm_1e7_20180727/2.unoverlap_tads
(computeMatrix scale-regions -S $ddir/cc.bw $ddir/6h.bw $ddir/12h.bw \
-R $taddir/icm.tad -p 16 -m 1200000 -a 600000 -b 600000 -bs 40000 --skipZeros -o $wdir/plotProfile/9.gz
plotProfile -m $wdir/plotProfile/9.gz -out $wdir/plotProfile/9.pdf --perGroup --regionsLabel "" --startLabel "start" --endLabel "end" \
-T "Insulation score around icm TAD" --legendLocation upper-right \
--samplesLabel "Cumulus Cell" "Activated 6h" "Activated 12h" \
--colors '#e41a1c' '#377eb8' '#4daf4a') &


#####################验证12h 的测序
wdir=~/workspace/8.NT-HiC/f.IS_ALL/b.test_12h
ln -s ~/workspace/8.NT-HiC/5.maps/1.sequence/2.iced/12h_rep*_40000_iced.matrix ./
for i in *.matrix; do k=$(basename $i _40000_iced.matrix)
awk -v OFS="\t" -v key=$k '{depth+=$3} END{print key,depth}' $i >> $wdir/depth.raw.tab &
done
sort -k1,1 $wdir/depth.raw.tab | awk '{total+=$2;n+=1;print} END{print "total",total;printf("%s\t%d\n","mean",total/n)}' > $wdir/depth.sorted.tab
for i in $wdir/*_40000_iced.matrix; do key=${i##*/};key=${key/_40000_iced.matrix}
depth=`awk -v key=$key '{if(key==$1)print $2}' $wdir/depth.sorted.tab`
awk -v f=1e8 -v depth=$depth '{print $1,$2,$3*f/depth}' $i > $wdir/temp/${i##*/} &
echo "normalize: "$i
done;
mv temp/* ./

bash ~/codes/hic_insulation_score.sh -w $wdir -d $wdir -r 40000 &

for i in $wdir/cat_insulation/*bedGraph; do k=$(basename $i .bedGraph)
#echo "track type=bedGraph name='$k--is120000--nt0--ids80000--ss40000--imiqrMean.insulation'" > $wdir/$k.bedGraph
(bedSort $i $wdir/plotProfile/0.bw/$k.bg
bedGraphToBigWig $wdir/plotProfile/0.bw/$k.bg ~/ann/mm10.chromSizes $wdir/plotProfile/0.bw/$k.bw )&
done

ddir=$wdir/plotProfile/0.bw/
taddir=~/workspace/8.NT-HiC/g.DI_ALL/10.except_norm_1e7_20180727/2.unoverlap_tads
(computeMatrix scale-regions -S $ddir/12h_rep1.bw $ddir/12h_rep1_seq1.bw $ddir/12h_rep1_seq2.bw \
-R $taddir/icm.tad -p 16 -m 1200000 -a 600000 -b 600000 -bs 40000 --skipZeros -o $wdir/plotProfile/9.gz
plotProfile -m $wdir/plotProfile/9.gz -out $wdir/plotProfile/9.pdf --perGroup --regionsLabel "" --startLabel "start" --endLabel "end" \
-T "Insulation score around icm TAD" --legendLocation upper-right \
) &


##########################################
###2018年8月3日 修改参数
wdir=~/workspace/8.NT-HiC/f.IS_ALL/9.except_norm_1e8_20180803
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
bash ~/codes/hic_insulation_score.sh -w $wdir -d $ddir -r 40000 &

taddir=~/workspace/8.NT-HiC/g.DI_ALL/10.except_norm_1e7_20180727/2.unoverlap_tads
n=1
(computeMatrix scale-regions -S $wdir/plotProfile/0.bw/cc.bw $wdir/plotProfile/0.bw/6h.bw $wdir/plotProfile/0.bw/12h.bw \
-R $taddir/icm.tad -p 16 -m 1200000 -a 600000 -b 600000 -bs 40000 --skipZeros -o $wdir/plotProfile/$n.gz
plotProfile -m $wdir/plotProfile/$n.gz -out $wdir/plotProfile/$n.png --perGroup --regionsLabel "" --startLabel "start" --endLabel "end" \
-T "Insulation score around icm TAD" --legendLocation upper-right \
--samplesLabel "Cumulus Cell" "Activated 6h" "Activated 12h" \
--colors '#e41a1c' '#377eb8' '#4daf4a') &

n=2
(computeMatrix scale-regions -S $wdir/plotProfile/0.bw/6h.bw $wdir/plotProfile/0.bw/12h.bw $wdir/plotProfile/0.bw/e2cell.bw $wdir/plotProfile/0.bw/l2cell.bw $wdir/plotProfile/0.bw/4cell.bw $wdir/plotProfile/0.bw/8cell.bw $wdir/plotProfile/0.bw/morula.bw $wdir/plotProfile/0.bw/icm.bw $wdir/plotProfile/0.bw/te.bw $wdir/plotProfile/0.bw/random.bw \
-R $taddir/icm.tad -p 16 -m 1200000 -a 600000 -b 600000 -bs 40000 --skipZeros -o $wdir/plotProfile/$n.gz
plotProfile -m $wdir/plotProfile/$n.gz -out $wdir/plotProfile/$n.pdf --perGroup --regionsLabel "" --startLabel "start" --endLabel "end" \
-T "Insulation score around icm TAD" --legendLocation upper-right \
--samplesLabel "Activated 6h" "Activated 12h" "early 2-cell" "late 2-cell" "4-cell" "8-cell" "morula" "ICM" "TE" "random" \
--colors '#a6cee3' '#1f78b4' '#b2df8a' '#33a02c' '#fb9a99' '#e31a1c' '#fdbf6f' '#ff7f00' '#cab2d6' '#6a3d9a') &

n=3
(computeMatrix scale-regions -S $wdir/plotProfile/0.bw/6h.bw $wdir/plotProfile/0.bw/12h.bw $wdir/plotProfile/0.bw/e2cell.bw $wdir/plotProfile/0.bw/l2cell.bw $wdir/plotProfile/0.bw/8cell.bw $wdir/plotProfile/0.bw/icm.bw \
-R $taddir/icm.tad -p 16 -m 1200000 -a 600000 -b 600000 -bs 40000 --skipZeros -o $wdir/plotProfile/$n.gz
plotProfile -m $wdir/plotProfile/$n.gz -out $wdir/plotProfile/$n.pdf --perGroup --regionsLabel "" --startLabel "start" --endLabel "end" \
-T "Insulation score around icm TAD" --legendLocation upper-right \
--samplesLabel "Activated 6h" "Activated 12h" "early 2-cell" "late 2-cell" "8-cell" "ICM" \
--colors '#fb9a99' '#e31a1c' '#fdbf6f' '#ff7f00' '#cab2d6' '#6a3d9a') &

######################################################
#以上8月1号到3号全作废
######################################################
wdir=~/workspace/8.NT-HiC/f.IS_ALL/11.except_norm_100M_20180803
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
bash ~/codes/hic_insulation_score.sh -w $wdir -d $ddir -r 40000 &

cp ~/workspace/8.NT-HiC/g.DI_ALL/d.map_is_to_tad/2.bw/random.bw /home/qszhu/workspace/8.NT-HiC/f.IS_ALL/11.except_norm_100M_20180803/plotProfile/0.bw/
taddir=~/workspace/8.NT-HiC/g.DI_ALL/10.except_norm_1e7_20180727/2.unoverlap_tads
n=2
(computeMatrix scale-regions -S $wdir/plotProfile/0.bw/6h.bw $wdir/plotProfile/0.bw/12h.bw $wdir/plotProfile/0.bw/e2cell.bw $wdir/plotProfile/0.bw/l2cell.bw $wdir/plotProfile/0.bw/4cell.bw $wdir/plotProfile/0.bw/8cell.bw $wdir/plotProfile/0.bw/morula.bw $wdir/plotProfile/0.bw/icm.bw $wdir/plotProfile/0.bw/te.bw $wdir/plotProfile/0.bw/random.bw \
-R $taddir/icm.tad -p 16 -m 1200000 -a 600000 -b 600000 -bs 40000 --skipZeros -o $wdir/plotProfile/$n.gz
plotProfile -m $wdir/plotProfile/$n.gz -out $wdir/plotProfile/$n.pdf --perGroup --regionsLabel "" --startLabel "start" --endLabel "end" \
-T "Insulation score around icm TAD" --legendLocation upper-right \
--samplesLabel "Activated 6h" "Activated 12h" "early 2-cell" "late 2-cell" "4-cell" "8-cell" "morula" "ICM" "TE" "random" \
--colors '#a6cee3' '#1f78b4' '#b2df8a' '#33a02c' '#fb9a99' '#e31a1c' '#fdbf6f' '#ff7f00' '#cab2d6' '#6a3d9a') &

n=3
(computeMatrix scale-regions -S $wdir/plotProfile/0.bw/6h.bw $wdir/plotProfile/0.bw/12h.bw $wdir/plotProfile/0.bw/e2cell.bw $wdir/plotProfile/0.bw/l2cell.bw $wdir/plotProfile/0.bw/8cell.bw $wdir/plotProfile/0.bw/icm.bw \
-R $taddir/icm.tad -p 16 -m 1200000 -a 600000 -b 600000 -bs 40000 --skipZeros -o $wdir/plotProfile/$n.gz
plotProfile -m $wdir/plotProfile/$n.gz -out $wdir/plotProfile/$n.pdf --perGroup --regionsLabel "" --startLabel "start" --endLabel "end" \
-T "Insulation score around icm TAD" --legendLocation upper-right \
--samplesLabel "Activated 6h" "Activated 12h" "early 2-cell" "late 2-cell" "8-cell" "ICM" \
--colors '#fb9a99' '#e31a1c' '#fdbf6f' '#ff7f00' '#cab2d6' '#6a3d9a') &

n=4
(computeMatrix scale-regions -S $wdir/plotProfile/0.bw/cc.bw $wdir/plotProfile/0.bw/6h.bw $wdir/plotProfile/0.bw/12h.bw $wdir/plotProfile/0.bw/e2cell.bw $wdir/plotProfile/0.bw/l2cell.bw \
-R $taddir/cc.tad -p 16 -m 1200000 -a 600000 -b 600000 -bs 40000 --skipZeros -o $wdir/plotProfile/$n.gz
plotProfile -m $wdir/plotProfile/$n.gz -out $wdir/plotProfile/$n.pdf --perGroup --regionsLabel "" --startLabel "start" --endLabel "end" \
-T "Insulation score around CC TAD" --legendLocation upper-right \
--samplesLabel "CC" "Activated 6h" "Activated 12h" "early 2-cell" "late 2-cell" \
--colors '#a6cee3' '#1f78b4' '#b2df8a' '#33a02c' '#fb9a99' '#e31a1c') &

##2018年8月13日 排除1h_rep3 2h_rep3
wdir=~/workspace/8.NT-HiC/f.IS_ALL/12.except_norm_1e8_20180813
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M


######################
#insulation score peaks test
######################
wdir=~/workspace/8.NT-HiC/f.IS_ALL/11.except_norm_100M_20180803
perl -I ~/1.HiC-software/1.cworld/lib/ ~/1.HiC-software/1.cworld/perl/matrix2insulation.pl --is 480000 --ids 120000 --nt 0.1 --bmoe 3 --im iqrMean -i $wdir/matrix/cc.chr2.mat

#######################
#XieWei parameter
#######################
wdir=~/workspace/8.NT-HiC/f.IS_ALL/13.except_xw_para_20180816
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
nohup bash ~/codes/hic_insulation_score.sh -w $wdir -d $ddir -r 40000 > logs/13.except_xw_para_20180816.log &

#######################
#parameter tes:
wdir=~/workspace/8.NT-HiC/f.IS_ALL/14.para_test_res40k_is2M_ids400k
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
bash ~/codes/hic_insulation_score.sh -w $wdir -d $ddir -r 40000 &

wdir=~/workspace/8.NT-HiC/f.IS_ALL/15.para_test_res100k_is2M_ids400k
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
bash ~/codes/hic_insulation_score.sh -w $wdir -d $ddir -r 100000 &

wdir=~/workspace/8.NT-HiC/f.IS_ALL/16.para_test_res100k_is1M_ids200k
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
bash ~/codes/hic_insulation_score.sh -w $wdir -d $ddir -r 100000 &





rm -r 1.all_except_norm_depth_20180412/ 2.norm_test/ 3.all_except_norm_mean_20180514/ 4.except_norm_mean_raw_20180531/
###############################
#2018年8月25日 12h_rep3 
###############################
wdir=~/workspace/8.NT-HiC/f.IS_ALL/2.except_res40k_is1M_ids200k_nt01   
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
nohup bash ~/codes/hic_insulation_score.sh -w $wdir -d $ddir -r 40000 > logs/20180827.log &


wdir=~/workspace/8.NT-HiC/f.IS_ALL/3.replicate_res40k_is1M_ids200k_nt01
ddir=~/workspace/8.NT-HiC/5.maps/2.replicate/4.norm_to_100M
bash ~/codes/hic_insulation_score.sh -w $wdir -d $ddir -r 40000


#######################
#para_test
ddir1=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
ddir2=~/workspace/8.NT-HiC/5.maps/2.replicate/4.norm_to_100M
wdir=~/workspace/8.NT-HiC/f.IS_ALL/d.para_test
nohup bash ~/codes/hic_insulation_score.sh -w $wdir/3.except_res100k_is1M_ids200k_nt01 -d $ddir1 -r 100000 > logs/1.log &
nohup bash ~/codes/hic_insulation_score.sh -w $wdir/4.rep_res100k_is1M_ids200k_nt01 -d $ddir2 -r 100000 > logs/2.log &

#bmoe test
wdir=~/workspace/8.NT-HiC/f.IS_ALL/3.replicate_res40k_is1M_ids200k_nt01
wc -l $wdir/1.cat_boundary/icm_rep*

for i in $(seq 0 20000 100000); do
awk -v i=$i '$5>=0.25{print $1,$2-i,$3+i}' $wdir/1.cat_boundary/icm_rep1.bed > icm_rep1.bed
awk -v i=$i '$5>=0.25{print $1,$2-i,$3+i}' $wdir/1.cat_boundary/icm_rep2.bed > icm_rep2.bed
echo -en "$i "`cat icm_rep1.bed | wc -l`" "`cat icm_rep2.bed | wc -l`" "`bedtools intersect -a icm_rep1.bed -b icm_rep2.bed -wa -u | wc -l `" "`bedtools intersect -b icm_rep1.bed -a icm_rep2.bed -wa -u | wc -l `" "
bedtools intersect -a icm_rep1.bed -b icm_rep2.bed -u | wc -l
done


#####################
#ICM icm boundary 交集
#######################
wdir=~/workspace/8.NT-HiC/f.IS_ALL/1.except_res40k_is1M_ids200k_nt025/1.cat_boundary
cut -f 1-3 icm.bed > test/icm.bed
cut -f 1-3 cc.bed > test/cc.bed
awk '$5>0.25{print $1,$2,$3}' ~/workspace/8.NT-HiC/3.public_data/f.XW_analysis/4.IS/1.cat_boundary/ICM.bed > ICM.bed



########################
#2018年12月23日 10M test
################
wdir=~/workspace/8.NT-HiC/f.IS_ALL/4.except_res500k_is10M_ids2000k_nt025
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
mkdir -p $wdir/logs
nohup bash ~/codes/hic_insulation_score.sh -w $wdir -d $ddir -r 500000 > $wdir/logs/20181223.log &

wdir=~/workspace/8.NT-HiC/f.IS_ALL/5.XW_res500k_is10M_ids2000k_nt025
ddir=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/6.norm_to_100M
mkdir -p $wdir/logs
nohup bash ~/codes/hic_insulation_score.sh -w $wdir -d $ddir -r 500000 > $wdir/logs/20181223.log &


###############################
#2019年10月25日 revise +e2cell_rep4 6h_rep6 tsa kdm4d partherno sertoli
###############################
wdir=~/workspace/8.NT-HiC/f.IS_ALL/10.except_res40k_is1M_ids200k_nt01-20191025   #这是最终使用的
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
mkd $wdir/logs
nohup bash ~/codes/hic_insulation_score.sh -w $wdir -d $ddir -r 40000 > $wdir/logs/20191025.log &



