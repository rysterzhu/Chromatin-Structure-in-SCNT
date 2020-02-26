##################random
awk '{srand($1);a=int(rand()*68148);print a,a+$2-$1,$3}' ~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M/12h_40000_iced.matrix | sort -k1n,1 -k2n,2 -S 5% > random_40000_iced.matrix

for j in $(seq 1 19) X; do chr=chr$j
(python ~/codes/convert_3col_to_matrix_for_cworld.py -i random_40000_iced.matrix -I $ConfigHP/40000_mm10.bed -c $chr -o random.$chr.mat
perl -I ~/1.HiC-software/1.cworld/lib/ ~/1.HiC-software/1.cworld/perl/matrix2insulation.pl --is 500000 --ids 120000 --nt 0 --nt 0.25 --bmoe 3 -i $o.$chr.mat
)&
done
######################################################################

taddir=~/workspace/8.NT-HiC/g.DI_ALL/10.except_norm_1e7_20180727/2.unoverlap_tads
wdir=~/workspace/8.NT-HiC/g.DI_ALL/d.map_is_to_tad
ddir=~/workspace/8.NT-HiC/f.IS_ALL/7.except_norm_1e8_20180801/cat_insulation

ln -s $ddir/* $wdir/1.bedGraph/
for i in $wdir/1.bedGraph/*bedGraph; do 
(bedSort $i $wdir/2.bw/$(basename $i)
bedGraphToBigWig $wdir/2.bw/$(basename $i) ~/ann/mm10.chromSizes $wdir/2.bw/$(basename $i bedGraph)bw)&
done

ddir=~/workspace/8.NT-HiC/g.DI_ALL/d.map_is_to_tad/2.bw
n=1
(computeMatrix scale-regions -S $ddir/cc.bw $ddir/12h.bw $ddir/e2cell.bw $ddir/l2cell.bw $ddir/8cell.bw $ddir/icm.bw $ddir/random.bw \
-R $taddir/icm.tad -p 16 -m 1200000 -a 600000 -b 600000 -bs 40000 --skipZeros -o $wdir/$n.gz
plotProfile -m $wdir/$n.gz -out $wdir/$n.png --perGroup --regionsLabel "" --startLabel "start" --endLabel "end" \
-T "Insulation score around ICM TAD" --legendLocation upper-right \
--samplesLabel "CC" "Activated 12h" "Early 2 Cell" "Late 2 Cell" "8 Cell" "ICM" "random" \
--colors '#66c2a5' '#fc8d62' '#8da0cb' '#e78ac3' '#a6d854' '#ffd92f' '#e5c494') &


awk '{n=$3-$2;print $1,$2-n/1,$3+n/2}' $taddir/icm.tad > $wdir/icm.tad
n=2
(computeMatrix scale-regions -S $ddir/cc.bw $ddir/12h.bw $ddir/e2cell.bw $ddir/l2cell.bw $ddir/8cell.bw $ddir/icm.bw $ddir/random.bw \
-R $wdir/icm.tad -p 16 -m 1200000 -a 600000 -b 600000 -bs 40000 --skipZeros -o $wdir/$n.gz
plotProfile -m $wdir/$n.gz -out $wdir/$n.png --perGroup --regionsLabel "" --startLabel "start" --endLabel "end" \
-T "Insulation score around ICM TAD" --legendLocation upper-right \
--samplesLabel "CC" "Activated 12h" "Early 2 Cell" "Late 2 Cell" "8 Cell" "ICM" "random" \
--colors '#66c2a5' '#fc8d62' '#8da0cb' '#e78ac3' '#a6d854' '#ffd92f' '#e5c494') &

n=3
(computeMatrix scale-regions -S $ddir/cc.bw $ddir/6h.bw $ddir/12h.bw $ddir/e2cell.bw $ddir/l2cell.bw $ddir/8cell.bw $ddir/icm.bw $ddir/random.bw \
-R $wdir/icm.tad -p 16 -m 1200000 -a 600000 -b 600000 -bs 40000 --skipZeros -o $wdir/$n.gz
plotProfile -m $wdir/$n.gz -out $wdir/$n.png --perGroup --regionsLabel "" --startLabel "start" --endLabel "end" \
-T "Insulation score around ICM TAD" --legendLocation upper-right \
--samplesLabel "CC" "Activated 6h" "Activated 12h" "Early 2 Cell" "Late 2 Cell" "8 Cell" "ICM" "random" \
--colors '#66c2a5' '#D70626' '#fc8d62' '#8da0cb' '#e78ac3' '#a6d854' '#ffd92f' '#e5c494') &

n=4
(computeMatrix scale-regions -S $ddir/*.bw \
-R $taddir/icm.tad -p 16 -m 1200000 -a 600000 -b 600000 -bs 40000 --skipZeros -o $wdir/$n.gz
plotProfile -m $wdir/$n.gz -out $wdir/$n.png --perGroup --regionsLabel "" --startLabel "start" --endLabel "end" \
-T "Insulation score around ICM TAD" --legendLocation upper-right) &


##############################################2018年8月20日
###2018年10月23日 change para
taddir=~/workspace/8.NT-HiC/g.DI_ALL/2.except_norm_100M_20180825/2.unoverlap_tads
wdir=~/workspace/8.NT-HiC/g.DI_ALL/d.map_is_to_tad
ddir=~/workspace/8.NT-HiC/f.IS_ALL/2.except_res40k_is1M_ids200k_nt01/2.cat_insulation  #不能用nt0.25的，某些is没有值
#create bigwig
ln -f -s $ddir/* $wdir/1.bedGraph/
cp -f ~/workspace/8.NT-HiC/5.maps/4.except/5.random/1.IS/random.insulation.bedGraph $wdir/1.bedGraph/random.bedGraph
for i in $wdir/1.bedGraph/*bedGraph; do 
(bedSort $i $wdir/2.bw/$(basename $i)
bedGraphToBigWig $wdir/2.bw/$(basename $i) ~/ann/mm10.chromSizes $wdir/2.bw/$(basename $i bedGraph)bw)&
done
#plot
taddir=~/workspace/8.NT-HiC/g.DI_ALL/2.except_norm_100M_20180825/2.unoverlap_tads
wdir=~/workspace/8.NT-HiC/g.DI_ALL/d.map_is_to_tad
ddir=~/workspace/8.NT-HiC/g.DI_ALL/d.map_is_to_tad/2.bw
n=5
(computeMatrix scale-regions -S $ddir/cc.bw $ddir/6h.bw $ddir/12h.bw $ddir/e2cell.bw $ddir/random.bw --smartLabels -R $taddir/cc.tad -p 16 -m 1200000 -a 600000 -b 600000 -bs 40000 --skipZeros -o $wdir/$n.gz
plotProfile -m $wdir/$n.gz -out $wdir/$n.pdf --perGroup --regionsLabel "CC TAD" --startLabel "start" --endLabel "end" \
-T "" --legendLocation upper-right \
--colors '#66c2a5' '#fc8d62' '#e78ac3' '#a6d854' '#e5c494' )&
n=6
(computeMatrix scale-regions -S $ddir/e2cell.bw $ddir/l2cell.bw $ddir/4cell.bw $ddir/8cell.bw $ddir/icm.bw $ddir/random.bw --smartLabels -R $taddir/icm.tad -p 16 -m 1200000 -a 600000 -b 600000 -bs 40000 --skipZeros -o $wdir/$n.gz
plotProfile -m $wdir/$n.gz -out $wdir/$n.pdf --perGroup --regionsLabel "ICM TAD" --startLabel "start" --endLabel "end" \
-T "" --legendLocation upper-right \
--colors '#66c2a5' '#fc8d62' '#8da0cb' '#e78ac3' '#a6d854' '#e5c494' )&

#test for all 
n=7
(computeMatrix scale-regions -S $ddir/cc.bw $ddir/6h.bw $ddir/12h.bw $ddir/e2cell.bw $ddir/random.bw --smartLabels -R $taddir/icm.tad -p 16 -m 1200000 -a 600000 -b 600000 -bs 40000 --skipZeros -o $wdir/$n.gz
plotProfile -m $wdir/$n.gz -out $wdir/$n.pdf --perGroup --regionsLabel "ICM TAD" --startLabel "start" --endLabel "end" \
-T "" --legendLocation upper-right \
--colors '#66c2a5' '#fc8d62' '#e78ac3' '#a6d854' '#e5c494' )&
n=8
(computeMatrix scale-regions -S $ddir/e2cell.bw $ddir/l2cell.bw $ddir/4cell.bw $ddir/8cell.bw $ddir/icm.bw $ddir/random.bw --smartLabels -R $taddir/cc.tad -p 16 -m 1200000 -a 600000 -b 600000 -bs 40000 --skipZeros -o $wdir/$n.gz
plotProfile -m $wdir/$n.gz -out $wdir/$n.pdf --perGroup --regionsLabel "CC TAD" --startLabel "start" --endLabel "end" \
-T "" --legendLocation upper-right \
--colors '#66c2a5' '#fc8d62' '#8da0cb' '#e78ac3' '#a6d854' '#e5c494' )&