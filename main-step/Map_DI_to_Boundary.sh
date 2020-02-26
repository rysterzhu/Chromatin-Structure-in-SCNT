ddir=~/workspace/8.NT-HiC/g.DI_ALL/2.except_norm_100M_20180825/4.directional_index
wdir=~/workspace/8.NT-HiC/f.IS_ALL/e.map_DI_to_boundary
mkdir -p $wdir/1.bw


for i in $ddir/*.bedGraph ~/workspace/8.NT-HiC/5.maps/4.except/5.random/2.DI/random.bedGraph;do k=$(basename $i .bedGraph)
(bedSort $i $wdir/1.bw/$k.bedGraph
bedGraphToBigWig $wdir/1.bw/$k.bedGraph ~/ann/mm10.chromSizes $wdir/1.bw/$k.bw )&
done

#2018年10月25日 改random
#nt 0.25
bdir=~/workspace/8.NT-HiC/f.IS_ALL/1.except_res40k_is1M_ids200k_nt025/1.cat_boundary
ddir=~/workspace/8.NT-HiC/f.IS_ALL/e.map_DI_to_boundary/1.bw
wdir=~/workspace/8.NT-HiC/f.IS_ALL/e.map_DI_to_boundary
n=CC_0.25
(computeMatrix reference-point -S $ddir/cc.bw $ddir/6h.bw $ddir/12h.bw $ddir/e2cell.bw $ddir/random.bw --smartLabels \
-R $bdir/cc.bed -p 16 -a 1400000 -b 1400000 --referencePoint center -bs 40000 --skipZeros -o $wdir/$n.gz
plotProfile -m $wdir/$n.gz -out $wdir/$n.pdf --perGroup --regionsLabel "" --refPointLabel boundary \
-T "Directionality Index around CC TAD boundary" --legendLocation upper-right \
--colors '#66c2a5' '#fc8d62' '#e78ac3' '#a6d854' '#e5c494' ) &


n=ICM_0.25
(computeMatrix reference-point -S $ddir/e2cell.bw $ddir/l2cell.bw $ddir/4cell.bw $ddir/8cell.bw $ddir/icm.bw $ddir/random.bw --smartLabels \
-R $bdir/icm.bed -p 16 -a 1400000 -b 1400000 --referencePoint center -bs 40000 --skipZeros -o $wdir/$n.gz
plotProfile -m $wdir/$n.gz -out $wdir/$n.pdf --perGroup --regionsLabel "" --refPointLabel boundary \
-T "Directionality Index around ICM TAD boundary" --legendLocation upper-right \
--colors '#66c2a5' '#fc8d62' '#8da0cb' '#e78ac3' '#a6d854' '#e5c494' ) &
#ICM有问题，总是比8cell低

n=morula_test
(computeMatrix reference-point -S $ddir/e2cell.bw $ddir/l2cell.bw $ddir/4cell.bw $ddir/8cell.bw $ddir/morula.bw \
-R $bdir/icm.bed -p 16 -a 1400000 -b 1400000 --referencePoint center -bs 40000 --skipZeros -o $wdir/$n.gz
plotProfile -m $wdir/$n.gz -out $wdir/$n.pdf --perGroup --regionsLabel "" --refPointLabel boundary \
-T "Directionality Index around ICM TAD boundary" --legendLocation upper-right \
--colors '#66c2a5' '#fc8d62' '#8da0cb' '#e78ac3' '#a6d854' '#e5c494') &
#morula和8cell差不多

n=te_test
(computeMatrix reference-point -S $ddir/e2cell.bw $ddir/l2cell.bw $ddir/4cell.bw $ddir/8cell.bw $ddir/morula.bw $ddir/icm.bw $ddir/te.bw \
-R $bdir/icm.bed -p 16 -a 1400000 -b 1400000 --referencePoint center -bs 40000 --skipZeros -o $wdir/$n.gz
plotProfile -m $wdir/$n.gz -out $wdir/$n.pdf --perGroup --regionsLabel "" --refPointLabel boundary \
-T "Directionality Index around ICM TAD boundary" --legendLocation upper-right ) &
#icm和te的di就是低

n=te_test2
(computeMatrix reference-point -S $ddir/e2cell.bw $ddir/l2cell.bw $ddir/4cell.bw $ddir/8cell.bw $ddir/morula.bw $ddir/icm.bw $ddir/te.bw \
-R $bdir/te.bed -p 16 -a 1400000 -b 1400000 --referencePoint center -bs 40000 --skipZeros -o $wdir/$n.gz
plotProfile -m $wdir/$n.gz -out $wdir/$n.pdf --perGroup --regionsLabel "" --refPointLabel boundary \
-T "Directionality Index around TE TAD boundary" --legendLocation upper-right ) &
#icm和te的di就是低

n=te_test4 #05.31的数据
bdir=~/workspace/8.NT-HiC/g.DI_ALL/7.except_norm_mean_raw_20180531/2.unoverlap_tads
(computeMatrix reference-point -S $ddir/e2cell.bw $ddir/l2cell.bw $ddir/4cell.bw $ddir/8cell.bw $ddir/morula.bw $ddir/icm.bw $ddir/te.bw \
-R $bdir/icm.tad -p 16 -a 1400000 -b 1400000 --referencePoint center -bs 40000 --skipZeros -o $wdir/$n.gz
plotProfile -m $wdir/$n.gz -out $wdir/$n.pdf --perGroup --regionsLabel "" --refPointLabel boundary \
-T "Directionality Index around ICM TAD boundary" --legendLocation upper-right ) &
#icm和te的di就是低

n=morula_test2
(computeMatrix reference-point -S $ddir/e2cell.bw $ddir/l2cell.bw $ddir/4cell.bw $ddir/8cell.bw $ddir/morula.bw $ddir/icm.bw $ddir/te.bw \
-R $bdir/morula.bed -p 16 -a 1400000 -b 1400000 --referencePoint center -bs 40000 --skipZeros -o $wdir/$n.gz
plotProfile -m $wdir/$n.gz -out $wdir/$n.pdf --perGroup --regionsLabel "" --refPointLabel boundary \
-T "Directionality Index around TE TAD boundary" --legendLocation upper-right ) &


n=ICM_test
(computeMatrix reference-point -S $ddir/cc.bw $ddir/6h.bw $ddir/12h.bw $ddir/e2cell.bw \
-R $bdir/icm.bed -p 16 -a 1400000 -b 1400000 --referencePoint center -bs 40000 --skipZeros -o $wdir/$n.gz
plotProfile -m $wdir/$n.gz -out $wdir/$n.pdf --perGroup --regionsLabel "" --refPointLabel boundary \
-T "Directionality Index around ICM TAD boundary" --legendLocation upper-right \
--samplesLabel "Cumulus Cell" "Activated 6h" "Activated 12h" "Early 2cell" \
--colors '#e41a1c' '#4daf4a' '#984ea3' '#e5c494') &


n=CC_test
(computeMatrix reference-point -S $ddir/e2cell.bw $ddir/l2cell.bw $ddir/4cell.bw $ddir/8cell.bw $ddir/icm.bw \
-R $bdir/cc.bed -p 16 -a 1400000 -b 1400000 --referencePoint center -bs 40000 --skipZeros -o $wdir/$n.gz
plotProfile -m $wdir/$n.gz -out $wdir/$n.pdf --perGroup --regionsLabel "" --refPointLabel boundary \
-T "Directionality Index around CC TAD boundary" --legendLocation upper-right \
--colors '#66c2a5' '#fc8d62' '#8da0cb' '#e78ac3' '#a6d854' '#e5c494') &

#############################################################
#############################################################
#map to DI TAD border

