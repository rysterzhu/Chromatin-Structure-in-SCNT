wdir=~/workspace/8.NT-HiC/f.IS_ALL/f.map_IS_to_regions
bdir=~/workspace/8.NT-HiC/b.RNA/2.ZhangYi_analysis/1.identify_RRR/regions
ddir=~/workspace/8.NT-HiC/f.IS_ALL/2.except_res40k_is1M_ids200k_nt01/2.cat_insulation
ddir2=~/workspace/8.NT-HiC/3.public_data/f.XW_analysis/4.IS/2.cat_insulation
cp $ddir/[ciel]*bedGraph $wdir/1.IS
cp $ddir/12h.bedGraph $wdir/1.IS
cp $ddir2/ICM.bedGraph $wdir/1.IS/NF-ICM.bedGraph
cp $ddir2/early_2cell.bedGraph $wdir/1.IS/NF-e2cell.bedGraph
cp $ddir2/late_2cell.bedGraph $wdir/1.IS/NF-l2cell.bedGraph
cp $ddir2/PN5_zygote.bedGraph $wdir/1.IS/NF-PN5.bedGraph

for i in $wdir/1.IS/*bedGraph; do o=${i##*/};o=${o%%.*}
(sort -k1,1 -k2n,2 $i > $wdir/1.IS/$o.sorted.bg
bedGraphToBigWig $wdir/1.IS/$o.sorted.bg ~/ann/chromSizes/mm10.chrom.sizes $wdir/1.IS/$o.bw) &
done
ss=(CC NT-12h NF-PN5 NT-e2cell NF-e2cell NT-ICM NF-ICM)

for i in $wdir/1.IS/*bw; do k1=$(basename $i .bw)
for j in $bdir/*RR.bed; do k2=$(basename $j .bed)
computeMatrix scale-regions -S $i -R $j -p 32 -m 400000 -bs 40000 --missingDataAsZero --averageTypeBins mean -o $wdir/3.RRR/${k2}.${k1}.gz \
--startLabel "" --endLabel "" -a 2000000 -b 2000000 &
#plotHeatmap -m $wdir/3.RRR/${k2}.${k1}.gz -o $wdir/3.RRR/${k2}.${k1}.png --dpi 360 -T "" --startLabel "" --endLabel ""
done 
done


for i in $wdir/1.IS/*bw; do k1=$(basename $i .bw)
computeMatrix scale-regions -S $i -R $bdir/RRR.bed $bdir/PRR.bed $bdir/FRR.bed -p 32 -m 400000 -bs 40000 --missingDataAsZero --averageTypeBins mean -o $wdir/2.RRR/${k1}.gz \
--startLabel "" --endLabel "" -a 2000000 -b 2000000 
plotHeatmap -m $wdir/2.RRR/${k1}.gz -o $wdir/2.RRR/${k1}.png --dpi 360 -T "" --regionsLabel FRR PRR RRR --startLabel "" --endLabel ""
done &

##################################
#map to boundary
wdir=~/workspace/8.NT-HiC/f.IS_ALL/f.map_IS_to_regions/4.boundary
cp $ddir/../1.cat_boundary/[ic]*bed $wdir

for i in $wdir/../1.IS/*bw; do k1=$(basename $i .bw); 
for j in $wdir/*bed; do k2=$(basename $j .bed); 
computeMatrix reference-point -S $i -R $j -p 32 -bs 40000 --missingDataAsZero --averageTypeBins mean -o $wdir/${k2}.${k1}.gz --referencePoint center -a 2000000 -b 2000000 & 
done; done


########################################################
#Map to boundary and TAD center, random
wdir=~/workspace/8.NT-HiC/f.IS_ALL/f.map_IS_to_regions/5.TAD
tdir=~/workspace/8.NT-HiC/g.DI_ALL/2.except_norm_100M_20180825/2.unoverlap_tads

for i in TSS center TES; do 
for j in $wdir/../1.IS/*bw; do k=$(basename $j .bw); 
(#computeMatrix reference-point -S $j -R $tdir/cc.tad -p 32 -bs 40000 --missingDataAsZero --averageTypeBins mean -o $wdir/$k.$i.gz --referencePoint $i -a 2000000 -b 2000000
plotProfile -m $wdir/$k.$i.gz -o $wdir/$k.$i.png --dpi 360 -T $k.$i )&
done
done

#random
awk 'NR==FNR{a[$1]=$2} NR>FNR{srand($2);s=rand()*a[$1];printf("%s\t%d\t%d\n",$1,s,s+40000);srand($3);s=rand()*a[$1];printf("%s\t%d\t%d\n",$1,s,s+40000)}' ~/ann/mm10.shortChromSizes $tdir/cc.tad | sort -k1V,1 -k2n,2 > $wdir/random.bed
i=random
for j in $wdir/../1.IS/*bw; do k=$(basename $j .bw); 
(computeMatrix reference-point -S $j -R $wdir/random.bed -p 32 -bs 40000 --missingDataAsZero --averageTypeBins mean -o $wdir/$k.$i.gz --referencePoint TSS -a 2000000 -b 2000000
plotProfile -m $wdir/$k.$i.gz -o $wdir/$k.$i.png --dpi 360 -T $k.$i )&
done

