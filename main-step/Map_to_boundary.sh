wdir=~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/1.SINE
bdir=~/workspace/8.NT-HiC/f.IS_ALL/13.except_xw_para_20180816/cat_boundary

#SINE bed统计10k bins内的counts，转bw
bedtools coverage -a $ConfigHP/10000_mm10.bed -b ~/ann/mm10.SINE.bed | awk '$5>0{print $1,$2,$3,$5}' | sort -k1,1 -k2n,2 > mm10.SINE.bg
bedGraphToBigWig mm10.SINE.bg ~/ann/mm10.chromSizes mm10.SINE.bw


n=1
computeMatrix reference-point -S mm10.SINE.bw -R $bdir/icm.bed \
-p 50 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/$n.gz -a 500000 -b 500000 --referencePoint center; 



plotHeatmap -m $wdir/$n.gz -o $wdir/$n.h.png --dpi 360 -T "" --colorMap "bwr" --refPointLabel "ICM boundary" 
plotProfile -m $wdir/$n.gz -o $wdir/$n.f.png --dpi 360 --refPointLabel "ICM boundary" 

n=2
computeMatrix reference-point -S mm10.SINE.bw -R `for i in ${keys[@]};do echo -n "\$bdir/$i.bed ";done` \
--smartLabels -p 50 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/$n.gz -a 1000000 -b 1000000 --referencePoint center; 
plotHeatmap -m $wdir/$n.gz -o $wdir/$n.h.png --dpi 360 -T "map to boundary" --colorMap "bwr" --refPointLabel "boundary" --perGroup -z ${keys[@]}
plotProfile -m $wdir/$n.gz -o $wdir/$n.f.png --dpi 360 --refPointLabel "boundary" -z ${keys[@]}


#NF
n=nf1
bdir2=~/workspace/8.NT-HiC/3.public_data/f.XW_analysis/4.IS/cat_boundary
wdir=~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/1.SINE
computeMatrix reference-point -S $wdir/mm10.SINE.bw -R `for i in ${kNF[@]};do echo -n "\$bdir2/$i.bed ";done` \
--smartLabels -p 50 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/$n.gz -a 1000000 -b 1000000 --referencePoint center; 
plotHeatmap -m $wdir/$n.gz -o $wdir/$n.h.png --dpi 360 -T "map to boundary" --colorMap "bwr" --refPointLabel "boundary" --perGroup -z ${nNF[@]}

#intersect of NT NF boundary
#
cat_boundary:
for i in *bed; do awk '{$2=$2+100000;$3=$3-100000;print $1,$2,$3}' $i > bmoe05/$i;done

keyt=("cc" "6h" "12h" "e2cell" "l2cell" "4cell" "8cell" "icm")
namet=("CC" "A-6h" "A-12h" "NT-e2cell" "NT-l2cell" "NT-4cell" "NT-8cell" "NT-ICM")

keyf=("Sperm" "PN3_zygote" "PN5_zygote" "early_2cell" "late_2cell" "4cell" "8cell" "ICM")
namef=("Sperm" "PN3-zygote" "PN5-zygote" "NF-e2cell" "NF-l2cell" "NF-4cell" "NF-8cell" "NF-ICM")
nametf=("CC-Sperm" "6h-PN3" "12h-PN5" "e2cell" "l2cell" "4cell" "8cell" "ICM")

bdir2=~/workspace/8.NT-HiC/3.public_data/f.XW_analysis/4.IS/1.cat_boundary
wdir=~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/1.SINE
bdir1=~/workspace/8.NT-HiC/f.IS_ALL/13.except_xw_para_20180816/1.cat_boundary
for i in {0..7}; do 
bedtools intersect -a $bdir1/${keyt[$i]}.bed -b $bdir2/${keyf[$i]}.bed -r -f 0.999 > $wdir/a.intersect_NT_NF/${nametf[$i]}.bed
bedtools intersect -a $bdir1/${keyt[$i]}.bed -b $bdir2/${keyf[$i]}.bed -r -f 0.999 -v > $wdir/a.intersect_NT_NF/${namet[$i]}.bed
bedtools intersect -a $bdir2/${keyf[$i]}.bed -b $bdir1/${keyt[$i]}.bed -r -f 0.999 -v > $wdir/a.intersect_NT_NF/${namef[$i]}.bed
done
for i in {0..7}; do 
wc -l $bdir1/${keyt[$i]}.bed $bdir2/${keyf[$i]}.bed $wdir/a.intersect_NT_NF/${namet[$i]}.bed $wdir/a.intersect_NT_NF/${namef[$i]}.bed $wdir/a.intersect_NT_NF/${nametf[$i]}.bed
done
  2852 /home/qszhu/workspace/8.NT-HiC/f.IS_ALL/13.except_xw_para_20180816/1.cat_boundary/icm.bed
  3000 /home/qszhu/workspace/8.NT-HiC/3.public_data/f.XW_analysis/4.IS/1.cat_boundary/ICM.bed
  1761 /home/qszhu/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/1.SINE/a.intersect_NT_NF/NT-ICM.bed
  1909 /home/qszhu/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/1.SINE/a.intersect_NT_NF/NF-ICM.bed
  1091 /home/qszhu/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/1.SINE/a.intersect_NT_NF/ICM.bed
#
wdir=~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/1.SINE/
for i in {0..7}; do 
(computeMatrix reference-point -S $wdir/mm10.SINE.bw -R $wdir/a.intersect_NT_NF/${namet[$i]}.bed $wdir/a.intersect_NT_NF/${nametf[$i]}.bed $wdir/a.intersect_NT_NF/${namef[$i]}.bed \
--smartLabels -p 8 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/a.intersect_NT_NF/${nametf[$i]}.gz -a 1000000 -b 1000000 --referencePoint center; 
plotHeatmap -m $wdir/a.intersect_NT_NF/${nametf[$i]}.gz -o $wdir/a.intersect_NT_NF/png/${nametf[$i]}.h.png --dpi 360 -T "" --colorMap "bwr" --refPointLabel "boundary" -z NT NT-NF NF --samplesLabel ${nametf[$i]}) &
done

#intersect of NT NF boundary: 1bp
for i in {0..7}; do 
bedtools intersect -a $bdir1/${keyt[$i]}.bed -b $bdir2/${keyf[$i]}.bed | sort -k1,1 -k2n,2| bedtools merge -d 120000 -i - > $wdir/b.intersect_NT_NF/${nametf[$i]}.bed
bedtools intersect -a $bdir1/${keyt[$i]}.bed -b $bdir2/${keyf[$i]}.bed -v > $wdir/b.intersect_NT_NF/${namet[$i]}.bed
bedtools intersect -a $bdir2/${keyf[$i]}.bed -b $bdir1/${keyt[$i]}.bed -v > $wdir/b.intersect_NT_NF/${namef[$i]}.bed
done
for i in {0..7}; do 
wc -l $bdir1/${keyt[$i]}.bed $bdir2/${keyf[$i]}.bed $wdir/b.intersect_NT_NF/${namet[$i]}.bed $wdir/b.intersect_NT_NF/${namef[$i]}.bed $wdir/b.intersect_NT_NF/${nametf[$i]}.bed
done

##########################################################
##Map HM ChIP-seq to boundary : K9me3 K4me3 K27me3
#log2 ratio of the total ChIP-seq signal over the 10kb window divided by the input signal of the window
ddir=~/workspace/9.NT-ChIP/2.public/3.bowtie2/c.merge_rep
wdir=~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/2.K9me3
bdir1=~/workspace/8.NT-HiC/f.IS_ALL/13.except_xw_para_20180816/1.cat_boundary
bdir2=~/workspace/8.NT-HiC/3.public_data/f.XW_analysis/4.IS/1.cat_boundary
bdir3=~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/1.SINE/a.intersect_NT_NF

for i in $ddir/*K9me3.sorted.bam; do key=${i##*/};key=${key%%.*}
nohup bamCompare -b1 $i -b2 $ddir/${key%%.*}.input.sorted.bam -o $wdir/0.bamCompare/$key.bw -of bigwig --scaleFactorsMethod None --operation log2 --pseudocount 1 \
--binSize 10000 -p 8 --normalizeUsing RPKM --ignoreDuplicates --samFlagInclude 2 > $wdir/logs/$key.log &
done 

###test:
i=$ddir/ICM.K27me3.sorted.bam;key=${i##*/};key=${key%%.*}
nohup bamCompare -b1 $i -b2 $ddir/$key.input.sorted.bam -o $wdir/test/$key.1.bedGraph -of bedgraph --scaleFactorsMethod None --operation log2 --pseudocount 1 \
--binSize 10000 -p 36 --normalizeUsing RPKM --ignoreDuplicates --samFlagInclude 2 > $wdir/test/$key.1.log &
nohup bamCompare -b1 $i -b2 $ddir/$key.input.sorted.bam -o $wdir/test/$key.2.bedGraph -of bedgraph --scaleFactorsMethod None --operation log2 --pseudocount 10 \
--binSize 10000 -p 36 --normalizeUsing RPKM --ignoreDuplicates --samFlagInclude 2 > $wdir/test/$key.2.log &
nohup bamCompare -b1 $i -b2 $ddir/$key.input.sorted.bam -o $wdir/test/$key.3.bedGraph -of bedgraph --scaleFactorsMethod None --operation log2 --pseudocount 1 \
--binSize 10000 -p 36 --normalizeUsing CPM --ignoreDuplicates --samFlagInclude 2 > $wdir/test/$key.3.log &
nohup bamCompare -b1 $i -b2 $ddir/$key.input.sorted.bam -o $wdir/test/$key.4.bedGraph -of bedgraph --scaleFactorsMethod None --operation log2 --pseudocount 1 \
--binSize 10000 --smoothLength 30000 -p 36 --normalizeUsing RPKM --ignoreDuplicates --samFlagInclude 2 > $wdir/test/$key.4.log &
#区别很大，但不知道哪个好

#####heatmap
keyt=("cc" "6h" "12h" "e2cell" "l2cell" "4cell" "8cell" "icm")
namet=("CC" "A-6h" "A-12h" "NT-e2cell" "NT-l2cell" "NT-4cell" "NT-8cell" "NT-ICM")

keyf=("Sperm" "PN3_zygote" "PN5_zygote" "early_2cell" "late_2cell" "4cell" "8cell" "ICM")
namef=("Sperm" "PN3-zygote" "PN5-zygote" "NF-e2cell" "NF-l2cell" "NF-4cell" "NF-8cell" "NF-ICM")
nametf=("CC-Sperm" "6h-PN3" "12h-PN5" "e2cell" "l2cell" "4cell" "8cell" "ICM")

chip=(zygote 2cell 4cell 8cell morula TE ICM)

#NT
for n in {0..7}; do 
(computeMatrix reference-point -S `for i in ${chip[@]}; do echo -n "\$wdir/0.bamCompare/$i.bw ";done` -R $bdir1/${keyt[$n]}.bed \
--smartLabels -p 8 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/1.NT/${keyt[$n]}.gz -a 1000000 -b 1000000 --referencePoint center; 
plotHeatmap -m $wdir/1.NT/${keyt[$n]}.gz -o $wdir/1.NT/${keyt[$n]}.h.png --dpi 360 -T ${namet[$n]} --colorMap "bwr" --refPointLabel "boundary" --samplesLabel ${chip[@]}
plotProfile -m $wdir/1.NT/${keyt[$n]}.gz -o $wdir/1.NT/${keyt[$n]}.f.png --dpi 360 -T ${namet[$n]} --refPointLabel "boundary" --samplesLabel ${chip[@]} )&
done
#NF
for n in {0..7}; do 
(computeMatrix reference-point -S `for i in ${chip[@]}; do echo -n "\$wdir/0.bamCompare/$i.bw ";done` -R $bdir2/${keyf[$n]}.bed \
--smartLabels -p 8 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/2.NF/${keyf[$n]}.gz -a 1000000 -b 1000000 --referencePoint center; 
plotHeatmap -m $wdir/2.NF/${keyf[$n]}.gz -o $wdir/2.NF/${keyf[$n]}.h.png --dpi 360 -T ${namef[$n]} --colorMap "bwr" --refPointLabel "boundary" --samplesLabel ${chip[@]}
plotProfile -m $wdir/2.NF/${keyf[$n]}.gz -o $wdir/2.NF/${keyf[$n]}.f.png --dpi 360 -T ${namef[$n]} --refPointLabel "boundary" --samplesLabel ${chip[@]} )&
done
###看起来还是zygote的NT和NF有一些差异

##intersect NT NF
for i in {0..7}; do 
(#computeMatrix reference-point -S `for i in ${chip[@]}; do echo -n "\$wdir/0.bamCompare/$i.bw ";done` -R $bdir3/${namet[$i]}.bed $bdir3/${nametf[$i]}.bed $bdir3/${namef[$i]}.bed \
#--smartLabels -p 8 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/3.NT-NF/${nametf[$i]}.gz -a 1000000 -b 1000000 --referencePoint center; 
plotHeatmap -m $wdir/3.NT-NF/${nametf[$i]}.gz -o $wdir/3.NT-NF/png/${nametf[$i]}.h.png --dpi 360 -T "" --colorMap "bwr" --refPointLabel "boundary" -z NT NT-NF NF --samplesLabel ${chip[@]}) &
done

#HM 按HM分，每种HM在所有时期的boundary上的图
for n in {0..6}; do 
(computeMatrix reference-point -S $wdir/0.bamCompare/${chip[$n]}.bw -R `for i in ${keyt[@]};do echo -n "\$bdir1/$i.bed ";done;for i in ${keyf[@]};do echo -n "\$bdir2/$i.bed ";done` \
--smartLabels -p 8 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/4.HM/${chip[$n]}.gz -a 1000000 -b 1000000 --referencePoint center; 
plotHeatmap -m $wdir/4.HM/${chip[$n]}.gz -o $wdir/4.HM/${chip[$n]}.h.png --dpi 360 -T ${chip[$n]} --colorMap "bwr" --refPointLabel "boundary" --perGroup -z ${namet[@]} ${namef[@]}
plotProfile -m $wdir/4.HM/${chip[$n]}.gz -o $wdir/4.HM/${chip[$n]}.f.png --dpi 360 -T ${chip[$n]} --refPointLabel "boundary" --perGroup -z ${namet[@]} ${namef[@]} )&
done

#################################################
##########################K27me3 K4me3
ddir=~/workspace/9.NT-ChIP/2.public/3.bowtie2/c.merge_rep
wdir=~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/4.K27me3   #3.K4me3
bdir1=~/workspace/8.NT-HiC/f.IS_ALL/13.except_xw_para_20180816/1.cat_boundary
bdir2=~/workspace/8.NT-HiC/3.public_data/f.XW_analysis/4.IS/1.cat_boundary
bdir3=~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/1.SINE/a.intersect_NT_NF

mkdir -p $wdir/0.bamCompare $wdir/1.NT $wdir/2.NF $wdir/3.NT-NF/png $wdir/logs $wdir/4.HM
for i in $ddir/*K27me3.sorted.bam; do key=${i##*/};key=${key%%.*}
nohup bamCompare -b1 $i -b2 $ddir/${key%%.*}.input.sorted.bam -o $wdir/0.bamCompare/$key.bw -of bigwig --scaleFactorsMethod None --operation log2 --pseudocount 1 \
--binSize 10000 -p 8 --normalizeUsing RPKM --ignoreDuplicates --samFlagInclude 2 > $wdir/logs/$key.log &
done 
wait
#####heatmap
keyt=("cc" "6h" "12h" "e2cell" "l2cell" "4cell" "8cell" "icm")
namet=("CC" "A-6h" "A-12h" "NT-e2cell" "NT-l2cell" "NT-4cell" "NT-8cell" "NT-ICM")

keyf=("Sperm" "PN3_zygote" "PN5_zygote" "early_2cell" "late_2cell" "4cell" "8cell" "ICM")
namef=("Sperm" "PN3-zygote" "PN5-zygote" "NF-e2cell" "NF-l2cell" "NF-4cell" "NF-8cell" "NF-ICM")
nametf=("CC-Sperm" "6h-PN3" "12h-PN5" "e2cell" "l2cell" "4cell" "8cell" "ICM")

chip=(zygote 2cell 4cell 8cell morula TE ICM)

#NT
for n in {0..7}; do 
(computeMatrix reference-point -S `for i in ${chip[@]}; do echo -n "\$wdir/0.bamCompare/$i.bw ";done` -R $bdir1/${keyt[$n]}.bed \
--smartLabels -p 8 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/1.NT/${keyt[$n]}.gz -a 1000000 -b 1000000 --referencePoint center; 
plotHeatmap -m $wdir/1.NT/${keyt[$n]}.gz -o $wdir/1.NT/${keyt[$n]}.h.png --dpi 360 -T ${namet[$n]} --colorMap "bwr" --refPointLabel "boundary" --samplesLabel ${chip[@]}
plotProfile -m $wdir/1.NT/${keyt[$n]}.gz -o $wdir/1.NT/${keyt[$n]}.f.png --dpi 360 -T ${namet[$n]} --refPointLabel "boundary" --samplesLabel ${chip[@]} )&
done
wait
#NF
for n in {0..7}; do 
(computeMatrix reference-point -S `for i in ${chip[@]}; do echo -n "\$wdir/0.bamCompare/$i.bw ";done` -R $bdir2/${keyf[$n]}.bed \
--smartLabels -p 8 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/2.NF/${keyf[$n]}.gz -a 1000000 -b 1000000 --referencePoint center; 
plotHeatmap -m $wdir/2.NF/${keyf[$n]}.gz -o $wdir/2.NF/${keyf[$n]}.h.png --dpi 360 -T ${namef[$n]} --colorMap "bwr" --refPointLabel "boundary" --samplesLabel ${chip[@]}
plotProfile -m $wdir/2.NF/${keyf[$n]}.gz -o $wdir/2.NF/${keyf[$n]}.f.png --dpi 360 -T ${namef[$n]} --refPointLabel "boundary" --samplesLabel ${chip[@]} )&
done
wait
#HM
for n in {0..6}; do 
(computeMatrix reference-point -S $wdir/0.bamCompare/${chip[$n]}.bw -R `for i in ${keyt[@]};do echo -n "\$bdir1/$i.bed ";done;for i in ${keyf[@]};do echo -n "\$bdir2/$i.bed ";done` \
--smartLabels -p 8 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/4.HM/${chip[$n]}.gz -a 1000000 -b 1000000 --referencePoint center; 
plotHeatmap -m $wdir/4.HM/${chip[$n]}.gz -o $wdir/4.HM/${chip[$n]}.h.png --dpi 360 -T ${chip[$n]} --colorMap "bwr" --refPointLabel "boundary" --perGroup -z ${namet[@]} ${namef[@]}
plotProfile -m $wdir/4.HM/${chip[$n]}.gz -o $wdir/4.HM/${chip[$n]}.f.png --dpi 360 -T ${chip[$n]} --refPointLabel "boundary" --perGroup -z ${namet[@]} ${namef[@]} )&
done
wait
###看起来还是zygote的NT和NF有一些差异

##intersect NT NF
for i in {0..7}; do 
(computeMatrix reference-point -S `for i in ${chip[@]}; do echo -n "\$wdir/0.bamCompare/$i.bw ";done` -R $bdir3/${namet[$i]}.bed $bdir3/${nametf[$i]}.bed $bdir3/${namef[$i]}.bed \
--smartLabels -p 8 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/3.NT-NF/${nametf[$i]}.gz -a 1000000 -b 1000000 --referencePoint center; 
plotHeatmap -m $wdir/3.NT-NF/${nametf[$i]}.gz -o $wdir/3.NT-NF/png/${nametf[$i]}.h.png --dpi 360 -T "" --colorMap "bwr" --refPointLabel "boundary" -z NT NT-NF NF --samplesLabel ${chip[@]}) &
done
##主要还是zygote到2cell的差异

#########################
#NT K9me3 
########################
ddir=~/workspace/9.NT-ChIP/1.align/[cdf]*/a.markDuplicates
wdir=~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/5.NT_K9me3

for i in 6h 14h 2cell; do 
samtools merge $wdir/00.bam/$i.sorted.bam $ddir/${i}_K9me3*bam &
samtools merge -@ 16 $wdir/00.bam/${i}_input.sorted.bam $ddir/${i}_input*bam &
done

bdir1=~/workspace/8.NT-HiC/f.IS_ALL/13.except_xw_para_20180816/1.cat_boundary
bdir2=~/workspace/8.NT-HiC/3.public_data/f.XW_analysis/4.IS/1.cat_boundary
bdir3=~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/1.SINE/a.intersect_NT_NF
keyt=("cc" "6h" "12h" "e2cell" "l2cell" "4cell" "8cell" "icm")
namet=("CC" "A-6h" "A-12h" "NT-e2cell" "NT-l2cell" "NT-4cell" "NT-8cell" "NT-ICM")
keyf=("Sperm" "PN3_zygote" "PN5_zygote" "early_2cell" "late_2cell" "4cell" "8cell" "ICM")
namef=("Sperm" "PN3-zygote" "PN5-zygote" "NF-e2cell" "NF-l2cell" "NF-4cell" "NF-8cell" "NF-ICM")
nametf=("CC-Sperm" "6h-PN3" "12h-PN5" "e2cell" "l2cell" "4cell" "8cell" "ICM")
chip=(6h 14h 2cell)

mkdir -p $wdir/0.bamCompare $wdir/1.NT $wdir/2.NF $wdir/3.NT-NF/png $wdir/logs $wdir/4.HM
for key in ${chip[@]}; do key=${i##*/};key=${key%%.*}
nohup bamCompare -b1 $wdir/00.bam/${key}.sorted.bam -b2 $wdir/00.bam/${key}_input.sorted.bam -o $wdir/0.bamCompare/$key.bw -of bigwig --scaleFactorsMethod None --operation log2 --pseudocount 1 \
--binSize 10000 -p 8 --normalizeUsing RPKM --ignoreDuplicates --samFlagInclude 2 > $wdir/logs/$key.log &
done 
wait
#HM
for n in {0..2}; do 
(computeMatrix reference-point -S $wdir/0.bamCompare/${chip[$n]}.bw -R `for i in ${keyt[@]};do echo -n "\$bdir1/$i.bed ";done;for i in ${keyf[@]};do echo -n "\$bdir2/$i.bed ";done` \
--smartLabels -p 8 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/4.HM/${chip[$n]}.gz -a 1000000 -b 1000000 --referencePoint center; 
plotHeatmap -m $wdir/4.HM/${chip[$n]}.gz -o $wdir/4.HM/${chip[$n]}.h.png --dpi 360 -T ${chip[$n]} --colorMap "bwr" --refPointLabel "boundary" --perGroup -z ${namet[@]} ${namef[@]}
plotProfile -m $wdir/4.HM/${chip[$n]}.gz -o $wdir/4.HM/${chip[$n]}.f.png --dpi 360 -T ${chip[$n]} --refPointLabel "boundary" --perGroup -z ${namet[@]} ${namef[@]} )&
done

#####NT K9 map to NT NF PN3
n=1
(computeMatrix reference-point -S $wdir/0.bamCompare/6h.bw -R $bdir1/6h.bed $bdir2/PN3_zygote.bed \
--smartLabels -p 8 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/4.HM/$n.gz -a 1000000 -b 1000000 --referencePoint center; 
plotProfile -m $wdir/4.HM/$n.gz -o $wdir/4.HM/$n.f.png --dpi 360 -T "NT 6h K9 map to NT&NF boundary" --refPointLabel "boundary" -z NT-6h NF-PN3 )&

n=2
(computeMatrix reference-point -S ~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/2.K9me3/0.bamCompare/zygote.bw -R $bdir1/6h.bed $bdir2/PN3_zygote.bed \
--smartLabels -p 8 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/4.HM/$n.gz -a 1000000 -b 1000000 --referencePoint center; 
plotProfile -m $wdir/4.HM/$n.gz -o $wdir/4.HM/$n.f.png --dpi 360 -T "NF zygote K9 map to NT&NF boundary" --refPointLabel "boundary" -z NT-6h NF-PN3 )&

bdir3=~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/1.SINE/b.intersect_NT_NF
n=3
computeMatrix reference-point -S $wdir/0.bamCompare/6h.bw -R $bdir3/A-6h.bed $bdir3/6h-PN3.bed $bdir3/PN3-zygote.bed \
--smartLabels -p 8 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/4.HM/$n.gz -a 1000000 -b 1000000 --referencePoint center; 
plotProfile -m $wdir/4.HM/$n.gz -o $wdir/4.HM/$n.f.png --dpi 360 -T "" --refPointLabel "boundary" -z NT both NF --samplesLabel "NT-6h K9me3"&
n=4
(computeMatrix reference-point -S ~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/2.K9me3/0.bamCompare/zygote.bw -R $bdir3/A-6h.bed $bdir3/6h-PN3.bed $bdir3/PN3-zygote.bed \
--smartLabels -p 8 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/4.HM/$n.gz -a 1000000 -b 1000000 --referencePoint center; 
plotProfile -m $wdir/4.HM/$n.gz -o $wdir/4.HM/$n.f.png --dpi 360 -T "" --refPointLabel "boundary" -z NT both NF --samplesLabel "NF-zygote K9me3")&

plotProfile -m $wdir/4.HM/3.gz -o $wdir/4.HM/3.f.png --dpi 360 -T "" --refPointLabel "boundary" -z NT-special both NF-special --samplesLabel "NT-6h K9me3" --yMin -0.07 --yMax 0.05 &
plotProfile -m $wdir/4.HM/4.gz -o $wdir/4.HM/4.f.png --dpi 360 -T "" --refPointLabel "boundary" -z NT-special both NF-special --samplesLabel "NF-zygote K9me3" --yMin -0.07 --yMax 0.05 &

n=5
(computeMatrix reference-point -S $wdir/0.bamCompare/2cell.bw -R $bdir3/NT-l2cell.bed $bdir3/l2cell.bed $bdir3/NF-l2cell.bed \
--smartLabels -p 8 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/4.HM/$n.gz -a 1000000 -b 1000000 --referencePoint center; 
plotProfile -m $wdir/4.HM/$n.gz -o $wdir/4.HM/$n.f.png --dpi 360 -T "" --refPointLabel "boundary" -z NT-special both NF-special --samplesLabel "NT-2cell K9me3")&
n=6
(computeMatrix reference-point -S ~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/2.K9me3/0.bamCompare/zygote.bw -R $bdir3/NT-l2cell.bed $bdir3/l2cell.bed $bdir3/NF-l2cell.bed \
--smartLabels -p 8 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/4.HM/$n.gz -a 1000000 -b 1000000 --referencePoint center; 
plotProfile -m $wdir/4.HM/$n.gz -o $wdir/4.HM/$n.f.png --dpi 360 -T "" --refPointLabel "boundary" -z NT-special both NF-special --samplesLabel "NF-2cell K9me3")&


##############################################
#2018年8月29日 Count the SINE class in boundary
################################################
wdir=~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/1.SINE/c.count_SINE_IN_boundary
mkdir $wdir/0.SINE
awk '{print $1,$2,$3 >> "0.SINE/"$4}' ~/ann/mm10.SINE.bed

#nt=0.1,bmoe=1,PN3 + 6h  #bmoe too low
awk -v b=2 '$5>0.1{print $1,$2-40000*b,$3+40000*b}' ~/workspace/8.NT-HiC/f.IS_ALL/2.except_res40k_is1M_ids200k_nt01/1.cat_boundary/l2cell.bed > NT.bed
awk -v b=2 '$5>0.1{print $1,$2-40000*(b-1),$3+40000*(b-1)}' ~/workspace/8.NT-HiC/3.public_data/f.XW_analysis/4.IS/1.cat_boundary/late_2cell.bed > NF.bed
bedtools intersect -a NT.bed -b NF.bed > inter.bed
bedtools intersect -a NT.bed -b inter.bed -v > NT-v.bed
bedtools intersect -a NF.bed -b inter.bed -v > NF-v.bed
cat NT.bed NF.bed | sort -k1,1 -k2n,2 | bedtools merge -i - > merge.bed

#########
echo -ne "SINE\ttotal\tNT\tNF\tinter\n" > count-6h.tab
for i in $wdir/0.SINE/*; do k=$(basename $i)
total=`awk '{n+=$3-$2} END{print n}' $i`
NT=`bedtools intersect -a NT-v.bed -b $i -c | awk '{n+=$4} END{print n}'`
NF=`bedtools intersect -a NF-v.bed -b $i -c | awk '{n+=$4} END{print n}'`
inter=`bedtools intersect -a inter.bed -b $i -c | awk '{n+=$4} END{print n}'`
echo -ne "${k}\t${total}\t${NT}\t${NF}\t${inter}\n" >> count-6h.tab
done
######## too low

############enrichment test
bedtools shuffle -i inter.bed -g ~/ann/mm10.chromSizes > r-inter.bed
bedtools intersect -a inter.bed -b ~/ann/mm10.SINE.bed -c | awk '{n+=$4} END{print n}'
bedtools intersect -a r-inter.bed -b ~/ann/mm10.SINE.bed -c | awk '{n+=$4} END{print n}'

bedtools shuffle -i inter.bed -g ~/ann/mm10.chromSizes | bedtools intersect -a - -b ~/ann/mm10.SINE.bed -c | \
awk -v t=`bedtools intersect -a inter.bed -b ~/ann/mm10.SINE.bed -c | awk '{n+=$4} END{print n}'` '{n+=$4} END{print log(t/n)/log(2)}'

bedtools shuffle -i NF.bed -g ~/ann/mm10.chromSizes | bedtools intersect -a - -b ~/ann/mm10.SINE.bed -c | \
awk -v t=`bedtools intersect -a NF.bed -b ~/ann/mm10.SINE.bed -c | awk '{n+=$4} END{print n}'` '{n+=$4} END{print log(t/n)/log(2)}'

bedtools shuffle -i NT.bed -g ~/ann/mm10.chromSizes | bedtools intersect -a - -b ~/ann/mm10.SINE.bed -c | \
awk -v t=`bedtools intersect -a NT.bed -b ~/ann/mm10.SINE.bed -c | awk '{n+=$4} END{print n}'` '{n+=$4} END{print log(t/n)/log(2)}'

bedtools shuffle -i merge.bed -g ~/ann/mm10.chromSizes | bedtools intersect -a - -b ~/ann/mm10.SINE.bed -c | \
awk -v t=`bedtools intersect -a merge.bed -b ~/ann/mm10.SINE.bed -c | awk '{n+=$4} END{print n}'` '{n+=$4} END{print log(t/n)/log(2)}'
###################good
#
seed=100
echo -ne "SINE\ttotal\tNT\tNF\tinter\tmerge\tenrich\n" > count-l2.tab
for i in $wdir/0.SINE/*; do k=$(basename $i)
total=`awk '{n+=$3-$2} END{print n}' $i`
NT=`bedtools intersect -a NT-v.bed -b $i -c | awk '{n+=$4} END{print n}'`
NF=`bedtools intersect -a NF-v.bed -b $i -c | awk '{n+=$4} END{print n}'`
inter=`bedtools intersect -a inter.bed -b $i -c | awk '{n+=$4} END{print n}'`
merge=`bedtools intersect -a merge.bed -b $i -c | awk '{n+=$4} END{print n}'`
enrich=`bedtools shuffle -seed $seed -i merge.bed -g ~/ann/mm10.chromSizes | bedtools intersect -a - -b $i -c | awk -v t=$merge '{n+=$4} END{print log(t/n)/log(2)}'`
echo -ne "${k}\t${total}\t${NT}\t${NF}\t${inter}\t${merge}\t${enrich}\n" >> count-l2.tab
done


###############################
#some map to cc icm ICM boundary
###############################
wdir=~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/6.cmm
ddir=~/workspace/8.NT-HiC/f.IS_ALL/2.except_res40k_is1M_ids200k_nt01/1.cat_boundary

#0.25的太少了，不行
awk '$5>0.1{print $1,$2-1,$3}' $ddir/cc.bed > $wdir/cc.bed
awk '$5>0.1{print $1,$2-1,$3}' $ddir/icm.bed > $wdir/icm.bed
awk '$5>0.1{print $1,$2+40000-1,$3-40000}' ~/workspace/8.NT-HiC/3.public_data/f.XW_analysis/4.IS/1.cat_boundary/ICM.bed > $wdir/ICM.bed

intersectBed -a icm.bed -b ICM.bed -wa -u > NT_NF.bed
intersectBed -a icm.bed -b ICM.bed -wa -v > NT_spe.bed

  3006 cc.bed
  2852 icm.bed
  3000 ICM.bed
  1091 NT_NF.bed
  1761 NT_spe.bed

#直接的SINE转bw
awk '{print $1,$2,$3,1}' ~/ann/mm10.SINE.bed | sort -k1,1 -k2n,2 |mergeBed -i - -o count -c 4> mm10.SINE.bg
bedGraphToBigWig mm10.SINE.bg ~/ann/mm10.chromSizes mm10.SINE.bw
#这种方式不行，还是用原来的方法
cp ../1.SINE/mm10.SINE.bw ./

(n=SINE
computeMatrix reference-point -S mm10.SINE.bw -R $wdir/*bed -p 32 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/$n.gz -a 500000 -b 500000 --referencePoint center; 
plotProfile -m $wdir/$n.gz -o $wdir/$n.f.png --dpi 360 --refPointLabel "boundary" --perGroup )&

mkdir -p $wdir/1.ggplot
for i in $wdir/*bed; do k=$(basename $i .bed)
computeMatrix reference-point -S mm10.SINE.bw -R $i -p 16 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/1.ggplot/SINE-$k.gz -a 500000 -b 500000 --referencePoint center &
done
###################K9
ddir=~/workspace/9.NT-ChIP/1.align/h.20180909/a.markDuplicates
wdir=~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/6.cmm/2.cc_K9
samtools merge -@ 8 $wdir/cc_K9me3.sorted.bam $ddir/cc_K9me3_*.sorted.bam & 
samtools merge -@ 8 $wdir/cc_input.sorted.bam $ddir/cc_input_rep[12].sorted.bam & 
wait
samtools index -@ 8 $wdir/cc_K9me3.sorted.bam &
samtools index -@ 8 $wdir/cc_input.sorted.bam &
wait
nohup bamCompare -b1 $wdir/cc_K9me3.sorted.bam -b2 $wdir/cc_input.sorted.bam -o $wdir/cc_K9me3_ratio.bw -of bigwig \
--scaleFactorsMethod None --operation ratio --pseudocount 1 \
--binSize 10000 -p 32 --normalizeUsing RPKM --ignoreDuplicates --samFlagInclude 2 &
wait
(n=K9_ratio
computeMatrix reference-point -S $wdir/cc_K9me3_ratio.bw -R $wdir/../*bed -p 32 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/$n.gz -a 500000 -b 500000 --referencePoint center; 
plotProfile -m $wdir/$n.gz -o $wdir/$n.f.png --dpi 360 --refPointLabel "boundary" --perGroup )&
#ratio和log2一样，单位不一样，使用ratio

for i in $wdir/../*bed; do k=$(basename $i .bed)
computeMatrix reference-point -S $wdir/cc_K9me3_ratio.bw -R $i -p 16 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/TSS-$k.gz -a 500000 -b 500000 --referencePoint center &
done

#########################TSS
wdir=~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/6.cmm/3.TSS

bedtools coverage -a $ConfigHP/10000_mm10.bed -b ~/ann/mm10_promoter_1000-1000.bed | awk '$5>0{print $1,$2,$3,$5}' | sort -k1,1 -k2n,2 > mm10.TSS.bg
bedGraphToBigWig mm10.TSS.bg ~/ann/mm10.chromSizes mm10.TSS.bw

(n=TSS
computeMatrix reference-point -S $wdir/mm10.TSS.bw -R $wdir/../*bed -p 32 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/$n.gz -a 500000 -b 500000 --referencePoint center; 
plotProfile -m $wdir/$n.gz -o $wdir/$n.f.png --dpi 360 --refPointLabel "boundary" --perGroup )&

for i in $wdir/../*bed; do k=$(basename $i .bed)
computeMatrix reference-point -S $wdir/mm10.TSS.bw -R $i -p 16 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/TSS-$k.gz -a 500000 -b 500000 --referencePoint center &
done

###################K4
ddir=~/workspace/9.NT-ChIP/1.align/h.20180909/a.markDuplicates
wdir=~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/6.cmm/4.cc_K4
samtools merge -@ 8 $wdir/cc_K4me3.sorted.bam $ddir/cc_K4me3_*.sorted.bam & 
samtools merge -@ 8 $wdir/cc_input.sorted.bam $ddir/cc_input_rep[12].sorted.bam & 
wait
samtools index -@ 8 $wdir/cc_K4me3.sorted.bam &
samtools index -@ 8 $wdir/cc_input.sorted.bam &
wait
nohup bamCompare -b1 $wdir/cc_K4me3.sorted.bam -b2 $wdir/cc_input.sorted.bam -o $wdir/cc_K4me3_ratio.bw -of bigwig \
--scaleFactorsMethod None --operation ratio --pseudocount 1 \
--binSize 10000 -p 8 --normalizeUsing RPKM --ignoreDuplicates --samFlagInclude 2 &
wait
(n=K4_ratio
computeMatrix reference-point -S $wdir/cc_K4me3_ratio.bw -R $wdir/../*bed -p 32 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/$n.gz -a 500000 -b 500000 --referencePoint center; 
plotProfile -m $wdir/$n.gz -o $wdir/$n.f.png --dpi 360 --refPointLabel "boundary" --perGroup )&

for i in $wdir/../*bed; do k=$(basename $i .bed)
computeMatrix reference-point -S $wdir/cc_K4me3_ratio.bw -R $i -p 16 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/TSS-$k.gz -a 500000 -b 500000 --referencePoint center &
done

################K27
wdir=~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/6.cmm/5.cc_K27
samtools merge -@ 8 $wdir/cc_K27me3.sorted.bam $ddir/cc_K27me3_*.sorted.bam & 
samtools merge -@ 8 $wdir/cc_input.sorted.bam $ddir/cc_input_rep[12].sorted.bam & 
wait
samtools index -@ 8 $wdir/cc_K27me3.sorted.bam &
samtools index -@ 8 $wdir/cc_input.sorted.bam &
wait
nohup bamCompare -b1 $wdir/cc_K27me3.sorted.bam -b2 $wdir/cc_input.sorted.bam -o $wdir/cc_K27me3_ratio.bw -of bigwig \
--scaleFactorsMethod None --operation ratio --pseudocount 1 \
--binSize 10000 -p 32 --normalizeUsing RPKM --ignoreDuplicates --samFlagInclude 2 &
wait

################LINE
wdir=~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/6.cmm/6.LINE
mkdir -p $wdir
bedtools coverage -a $ConfigHP/10000_mm10.bed -b ~/ann/mm10.LINE.bed | awk '$5>0{print $1,$2,$3,$5}' | sort -k1,1 -k2n,2 > $wdir/mm10.LINE.bg
bedGraphToBigWig $wdir/mm10.LINE.bg ~/ann/mm10.chromSizes $wdir/mm10.LINE.bw


################LTR
wdir=~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/6.cmm/7.LTR
mkdir -p $wdir
bedtools coverage -a $ConfigHP/10000_mm10.bed -b ~/ann/mm10.LTR.bed | awk '$5>0{print $1,$2,$3,$5}' | sort -k1,1 -k2n,2 > $wdir/mm10.LTR.bg
bedGraphToBigWig $wdir/mm10.LTR.bg ~/ann/mm10.chromSizes $wdir/mm10.LTR.bw


################CGI
wdir=~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/6.cmm/8.CGI
mkdir -p $wdir
bedtools coverage -a $ConfigHP/10000_mm10.bed -b ~/ann/mm10.CGI.bed | awk '$5>0{print $1,$2,$3,$5}' | sort -k1,1 -k2n,2 > $wdir/mm10.CGI.bg
bedGraphToBigWig $wdir/mm10.CGI.bg ~/ann/mm10.chromSizes $wdir/mm10.CGI.bw

###############NF-ICM K4 K9 K27
ddir=~/workspace/9.NT-ChIP/2.public/3.bowtie2/c.merge_rep
wdir=~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/6.cmm/9.NF-ICM_HM
for i in $ddir/ICM.K*bam; do k=${i##*ICM.}; k=${k%%.*}
nohup bamCompare -b1 $i -b2 $ddir/ICM.input.sorted.bam -o $wdir/NF-$k.bw -of bigwig \
--scaleFactorsMethod None --operation ratio --pseudocount 1 \
--binSize 10000 -p 16 --normalizeUsing RPKM --ignoreDuplicates --samFlagInclude 2 &
done

##############################
#从第5个开始，全部统一做了
##############################
nt=0.1
wdir=~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/6.cmm/a.all
ddir=~/workspace/8.NT-HiC/f.IS_ALL/2.except_res40k_is1M_ids200k_nt01/1.cat_boundary
awk -v nt=$nt '$5>nt{print $1,$2-1,$3}' $ddir/cc.bed > $wdir/CC.bed
awk -v nt=$nt '$5>nt{print $1,$2-1,$3}' $ddir/icm.bed > $wdir/NT-ICM.bed
awk -v nt=$nt '$5>nt{print $1,$2+40000-1,$3-40000}' ~/workspace/8.NT-HiC/3.public_data/f.XW_analysis/4.IS/1.cat_boundary/ICM.bed > $wdir/NF-ICM.bed
awk -v nt=$nt '$5>nt{print $1,$2-1,$3}' $ddir/te.bed > $wdir/NT-TE.bed
intersectBed -a $wdir/NT-ICM.bed -b $wdir/NF-ICM.bed -wa -u > $wdir/NT-NF.bed
intersectBed -a $wdir/NT-ICM.bed -b $wdir/NF-ICM.bed -wa -v > $wdir/NT-spe.bed


for i in $wdir/*bw; do k1=$(basename $i .bw)
for j in $wdir/*bed; do k2=$(basename $j .bed)
computeMatrix reference-point -S $i -R $j -p 16 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/${k2}.${k1}.gz -a 500000 -b 500000 --referencePoint center 
done
done &
wait
Rscript $wdir/../plot.R &
wait

#nt=0.25
nt=0.25
wdir=~/workspace/8.NT-HiC/f.IS_ALL/b.map_to_boundary/6.cmm/b.all
mkdir -p $wdir/0.pdfs $wdir/1.pngs
ddir=~/workspace/8.NT-HiC/f.IS_ALL/2.except_res40k_is1M_ids200k_nt01/1.cat_boundary
awk -v nt=$nt '$5>nt{print $1,$2-1-40000,$3+40000}' $ddir/cc.bed > $wdir/CC.bed
awk -v nt=$nt '$5>nt{print $1,$2-1-40000,$3+40000}' $ddir/icm.bed > $wdir/NT-ICM.bed
awk -v nt=$nt '$5>nt{print $1,$2-1,$3}' ~/workspace/8.NT-HiC/3.public_data/f.XW_analysis/4.IS/1.cat_boundary/ICM.bed > $wdir/NF-ICM.bed
awk -v nt=$nt '$5>nt{print $1,$2-1-40000,$3+40000}' $ddir/te.bed > $wdir/NT-TE.bed
awk -v nt=$nt '$5>nt{print $1,$2-1,$3}' ~/workspace/8.NT-HiC/3.public_data/f.XW_analysis/4.IS/1.cat_boundary/Sperm.bed > $wdir/Sperm.bed
intersectBed -a $wdir/NT-ICM.bed -b $wdir/NF-ICM.bed -wa -u -r -f 0.1 > $wdir/NT-NF.bed
intersectBed -a $wdir/NT-ICM.bed -b $wdir/NF-ICM.bed -wa -v -r -f 0.1 > $wdir/NT-spe.bed
intersectBed -b $wdir/NT-ICM.bed -a $wdir/NF-ICM.bed -wa -v -r -f 0.1 > $wdir/NF-spe.bed
#wc -l *bed #Figure4A，NT-ICM和NF-ICM的boundary的交集和特有的数量
intersectBed -b $wdir/CC.bed -a $wdir/NT-spe.bed -wa -u | wc -l #148
intersectBed -b $wdir/CC.bed -a $wdir/NT-spe.bed -wa -v | wc -l #115
intersectBed -b $wdir/CC.bed -a $wdir/NF-spe.bed -wa -u | wc -l #259
intersectBed -b $wdir/CC.bed -a $wdir/NF-spe.bed -wa -v | wc -l #241
#NT和NF特异的boundary从CC继承下来的差别，p=0.25

intersectBed -b $wdir/Sperm.bed -a $wdir/NT-spe.bed -wa -u | wc -l #53
intersectBed -b $wdir/Sperm.bed -a $wdir/NT-spe.bed -wa -v | wc -l #210
#NT特异的boundary从CC和Sperm继承下来的差别,p=2.2e-16



for i in $wdir/../a.all/*bw; do k1=$(basename $i .bw)
for j in $wdir/*bed; do k2=$(basename $j .bed)
computeMatrix reference-point -S $i -R $j -p 16 -bs 10000 --missingDataAsZero --averageTypeBins mean -o $wdir/${k2}.${k1}.gz -a 500000 -b 500000 --referencePoint center 
done &
done 

wait
Rscript $wdir/../plot.R $wdir &
wait

###test
intersectBed -a CC.bed -b NF-ICM.bed -wa -v -r -f 0.1 > CC-spe.bed
intersectBed -b CC.bed -a NF-ICM.bed -wa -v -r -f 0.1 > NF-spe.bed
intersectBed -b NT-ICM.bed -a CC-spe.bed -wa -u -r -f 0.1 | wc -l #148
intersectBed -b NT-ICM.bed -a CC-spe.bed -wa -v -r -f 0.1 | wc -l #115
intersectBed -b NT-ICM.bed -a NF-spe.bed -wa -u -r -f 0.1 | wc -l #259
intersectBed -b NT-ICM.bed -a NF-spe.bed -wa -v -r -f 0.1 | wc -l #241