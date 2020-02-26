[~/workspace/8.NT-HiC/l.IS_analysis/1.cc_2cell_20180703]$ln -s ~/workspace/8.NT-HiC/f.IS_ALL/4.except_norm_mean_raw_20180531/cat_boundary/cc.bed ./
[~/workspace/8.NT-HiC/l.IS_analysis/1.cc_2cell_20180703]$ln -s ~/workspace/8.NT-HiC/f.IS_ALL/4.except_norm_mean_raw_20180531/cat_boundary/e2cell.bed ./
[~/workspace/8.NT-HiC/l.IS_analysis/1.cc_2cell_20180703]$ln -s ~/workspace/8.NT-HiC/3.public_data/f.XW_analysis/4.IS/cat_boundary/early_2cell.bed ./

bedtools intersect -a cc.bed -b e2cell.bed -r -e -f 0.5 -wa | bedtools intersect -a - -b early_2cell.bed -r -e -f 0.5 -wa -v > temp1.bed
bedtools intersect -a early_2cell.bed -b cc.bed -r -e -f 0.5 -wa -v | bedtools intersect -a - -b e2cell.bed -r -e -f 0.5 -wa -v > temp2.bed


bedtools intersect -a early_2cell.bed -b cc.bed -r -e -f 0.5 -v | bedtools intersect -a - -b e2cell.bed -r -e -f 0.5 -v > temp1.bed
bedtools intersect -a early_2cell.bed -b cc.bed -r -e -f 0.5 -wa | bedtools intersect -a - -b e2cell.bed -r -e -f 0.5 -v > temp2.bed
bedtools intersect -a early_2cell.bed -b cc.bed -r -e -f 0.5 -wa | bedtools intersect -a - -b e2cell.bed -r -e -f 0.5 -wa > temp3.bed

bedtools intersect -a e2cell.bed -b cc.bed -r -e -f 0.5 -v | bedtools intersect -a - -b early_2cell.bed -r -e -f 0.5 -v > temp4.bed
bedtools intersect -a e2cell.bed -b cc.bed -r -e -f 0.5 -wa | bedtools intersect -a - -b early_2cell.bed -r -e -f 0.5 -v > temp5.bed
bedtools intersect -a e2cell.bed -b early_2cell.bed -r -e -f 0.5 -wa | bedtools intersect -a - -b cc.bed -r -e -f 0.5 -v > temp6.bed
#280kb

#####plot Profile of K9me3
k9me3=~/workspace/9.NT-ChIP/2.public/3.bowtie2/c.merge_rep/bamCoverage/2cell.K9me3.bw

n=1
computeMatrix scale-regions -S $k9me3 -R cc.bed e2cell.bed early_2cell.bed temp1.bed temp2.bed temp3.bed temp4.bed temp5.bed temp6.bed \
-p 64 -m 280000 -bs 4000 --missingDataAsZero --averageTypeBins mean -o 1.plotK9me3/$n.gz --startLabel "" --endLabel "" -a 200000 -b 200000 
plotHeatmap -m 1.plotK9me3/$n.gz -o 1.plotK9me3/$n.h.png --dpi 360 --regionsLabel CC NT NF "NF-CC-NT" "NF+CC-NT" "NF+CC+NT" "NT-CC-NF" "NT+CC-NF" "NT+NF-CC" \
-T "" --startLabel "" --endLabel "" --colorMap "bwr" --perGroup

for i in [!t]*bed; do k=${i%%.*}
annotatePeaks.pl $i mm10 -annStats 2.annotatePeaks/$k.txt > 2.annotatePeaks/$k.tab &
done