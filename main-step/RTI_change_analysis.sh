ddir=~/workspace/9.NT-ChIP/3.deeptools/1.bamCoverage/h.20180909
wdir=~/workspace/8.NT-HiC/g.DI_ALL/i.RTI_change_analysis
n=all.RTI
multiBigwigSummary BED-file --BED $wdir/$n.bed -b $ddir/cc_input_rep1.bw $ddir/cc_K9me3_rep1.bw \
$ddir/cc_input_rep2.bw $ddir/cc_K9me3_rep2.bw -o $wdir/$n.npz --smartLabels -p 32 --outRawCounts $wdir/$n
awk 'FNR>1{printf("%s\t%s\t%s\t%.6f\t%.6f\n",$1,$2,$3,$5/$4,$7/$6)}' $n > k9me3-all.tab #ratio
#

######################CC TAD
ddir=~/workspace/9.NT-ChIP/3.deeptools/1.bamCoverage/h.20180909
wdir=~/workspace/8.NT-HiC/g.DI_ALL/i.RTI_change_analysis

ln -s /home/qszhu/workspace/8.NT-HiC/g.DI_ALL/2.except_norm_100M_20180825/2.unoverlap_tads/cc.tad CC.tad

n=CC
multiBigwigSummary BED-file --BED $wdir/$n.tad -b $ddir/cc_input_rep1.bw $ddir/cc_K9me3_rep1.bw \
$ddir/cc_input_rep2.bw $ddir/cc_K9me3_rep2.bw -o $wdir/$n.npz --smartLabels -p 32 --outRawCounts $wdir/$n
awk 'FNR>1{printf("%s\t%s\t%s\t%.6f\t%.6f\n",$1,$2,$3,$5/$4,$7/$6)}' $n > k9me3-CC-TAD.tab #ratio

~/workspace/8.NT-HiC/g.DI_ALL/i.RTI_change_analysis/genes:
intersectBed -a ~/ann/mm10_promoter_1000-1000.bed -b cluster.bed -r -e -f 0.5 -wo | cut -f 5,10 | sort -k2n,2 -k1,1 |uniq> temp.tab

mv k9me3-CC-TAD.tab  k9me3-CC-TAD2.tab 

#############用multiBamSummary
wdir=~/workspace/8.NT-HiC/g.DI_ALL/i.RTI_change_analysis/1.multiBamSummary
ddir1=~/workspace/9.NT-ChIP/1.align/0.input_usable
ddir2=~/workspace/9.NT-ChIP/1.align/1.K9me3_usable
ln -s /home/qszhu/workspace/8.NT-HiC/g.DI_ALL/2.except_norm_100M_20180825/2.unoverlap_tads/cc.tad CC.tad
n=CC
multiBamSummary BED-file --BED $wdir/$n.tad -b $ddir1/cc_input_rep1.sorted.bam $ddir2/cc_K9me3_rep1.sorted.bam $ddir1/cc_input_rep2.sorted.bam $ddir2/cc_K9me3_rep2.sorted.bam \
-o $wdir/$n.npz --smartLabels -p 32 --outRawCounts $wdir/$n --ignoreDuplicates --samFlagInclude 2 

awk -v cc_input_rep1=60915922 -v cc_input_rep2=58549002 -v cc_K9me3_rep1=63109442 -v cc_K9me3_rep2=67676976 'NR==1{print}
NR>1{l=($3-$2)/1000
a=($4*1e6)/(l*cc_input_rep1);b=($5*1e6)/(l*cc_K9me3_rep1);c=($6*1e6)/(l*cc_input_rep2);d=($7*1e6)/(l*cc_K9me3_rep2);
printf("%s\t%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f\n",$1,$2,$3,a,b,c,d)
}' $n > ../k9me3-CC-TAD.tab

awk '{print $1,1,$2}' ~/ann/hic-pro/chrom_mm10.sizes > test.tad
n=test
multiBamSummary BED-file --BED $wdir/$n.tad -b $ddir1/cc_input_rep1.sorted.bam $ddir2/cc_K9me3_rep1.sorted.bam $ddir1/cc_input_rep2.sorted.bam $ddir2/cc_K9me3_rep2.sorted.bam \
-o $wdir/$n.npz --smartLabels -p 32 --outRawCounts $wdir/$n --ignoreDuplicates --samFlagInclude 2 
awk 'NR==FNR{A+=$4;B+=$5;C+=$6;D+=$7} 
NR>FNR&&FNR==1{printf("chr\tstart\tend\tcc_input_rep1\tcc_K9me3_rep1\tcc_input_rep2\tcc_K9me3_rep2\tcc_input\tcc_K9me3\n")}
NR>FNR&&FNR>1{l=($3-$2)/1000
a=($4*1e6)/(l*A);b=($5*1e6)/(l*B);c=($6*1e6)/(l*C);d=($7*1e6)/(l*D);
e=(($4+$6)*1e6)/(l*(A+C));f=(($5+$7)*1e6)/(l*(B+D))
printf("%s\t%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",$1,$2,$3,a,b,c,d,e,f)
}' test CC > ../k9me3-CC-TAD.tab


awk -v cc_input_rep1=41619996 -v cc_input_rep2=33806301 -v cc_K9me3_rep1=42891472 -v cc_K9me3_rep2=67676976 'NR==1{print}
NR>1{l=($3-$2)/1000
a=($4*1e6)/(l*cc_input_rep1);b=($5*1e6)/(l*cc_K9me3_rep1);c=($6*1e6)/(l*cc_input_rep2);d=($7*1e6)/(l*cc_K9me3_rep2);
printf("%s\t%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f\n",$1,$2,$3,a,b,c,d)
}' $n > ../k9me3-CC-TAD.tab


#################kdm4d-2cell RTI in CC K9me3 mark TAD
ddir1=~/workspace/8.NT-HiC/5.maps/2.replicate/7.ICED_norm
ddir2=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/7.replicate/8.norm_to_100M
wdir=~/workspace/8.NT-HiC/g.DI_ALL/i.RTI_change_analysis/2.kdm4d-2cell
for i in $ddir1/e2cell*_40000_iced.matrix $ddir1/kdm4d-e2cell*_40000_iced.matrix $ddir2/early_2cell*_40000_iced.matrix; do key=$(basename $i _40000_iced.matrix)
python ~/codes/score_of_tad.py -i $i -I ~/workspace/8.NT-HiC/g.DI_ALL/i.RTI_change_analysis/CC-K9_CC-TAD_flag3.bed -b $ConfigHP/40000_mm10.bed -o $wdir/${key}_RTI.txt -c -1 -r 40000 &
done;wait
awk 'NR==FNR{flag[$1"\t"$2]=$4} NR!=FNR&&FNR>1{gsub("(^.*/)|(_RTI.txt)","",FILENAME);print FILENAME,flag[$1"\t"$2],$0}' ~/workspace/8.NT-HiC/g.DI_ALL/i.RTI_change_analysis/CC-K9_CC-TAD_flag3.bed $wdir/*RTI.txt > $wdir/all.RTI

#merge-rep
ddir1=~/workspace/8.NT-HiC/5.maps/4.except/1.raw
ddir2=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/3.raw
wdir=~/workspace/8.NT-HiC/g.DI_ALL/i.RTI_change_analysis/2.kdm4d-2cell/1.merge-rep
tad=~/workspace/8.NT-HiC/g.DI_ALL/i.RTI_change_analysis/CC-K9_CC-TAD_flag3.bed
mkd $wdir
for i in $ddir1/cc_40000.matrix $ddir1/*e2cell_40000.matrix $ddir2/early_2cell*_40000.matrix; do k=$(basename $i _40000.matrix) 
python ~/codes/score_of_tad.py -i $i -I $tad -b $ConfigHP/40000_mm10.bed -o $wdir/${k}_RTI.txt -c -1 -r 40000 &
done;wait
awk 'NR==FNR{flag[$1"\t"$2]=$4} NR!=FNR&&FNR>1{gsub("(^.*/)|(_RTI.txt)","",FILENAME);print FILENAME,flag[$1"\t"$2],$0}' $tad $wdir/*RTI.txt > $wdir/all.RTI

#因数据更新，使用新的数据重做RTI in marked and unmarked K9me3 TAD
ddir=~/workspace/8.NT-HiC/g.DI_ALL/a.Consolidation_analysis/6.NF_20180803/2.CC
awk 'NR==FNR{flag[$1"\t"$2]=$4} NR!=FNR&&FNR>1{gsub("(^.*/)|(_RTI.txt)","",FILENAME);print FILENAME,flag[$1"\t"$2],$0}' $tad $wdir/cc_RTI.txt $wdir/e2cell_RTI.txt $ddir/A-6h_CC_RTI.txt $ddir/A-12h_CC_RTI.txt > $wdir/NT.RTI



