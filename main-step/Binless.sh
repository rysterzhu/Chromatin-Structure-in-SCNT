read name / chr_reads1 / pos_reads1 / strand_reads1 / chr_reads2 / pos_reads2 / strand_reads2 / fragment_size [/ allele_specific_tag]

id: ID of the mapped pair (alphanumeric without spaces)
chr1: chromosome of read 1 (alphanumeric without spaces)
pos1: position of the 5 end of read 1 (integer)
strand1: whether the read maps on the forward (1) or reverse (0) direction
length1: how many bases were mapped (integer)
re.up1: upstream restriction site closest to pos1 (integer <= pos1)
re.dn1: downstream restriction site closest to pos1 (integer > pos1)
chr2
pos2
strand2
length2
re.up2
re.dn2

E00572:82:HCNJ7CCXY:7:1108:14458:35326	chr1	3000761	+	chr1	15538042	+	216	HIC_chr1_2	HIC_chr1_33255	42	42
E00572:82:HCNJ7CCXY:7:1111:20121:30685	chr1	3000788	+	chr9	112342636	-	313	HIC_chr1_2	HIC_chr9_273451	35	42
E00572:82:HCNJ7CCXY:7:2114:5000:68799	chr1	3000837	-	chr1	17192496	+	269	HIC_chr1_3	HIC_chr1_37650	42	23
E00572:82:HCNJ7CCXY:7:1102:11688:44345	chr1	3000867	-	chr9	115962456	-	354	HIC_chr1_3	HIC_chr9_282167	38	42

#convert all valid pairs to TADbit tsv
ddir=~/workspace/8.NT-HiC/2.merge_sample/2.replicate/1.allValidPairs
wdir=~/workspace/8.NT-HiC/8.softeware_test/i.binless/1.replicates
mkd $wdir
for i in $ddir/*; do k=$(basename $i _allValidPairs)
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,int($8/2),a[$9],b[$9],$5,$6,$7,int($8/2),a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $i> $wdir/$k.tsv &
done 


ddir=~/workspace/8.NT-HiC/2.merge_sample/4.except/1.allValidPairs
wdir=~/workspace/8.NT-HiC/8.softeware_test/i.binless/1.replicates
mkd $wdir
for i in $ddir/*; do k=$(basename $i _allValidPairs)
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,int($8/2),a[$9],b[$9],$5,$6,$7,int($8/2),a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $i> $wdir/$k.tsv &
done 


ddir=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/1.allValidPairs
wdir=~/workspace/8.NT-HiC/8.softeware_test/i.binless/3.Xiewei
mkd $wdir
for i in $ddir/*; do k=$(basename $i _allValidPairs)
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,int($8/2),a[$9],b[$9],$5,$6,$7,int($8/2),a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $i> $wdir/$k.tsv &
done 


awk '$2=="chr7"&&$8=="chr7"{print}' 3.Xiewei/early_2cell.tsv > 4.e2cell/NF-chr7.tsv &
awk '$2=="chr7"&&$8=="chr7"{print}' 2.samples/e2cell.tsv > 4.e2cell/NT-chr7.tsv &



###############################test kdm4d
#一种是allValidPairs不分重复
#一种是allValidPairs分重复
#一种是bwt2Pairs分重复
ddir=~/workspace/8.NT-HiC/1.align/d.bwt2Pairs
wdir=~/workspace/8.NT-HiC/8.softeware_test/i.binless/5.test


awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr7"&&$5=="chr7"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed ~/workspace/8.NT-HiC/2.merge_sample/3.sample/1.allValidPairs/kdm4d-e2cell_allValidPairs > $wdir/kdm4d-chr7.tsv

awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr7"&&$5=="chr7"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed ~/workspace/8.NT-HiC/2.merge_sample/2.replicate/1.allValidPairs/kdm4d-e2cell_rep1_allValidPairs > $wdir/kdm4d_rep1-chr7.tsv &
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr7"&&$5=="chr7"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed ~/workspace/8.NT-HiC/2.merge_sample/2.replicate/1.allValidPairs/kdm4d-e2cell_rep2_allValidPairs > $wdir/kdm4d_rep2-chr7.tsv &


awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr7"&&$5=="chr7"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir/kdm4d-e2cell_rep1* > $wdir/kdm4d_rep1-bwt2-chr7.tsv &
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr7"&&$5=="chr7"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir/kdm4d-e2cell_rep2* > $wdir/kdm4d_rep2-bwt2-chr7.tsv &

####################################
ddir1=~/workspace/8.NT-HiC/1.align/a.all_validPairs
ddir2=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/0.validPairs
wdir=~/workspace/8.NT-HiC/8.softeware_test/i.binless/1.e2cell-replicates-chr7

awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr7"&&$5=="chr7"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir2/early_2cell_rep[123]* > $wdir/NF-rep6.tsv &
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr7"&&$5=="chr7"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir2/early_2cell_rep4* > $wdir/NF-rep4.tsv &
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr7"&&$5=="chr7"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir2/early_2cell_rep5* > $wdir/NF-rep5.tsv &

awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr7"&&$5=="chr7"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/e2cell_rep3* > $wdir/NT-rep3.tsv &
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr7"&&$5=="chr7"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/e2cell_rep4* > $wdir/NT-rep4.tsv &

awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr7"&&$5=="chr7"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/kdm4d-e2cell_rep1* > $wdir/kdm4d-rep1.tsv &
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr7"&&$5=="chr7"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/kdm4d-e2cell_rep2* > $wdir/kdm4d-rep2.tsv &

awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr7"&&$5=="chr7"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/partheno-e2cell_rep1* > $wdir/partheno-rep1.tsv &
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr7"&&$5=="chr7"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/partheno-e2cell_rep2* > $wdir/partheno-rep2.tsv &



ddir1=~/workspace/8.NT-HiC/2.merge_sample/3.sample/1.allValidPairs
ddir2=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/1.allValidPairs
wdir=~/workspace/8.NT-HiC/p.loops/a.binless/2.e2cell-sample-chr7
mkd $wdir
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr7"&&$5=="chr7"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir2/early_2cell* > $wdir/NF.tsv &
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr7"&&$5=="chr7"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/e2cell* > $wdir/NT.tsv &
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr7"&&$5=="chr7"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/kdm4d-e2cell* > $wdir/kdm4d.tsv &
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr7"&&$5=="chr7"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/cc* > $wdir/cc.tsv &
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr7"&&$5=="chr7"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/tsa* > $wdir/tsa.tsv &

ddir1=~/workspace/8.NT-HiC/2.merge_sample/3.sample/1.allValidPairs
ddir2=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/1.allValidPairs
wdir=~/workspace/8.NT-HiC/p.loops/a.binless/3.e2cell-sample-chr16
mkd $wdir
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr16"&&$5=="chr16"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir2/early_2cell* > $wdir/NF.tsv &
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr16"&&$5=="chr16"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/e2cell* > $wdir/NT.tsv &
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr16"&&$5=="chr16"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/kdm4d-e2cell* > $wdir/kdm4d.tsv &
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr16"&&$5=="chr16"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/cc* > $wdir/cc.tsv &
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr16"&&$5=="chr16"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/tsa* > $wdir/tsa.tsv &

#########################################################################
ddir1=~/workspace/8.NT-HiC/1.align/a.all_validPairs
ddir2=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/0.validPairs
wdir=~/workspace/8.NT-HiC/p.loops/a.binless/7.e2cell-replicates-chr16

awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr16"&&$5=="chr16"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir2/early_2cell_rep[123]* > $wdir/NF-rep6.tsv &
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr16"&&$5=="chr16"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir2/early_2cell_rep4* > $wdir/NF-rep4.tsv &
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr16"&&$5=="chr16"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir2/early_2cell_rep5* > $wdir/NF-rep5.tsv &

awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr16"&&$5=="chr16"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/cc_rep1* > $wdir/cc-rep1.tsv &
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr16"&&$5=="chr16"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/cc_rep2* > $wdir/cc-rep2.tsv &
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr16"&&$5=="chr16"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/cc_rep3* > $wdir/cc-rep3.tsv &
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr16"&&$5=="chr16"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/cc_rep4* > $wdir/cc-rep4.tsv &


awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr16"&&$5=="chr16"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/kdm4d-e2cell_rep1* > $wdir/kdm4d-rep1.tsv &
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr16"&&$5=="chr16"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/kdm4d-e2cell_rep2* > $wdir/kdm4d-rep2.tsv &

awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr16"&&$5=="chr16"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/tsa-e2cell_rep1* > $wdir/tsa-rep1.tsv &
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr16"&&$5=="chr16"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/tsa-e2cell_rep2* > $wdir/tsa-rep2.tsv &

awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr16"&&$5=="chr16"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/e2cell_rep3* > $wdir/NT-rep3.tsv &
awk 'NR==FNR{a[$4]=$2;b[$4]=$3} NR>FNR&&$2=="chr16"&&$5=="chr16"{if($4=="+"){$4=1}else{$4=0};if($7=="+"){$7=1}else{$7=0}
print $1,$2,$3,$4,107,a[$9],b[$9],$5,$6,$7,107,a[$10],b[$10]}' $ConfigHP/MboI_resfrag_mm10.bed $ddir1/e2cell_rep4* > $wdir/NT-rep4.tsv &


