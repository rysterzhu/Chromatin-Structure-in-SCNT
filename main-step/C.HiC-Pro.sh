##2018年1月19日 05h/  1h/  2h/  4cell/  6h/  8cell/  cc/  e2cell/  morula/
#########
ddir=~/workspace/8.NT-HiC/0.cutadapt-qc/4.fastp/0.all
wdir=~/workspace/8.NT-HiC/c.HiC-Pro
for i in $ddir/*_rep1_seq1.R1.fastq.gz;do i=${i##*/};mkdir $wdir/0.input/${i%%_*}; done
for i in $ddir/*_seq1.R1.fastq.gz;do i=${i##*/};mkdir $wdir/0.input/${i%%_*}/${i%%_seq*}; done

for i in $ddir/*gz; do o=${i##*/}; ln -s $i $wdir/0.input/${o%%_*}/${o%%_seq*}/$o; done

#p4 每个sample 8个线程
for i in $wdir/0.input/*; do s=${i##*/}
nohup HiC-Pro -i $i -o $wdir/$s -c $ConfigHP/config-hicpro-p4.txt > $wdir/logs/nohup.${s}.log &
done
#太慢了

ln -s $wdir/*/hic_results/matrix/*/*/*/*_iced.matrix $wdir/3.matrix/
ln -s $wdir/*/hic_results/data/*/*allValidPairs $wdir/2.allValidPairs/

###################################统计日志到tab
cp */hic_results/data/*/*stat stats/
cp */bowtie_results/bwt2/*/*stat stats/
cd stats
echo -e "Sample\t\
Total_pairs_processed\t\
Reported_pairs\t%map_pairs\t\
valid_interaction\t%valid_interaction\t\
valid_interaction_rmdup\t%Duplicates\t\
%final_valid\t\
trans_interaction\t%trans\t\
cis_interaction\t%cis\t\
cis_shortRange\tcis_longRange" > mergestats

for i in *mpairstat; do 
k=${i%%.*}
j=${k}_allValidPairs.mergestat
awk -v k=$k 'BEGIN{}
$1=="Total_pairs_processed"{a=($2/1e6)}
$1=="Reported_pairs"{b=($2/1e6)}
$1=="valid_interaction"{c=($2/1e6)}
$1=="valid_interaction_rmdup"{d=($2/1e6)}
$1=="trans_interaction"{e=($2/1e6)}
$1=="cis_interaction"{f=($2/1e6)}
$1=="cis_shortRange"{g=($2/1e6)}
$1=="cis_longRange"{h=($2/1e6)}
END{printf("%s\t%.2f\t%.2f\t%.2f%\t%.2f\t%.2f%\t%.2f\t%.2f%\t%.2f%\t%.2f\t%.2f%\t%.2f\t%.2f%\t%.2f\t%.2f\n",k,a,b,100*b/a,c,100*c/b,d,(c-d)*100/c,100*d/a, e,100*e/d,f,100*f/d,g,h)}
' $i $j >> mergestats
done

cut -f 1,7,10- mergestats | awk -v OFS=" | " -v FS="\t" 'NR==1{$1=$1;print $0;for(i=1;i<=NF;i++){$i="--"};print $0;} NR>1{$1=$1;print $0}' | xsel -ib
awk -v OFS=" | " -v FS="\t" 'NR==1{$1=$1;print $0;for(i=1;i<=NF;i++){$i="--"};print $0;} NR>1{$1=$1;print $0}' mergestats | xsel -ib

#############################################################################################################################
########################################################all repicates PCA 2018年1月24日
ddir=~/workspace/8.NT-HiC/c.HiC-Pro/3.matrix
wdir=~/workspace/8.NT-HiC/9.PCA/3.all_replicates
#取100k分辨率的matrix，分出interchrom和intrachrom
for i in $ddir/*100000_iced.matrix; do o=${i##*/};#o=${o/_100000_iced.matrix}
awk -v out=$o -v wdir=$wdir 'NR==FNR{if($4>a[$1]){a[$1]=$4}} NR>FNR{for(i in a){if($1<=a[i]&&$2>a[i]){print $0 > wdir"/inter."out;next}};
print $0 > wdir"/intra."out}' $ConfigHP/100000_mm10.bed $i &
done
#取intrachrom的5M内的互作
for i in $wdir/intra*matrix; do o=${i##*intra.}; o=${o/_100000_iced.matrix/.txt};
awk '$2-$1<=50&&$2-$1>0{print $1"-"$2,$3}' $i > $wdir/$o &
done
#normalize 所有regions的范围是一定的，所以只需要除以总的interaction和即可
for i in $wdir/*.txt; do o=${i##*/}
awk 'NR==FNR{a+=$3} NR!=FNR{$2=$2*1e7/a;print }' $ddir/${i/.txt/_100000_iced.matrix} $i > $wdir/${o/.txt/.norm} & done





##############################################################################################
##########2018年3月22日   
#############################################################################################
ddir=~/workspace/8.NT-HiC/0.cutadapt-qc/4.fastp/9.nt_6h_l2_4cell_20180233
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/a.20180322

for i in $ddir/*_rep1_seq1.R1.fastq.gz;do i=${i##*/};mkdir $wdir/0.input/${i%%_*}; done
for i in $ddir/*_seq1.R1.fastq.gz;do i=${i##*/};mkdir $wdir/0.input/${i%%_*}/${i%%_seq*}; done

for i in $ddir/*gz; do o=${i##*/}; ln -s $i $wdir/0.input/${o%%_*}/${o%%_seq*}/$o; done

for i in $wdir/0.input/*; do s=${i##*/}
nohup HiC-Pro -i $i -o $wdir/$s -c $ConfigHP/config-hicpro.txt > $wdir/logs/nohup.${s}.log &
done


ln -s $wdir/*/hic_results/matrix/*/*/*/*_iced.matrix $wdir/3.matrix/
ln -s $wdir/*/hic_results/data/*/*allValidPairs $wdir/2.allValidPairs/


#####library twice test each seq
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/a.20180322/a.library_twice_test
nohup HiC-Pro -s proc_hic -s quality_checks -i $wdir/0.input/bwt2 -o $wdir -c $ConfigHP/config-hicpro.txt >> $wdir/nohup.20180323.log &


##################2018年4月5日 norm to depth
bash ~/codes/hic_norm_to_depth.sh -r 100000 -d ~/workspace/8.NT-HiC/c.HiC-Pro/3.matrix/ -w ~/workspace/8.NT-HiC/c.HiC-Pro/4.norm_to_depth/
bash ~/codes/hic_norm_to_depth.sh -r 100000 -d ~/workspace/8.NT-HiC/c.HiC-Pro/a.20180322/3.matrix/ -w ~/workspace/8.NT-HiC/c.HiC-Pro/a.20180322/4.norm_to_depth/



##############################################################################################
##########2018年4月10日
#############################################################################################
ddir=~/workspace/8.NT-HiC/0.cutadapt-qc/4.fastp/10.6h_4cell_morula_20180409
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/b.20180410

for i in $ddir/*_rep1_seq1.R1.fastq.gz;do i=${i##*/};mkdir -p $wdir/0.input/${i%%_*}; done
for i in $ddir/*_seq1.R1.fastq.gz;do i=${i##*/};mkdir -p $wdir/0.input/${i%%_*}/${i%%_seq*}; done
#每个rep分开做HiC-Pro，各自目录下有不同的seq；
for i in $ddir/*gz; do o=${i##*/}; ln -s $i $wdir/0.input/${o%%_*}/${o%%_seq*}/$o; done

mkdir -p $wdir/logs
for i in $wdir/0.input/*; do s=${i##*/}
nohup HiC-Pro -i $i -o $wdir/$s -c $ConfigHP/config-hicpro-p8.txt > $wdir/logs/nohup.${s}.log &
done

mkdir -p $wdir/3.matrix/ $wdir/2.allValidPairs/
ln -s $wdir/*/hic_results/matrix/*/*/*/*_iced.matrix $wdir/3.matrix/
ln -s $wdir/*/hic_results/data/*/*allValidPairs $wdir/2.allValidPairs/

#statistic logs
bash ~/codes/hic_multiqc.sh -w $wdir
awk -v OFS=" | " -v FS="\t" 'NR==1{$1=$1;print $0;for(i=1;i<=NF;i++){$i="--"};print $0;} NR>1{$1=$1;print $0}' mergestats | xsel -ib

#norm to depth for PCA
bash ~/codes/hic_norm_to_depth.sh -r 100000 -d $wdir/3.matrix/ -w $wdir/4.norm_to_depth

#prepare for PCA
bash ~/codes/hic_pca.sh -d $wdir/4.norm_to_depth -w ~/workspace/8.NT-HiC/9.PCA/5.norm_to_depth_allRep -r 100000

##############################################################################################
##########2018年4月26日
#############################################################################################
ddir=~/workspace/8.NT-HiC/0.cutadapt-qc/4.fastp/11.4cell_morula_20180426
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/c.20180426

#for i in $ddir/*_rep1_seq1.R1.fastq.gz;do i=${i##*/};mkdir -p $wdir/0.input/${i%%_*}; done
for i in $ddir/*_seq2.R1.fastq.gz;do i=${i##*/};mkdir -p $wdir/0.input/${i%%_*}/${i%%_seq*}; done
#每个rep分开做HiC-Pro，各自目录下有不同的seq；
for i in $ddir/*gz; do o=${i##*/}; ln -s $i $wdir/0.input/${o%%_*}/${o%%_seq*}/$o; done

mkdir -p $wdir/logs
for i in $wdir/0.input/*; do s=${i##*/}
nohup HiC-Pro -i $i -o $wdir/$s -c $ConfigHP/config-hicpro.txt > $wdir/logs/nohup.${s}.log &
done

########merge for rep
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/1.merge_sample/

mkdir -p $wdir/1.validPairs/4cell_rep3
mkdir -p $wdir/1.validPairs/morula_rep4

ln -s ~/workspace/8.NT-HiC/1.align/*/1.validPairs/4cell_rep3*.validPairs $wdir/1.validPairs/4cell_rep3
ln -s ~/workspace/8.NT-HiC/1.align/*/1.validPairs/morula_rep4*.validPairs $wdir/1.validPairs/morula_rep4

(nohup HiC-Pro -i $wdir/1.validPairs -o $wdir -c $ConfigHP/config-hicpro.txt -s merge_persample >> $wdir/nohup.log
nohup HiC-Pro -i $wdir/hic_results/data -o $wdir -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.log )&

#######6h 2018年5月2日
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/1.merge_sample/
rm -r $wdir/1.validPairs/*
mkdir -p $wdir/1.validPairs/6h
ln -s ~/workspace/8.NT-HiC/1.align/*/1.validPairs/6h*.validPairs $wdir/1.validPairs/6h
(nohup HiC-Pro -i $wdir/1.validPairs -o $wdir -c $ConfigHP/config-hicpro.txt -s merge_persample >> $wdir/nohup.log
nohup HiC-Pro -i $wdir/hic_results/data -o $wdir -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.log )&


##################################4cell_rep2_seq3 补测二次建库的第一次(seq1)
ddir=~/workspace/8.NT-HiC/0.Pre/4cell/2.cutadapt
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/d.20180501
mkdir -p $wdir/0.input/4cell_rep2_seq3
ln -s $ddir/4cell_rep2_seq3*gz $wdir/0.input/4cell_rep2_seq3
nohup HiC-Pro -i $wdir/0.input -o $wdir/4cell_rep2_seq3 -c $ConfigHP/config-hicpro.txt > $wdir/nohup.log &

ddir=~/workspace/8.NT-HiC/1.align/4cell/1.validPairs
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/d.20180501
key=4cell_rep2
mkdir -p $wdir/1.validPairs/$key
ln -s $ddir/${key}* $wdir/1.validPairs/$key
(nohup HiC-Pro -i $wdir/1.validPairs -o $wdir/$key -c $ConfigHP/config-hicpro.txt -s merge_persample >> $wdir/nohup.log
nohup HiC-Pro -i $wdir/$key/hic_results/data -o $wdir/$key -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.log )&

key=4cell
mkdir -p $wdir/2.samplePairs/$key
ln -s $ddir/${key}* $wdir/2.samplePairs/$key
(nohup HiC-Pro -i $wdir/2.samplePairs -o $wdir/$key -c $ConfigHP/config-hicpro.txt -s merge_persample >> $wdir/nohup.log
nohup HiC-Pro -i $wdir/$key/hic_results/data -o $wdir/$key -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.log )&

##############test 4cell_rep2 是第二次建库好，还是加测好； 分别合并seq1 seq2和seq1 seq3
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/1.merge_sample/
rm -r $wdir/1.validPairs/*
mkdir -p $wdir/1.validPairs/4cell_rep2_test2
mkdir -p $wdir/1.validPairs/4cell_rep2_test3
ln -s ~/workspace/8.NT-HiC/1.align/*/1.validPairs/4cell_rep2_seq[12]*.validPairs $wdir/1.validPairs/4cell_rep2_test2
ln -s ~/workspace/8.NT-HiC/1.align/*/1.validPairs/4cell_rep2_seq[13]*.validPairs $wdir/1.validPairs/4cell_rep2_test3
nohup HiC-Pro -i $wdir/1.validPairs -o $wdir -c $ConfigHP/config-hicpro.txt -s merge_persample >> $wdir/nohup.log &


############################2018-05-05
######3h 6h 12h e2cell l2cell 
# 12h_rep1_seq1
# 3h_rep1_seq1
# 6h_rep5_seq2
# e2cell_rep3_seq1
# l2cell_rep2_seq1

ddir=~/workspace/8.NT-HiC/0.Pre/n.date/1.20180505/2.cutadapt
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/e.20180505

for i in $ddir/*.R1.fastq.gz; do key=${i##*/};key=${key%%.*}
mkdir -p $wdir/0.input/$key
ln -s $i $wdir/0.input/$key
ln -s ${i/R1/R2} $wdir/0.input/$key
done
nohup HiC-Pro -i $wdir/0.input -o $wdir/1.sequence -c $ConfigHP/config-hicpro.txt > $wdir/nohup.20180506.log &


nohup HiC-Pro -i $wdir/1.sequence/bowtie_results/bwt2 -o $wdir/1.sequence -c $ConfigHP/config-hicpro.txt -s proc_hic quality_checks merge_persample build_contact_maps ice_norm > $wdir/nohup.20180511.log &

##程序因磁盘空间满了出错，从sort开始跑

for i in $ddir/*.R1.fastq.gz; do key=${i##*/};key=${key%_*}; 
mkdir -p $wdir/2.replicate_valid/$key; 
ln -s ~/workspace/8.NT-HiC/1.align/*/1.validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/2.replicate_valid/$key; 
done
nohup HiC-Pro -i $wdir/2.replicate_valid -o $wdir/3.replicate -c $ConfigHP/config-hicpro.txt -s merge_persample > $wdir/nohup.20180512.rep.log &

for i in $ddir/*.R1.fastq.gz; do 
key=${i##*/};key=${key%_rep*}; 
mkdir -p $wdir/4.sample_valids/$key; 
ln -s ~/workspace/8.NT-HiC/1.align/*/1.validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/4.sample_valids/$key; 
done
nohup HiC-Pro -i $wdir/1.sequence/hic_results/data -o $wdir/1.sequence -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.20180512.log &

nohup HiC-Pro -i $wdir/3.replicate/hic_results/data -o $wdir/3.replicate -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.20180512.rep.log &

mv $wdir/1.sequence/hic_results/data/*/*allValidPairs ~/workspace/8.NT-HiC/2.merge_sample/1.sequence/1.allValidPairs/
mv $wdir/1.sequence/hic_results/data/*/*mergestat ~/workspace/8.NT-HiC/2.merge_sample/1.sequence/2.stats/


for i in 1.sequence/hic_results/data/*; do k=${i##*/}; 
mv ~/workspace/8.NT-HiC/2.merge_sample/1.sequence/1.allValidPairs/${k}_allValidPairs $i/;
done
nohup HiC-Pro -i $wdir/1.sequence/hic_results/data -o $wdir/1.sequence -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.20180512.log &
nohup HiC-Pro -i $wdir/5.sample/hic_results/data -o $wdir/5.sample -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.20180512.sample.log &

fs=(1.sequence 3.replicate 5.sample)
gs=(1.sequence 2.replicate 3.sample)
for i in {1..2}; do f=${fs[i]};g=${gs[i]}
ln $wdir/$f/hic_results/data/*/*allValidPairs ~/workspace/8.NT-HiC/2.merge_sample/$g/1.allValidPairs/
ln $wdir/$f/hic_results/data/*/*mergestat ~/workspace/8.NT-HiC/2.merge_sample/$g/2.stats/
ln $wdir/$f/hic_results/matrix/*/raw/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/1.raw
ln $wdir/$f/hic_results/matrix/*/iced/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/2.iced
ln $wdir/$f/hic_results/matrix/*/iced/*/*biases ~/workspace/8.NT-HiC/5.maps/$g/3.biases
done
#建立硬链接，这样即使删除了原文件也可

###############4.except 6h_rep1 6h_rep2 6h_rep4(5-16补充)
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/e.20180505
mkdir -p $wdir/6.except_valids/6h
ln -s ~/workspace/8.NT-HiC/1.align/a.all_validPairs/6h_rep[35]*._mm10.bwt2pairs.validPairs $wdir/6.except_valids/6h
(nohup HiC-Pro -i $wdir/6.except_valids -o $wdir/7.except -c $ConfigHP/config-hicpro.txt -s merge_persample > $wdir/nohup.except.log
nohup HiC-Pro -i $wdir/7.except/hic_results/data -o $wdir/7.except -c $ConfigHP/config-hicpro.txt -s build_contact_maps -s ice_norm >> $wdir/nohup.except.log) &

f=7.except
g=4.except
ln -f $wdir/$f/hic_results/data/*/*allValidPairs ~/workspace/8.NT-HiC/2.merge_sample/$g/1.allValidPairs/
ln -f $wdir/$f/hic_results/data/*/*mergestat ~/workspace/8.NT-HiC/2.merge_sample/$g/2.stats/
ln -f $wdir/$f/hic_results/matrix/*/raw/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/1.raw
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/2.iced
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*biases ~/workspace/8.NT-HiC/5.maps/$g/3.biases

#最好是用软连接，这样可以看到原始文件位置
ln -s -f ~/workspace/8.NT-HiC/2.merge_sample/3.sample/1.allValidPairs/[!6m]* ~/workspace/8.NT-HiC/2.merge_sample/4.except/1.allValidPairs/
ln -s -f ~/workspace/8.NT-HiC/5.maps/3.sample/1.raw/[!6m]* ~/workspace/8.NT-HiC/5.maps/4.except/1.raw/
ln -s -f ~/workspace/8.NT-HiC/5.maps/3.sample/2.iced/[!6m]* ~/workspace/8.NT-HiC/5.maps/4.except/2.iced/
ln -s -f ~/workspace/8.NT-HiC/5.maps/3.sample/3.biases/[!6m]* ~/workspace/8.NT-HiC/5.maps/4.except/3.biases/


############################2018-05-05
######2h
ddir=~/workspace/8.NT-HiC/0.Pre/n.date/2.20180514/2.cutadapt
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/f.20180514
for i in $ddir/*.R1.fastq.gz; do key=${i##*/};key=${key%%.*}
mkdir -p $wdir/0.input/$key
ln -s $i $wdir/0.input/$key
ln -s ${i/R1/R2} $wdir/0.input/$key
done
nohup HiC-Pro -i $wdir/0.input -o $wdir/1.sequence -c $ConfigHP/config-hicpro.txt > $wdir/nohup.20180515.log &

ln $wdir/1.sequence/hic_results/data/*/*validPairs ~/workspace/8.NT-HiC/1.align/a.all_validPairs/
ln $wdir/1.sequence/bowtie_results/bwt2/*/*pairstat ~/workspace/8.NT-HiC/1.align/b.all_stats/

#only 1 sequence; pass
for i in $ddir/*.R1.fastq.gz; do key=${i##*/};key=${key%_*}; 
mkdir -p $wdir/2.replicate_valid/$key; 
ln -s ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/2.replicate_valid/$key; 
done
nohup HiC-Pro -i $wdir/2.replicate_valid -o $wdir/3.replicate -c $ConfigHP/config-hicpro.txt -s merge_persample > $wdir/nohup.rep.log &

for i in $ddir/*.R1.fastq.gz; do 
key=${i##*/};key=${key%_rep*}; 
mkdir -p $wdir/4.sample_valids/$key; 
ln -s ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/4.sample_valids/$key; 
done
(nohup HiC-Pro -i $wdir/4.sample_valids -o $wdir/5.sample -c $ConfigHP/config-hicpro.txt -s merge_persample > $wdir/nohup.sample.log
nohup HiC-Pro -i $wdir/5.sample/hic_results/data -o $wdir/5.sample -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.sample.log )&


####################################
fs=(1.sequence 1.sequence 5.sample)
gs=(1.sequence 2.replicate 3.sample)
for i in {1..2}; do f=${fs[i]};g=${gs[i]}
ln -f $wdir/$f/hic_results/data/*/*allValidPairs ~/workspace/8.NT-HiC/2.merge_sample/$g/1.allValidPairs/
ln -f $wdir/$f/hic_results/data/*/*mergestat ~/workspace/8.NT-HiC/2.merge_sample/$g/2.stats/
ln -f $wdir/$f/hic_results/matrix/*/raw/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/1.raw
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/2.iced
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*biases ~/workspace/8.NT-HiC/5.maps/$g/3.biases
done

#最好是用软连接，这样可以看到原始文件位置
# key=2h
# ln -s -f ~/workspace/8.NT-HiC/2.merge_sample/3.sample/1.allValidPairs/$key* ~/workspace/8.NT-HiC/2.merge_sample/4.except/1.allValidPairs/
# ln -s -f ~/workspace/8.NT-HiC/5.maps/3.sample/1.raw/$key* ~/workspace/8.NT-HiC/5.maps/4.except/1.raw/
# ln -s -f ~/workspace/8.NT-HiC/5.maps/3.sample/2.iced/$key* ~/workspace/8.NT-HiC/5.maps/4.except/2.iced/
# ln -s -f ~/workspace/8.NT-HiC/5.maps/3.sample/3.biases/$key* ~/workspace/8.NT-HiC/5.maps/4.except/3.biases/
#新加rep不需要这些，本来就有指向

############################2018-05-27   12h_rep1_seq2 1h_rep3_seq1 2h_rep3_seq2 icm_rep1_seq1 te_rep1_seq1
ddir=~/workspace/8.NT-HiC/0.Pre/n.date/3.20180527/2.cutadapt
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/g.20180527
for i in $ddir/*.R1.fastq.gz; do key=$(basename $i .R1.fastq.gz)
mkdir -p $wdir/0.input/$key
ln -s $i $wdir/0.input/$key
ln -s ${i/R1/R2} $wdir/0.input/$key
done
nohup HiC-Pro -i $wdir/0.input -o $wdir/1.sequence -c $ConfigHP/config-hicpro-p32.txt > $wdir/nohup.20180515.log &

#连接validpairs 和stats到1.align
ln $wdir/1.sequence/hic_results/data/*/*validPairs ~/workspace/8.NT-HiC/1.align/a.all_validPairs/
ln $wdir/1.sequence/bowtie_results/bwt2/*/*bwt2pairs.pairstat ~/workspace/8.NT-HiC/1.align/b.all_stats/
#for this report

#2.合并seq到rep，为了方便，即使只有一个seq，也做一遍
for i in $ddir/*.R1.fastq.gz; do key=${i##*/};key=${key%_*}; 
mkdir -p $wdir/2.replicate_valid/$key; 
ln -s ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/2.replicate_valid/$key; 
done
(nohup HiC-Pro -i $wdir/2.replicate_valid -o $wdir/3.replicate -c $ConfigHP/config-hicpro.txt -s merge_persample > $wdir/nohup.rep.log
nohup HiC-Pro -i $wdir/3.replicate/hic_results/data -o $wdir/3.replicate -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.rep.log)&

#3.合并rep到sample，可以与上一步同步进行
for i in $ddir/*.R1.fastq.gz; do 
key=${i##*/};key=${key%_rep*}; 
mkdir -p $wdir/4.sample_valids/$key; 
ln -s ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/4.sample_valids/$key; 
done
(nohup HiC-Pro -i $wdir/4.sample_valids -o $wdir/5.sample -c $ConfigHP/config-hicpro.txt -s merge_persample > $wdir/nohup.sample.log
nohup HiC-Pro -i $wdir/5.sample/hic_results/data -o $wdir/5.sample -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.sample.log )&
#4.连接所有数据到2.mergesample和5.map
fs=(1.sequence 3.replicate 5.sample)
gs=(1.sequence 2.replicate 3.sample)
for i in {0..2}; do f=${fs[i]};g=${gs[i]}
ln -f $wdir/$f/hic_results/data/*/*allValidPairs ~/workspace/8.NT-HiC/2.merge_sample/$g/1.allValidPairs/
ln -f $wdir/$f/hic_results/data/*/*mergestat ~/workspace/8.NT-HiC/2.merge_sample/$g/2.stats/
ln -f $wdir/$f/hic_results/matrix/*/raw/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/1.raw
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/2.iced
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*biases ~/workspace/8.NT-HiC/5.maps/$g/3.biases
done
#5.统计日志到表格
Rmd HiC-Pro_stat.Rmd
#6.有两个新的sample: icm te，没有改变6h和morula
ln -s -f ~/workspace/8.NT-HiC/2.merge_sample/3.sample/1.allValidPairs/[!6m]* ~/workspace/8.NT-HiC/2.merge_sample/4.except/1.allValidPairs/
ln -s -f ~/workspace/8.NT-HiC/5.maps/3.sample/1.raw/[!6m]* ~/workspace/8.NT-HiC/5.maps/4.except/1.raw/
ln -s -f ~/workspace/8.NT-HiC/5.maps/3.sample/2.iced/[!6m]* ~/workspace/8.NT-HiC/5.maps/4.except/2.iced/
ln -s -f ~/workspace/8.NT-HiC/5.maps/3.sample/3.biases/[!6m]* ~/workspace/8.NT-HiC/5.maps/4.except/3.biases/

#####################################################################################################################################################################
############################2018-06-07   CC_rep3 CC_rep4 te_rep2 12h_rep2
#####################################################################################################################################################################
ddir=~/workspace/8.NT-HiC/0.Pre/n.date/4.20180607/2.cutadapt
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/h.20180607
for i in $ddir/*.R1.fastq.gz; do key=$(basename $i .R1.fastq.gz)
mkdir -p $wdir/0.input/$key
ln -s $i $wdir/0.input/$key
ln -s ${i/R1/R2} $wdir/0.input/$key
done
nohup HiC-Pro -i $wdir/0.input -o $wdir/1.sequence -c $ConfigHP/config-hicpro-p32.txt > $wdir/nohup.20180607.log &

#连接validpairs 和stats到1.align
ln $wdir/1.sequence/hic_results/data/*/*validPairs ~/workspace/8.NT-HiC/1.align/a.all_validPairs/
ln $wdir/1.sequence/bowtie_results/bwt2/*/*bwt2pairs.pairstat ~/workspace/8.NT-HiC/1.align/b.all_stats/
#for this report

#2.合并seq到rep，为了方便，即使只有一个seq，也做一遍
for i in $ddir/*.R1.fastq.gz; do key=${i##*/};key=${key%_*}; 
mkdir -p $wdir/2.replicate_valid/$key; 
ln -s ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/2.replicate_valid/$key; 
done
(nohup HiC-Pro -i $wdir/2.replicate_valid -o $wdir/3.replicate -c $ConfigHP/config-hicpro.txt -s merge_persample > $wdir/nohup.rep.log
nohup HiC-Pro -i $wdir/3.replicate/hic_results/data -o $wdir/3.replicate -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.rep.log)&

#3.合并rep到sample，可以与上一步同步进行
for i in $ddir/*.R1.fastq.gz; do 
key=${i##*/};key=${key%_rep*}; 
mkdir -p $wdir/4.sample_valids/$key; 
ln -s ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/4.sample_valids/$key; 
done
(nohup HiC-Pro -i $wdir/4.sample_valids -o $wdir/5.sample -c $ConfigHP/config-hicpro.txt -s merge_persample > $wdir/nohup.sample.log
nohup HiC-Pro -i $wdir/5.sample/hic_results/data -o $wdir/5.sample -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.sample.log )&

#4.连接所有数据到2.mergesample和5.map
fs=(1.sequence 3.replicate 5.sample)
gs=(1.sequence 2.replicate 3.sample)
for i in {0..2}; do f=${fs[i]};g=${gs[i]}
ln -f $wdir/$f/hic_results/data/*/*allValidPairs ~/workspace/8.NT-HiC/2.merge_sample/$g/1.allValidPairs/
ln -f $wdir/$f/hic_results/data/*/*mergestat ~/workspace/8.NT-HiC/2.merge_sample/$g/2.stats/
ln -f $wdir/$f/hic_results/matrix/*/raw/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/1.raw
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/2.iced
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*biases ~/workspace/8.NT-HiC/5.maps/$g/3.biases
done

#5.统计日志到表格
Rmd HiC-Pro_stat.Rmd
#6.排除某些重复；有两个新的sample: cc 12h，没有改变6h和morula
ln -s -f ~/workspace/8.NT-HiC/2.merge_sample/3.sample/1.allValidPairs/[!6m]* ~/workspace/8.NT-HiC/2.merge_sample/4.except/1.allValidPairs/
ln -s -f ~/workspace/8.NT-HiC/5.maps/3.sample/1.raw/[!6m]* ~/workspace/8.NT-HiC/5.maps/4.except/1.raw/
ln -s -f ~/workspace/8.NT-HiC/5.maps/3.sample/2.iced/[!6m]* ~/workspace/8.NT-HiC/5.maps/4.except/2.iced/
ln -s -f ~/workspace/8.NT-HiC/5.maps/3.sample/3.biases/[!6m]* ~/workspace/8.NT-HiC/5.maps/4.except/3.biases/

#####################################################################################################################################################################
#################2018-06-25 06-09   te_rep1_seq2 1h_rep3_seq2 icm_rep2_seq1 te_rep3_seq1
#####################################################################################################################################################################
ddir=~/workspace/8.NT-HiC/0.Pre/n.date/5.20180624/2.cutadapt
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/i.20180624
for i in $ddir/*.R1.fastq.gz; do key=$(basename $i .R1.fastq.gz)
mkdir -p $wdir/0.input/$key
ln -s $i $wdir/0.input/$key
ln -s ${i/R1/R2} $wdir/0.input/$key
done
nohup HiC-Pro -i $wdir/0.input -o $wdir/1.sequence -c $ConfigHP/config-hicpro-p32.txt > $wdir/nohup.20180624.log &

#连接validpairs 和stats到1.align
ln $wdir/1.sequence/hic_results/data/*/*validPairs ~/workspace/8.NT-HiC/1.align/a.all_validPairs/
ln $wdir/1.sequence/bowtie_results/bwt2/*/*bwt2pairs.pairstat ~/workspace/8.NT-HiC/1.align/b.all_stats/
#for this report

#2.合并seq到rep，为了方便，即使只有一个seq，也做一遍
for i in $ddir/*.R1.fastq.gz; do key=${i##*/};key=${key%_*}; 
mkdir -p $wdir/2.replicate_valid/$key; 
ln -s ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/2.replicate_valid/$key; 
done
(nohup HiC-Pro -i $wdir/2.replicate_valid -o $wdir/3.replicate -c $ConfigHP/config-hicpro.txt -s merge_persample > $wdir/nohup.rep.log
nohup HiC-Pro -i $wdir/3.replicate/hic_results/data -o $wdir/3.replicate -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.rep.log)&

#3.合并rep到sample，可以与上一步同步进行
for i in $ddir/*.R1.fastq.gz; do 
key=${i##*/};key=${key%_rep*}; 
mkdir -p $wdir/4.sample_valids/$key; 
ln -s -f ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/4.sample_valids/$key; 
done  #如果有两个同一样本的重复，有File exists错误
(nohup HiC-Pro -i $wdir/4.sample_valids -o $wdir/5.sample -c $ConfigHP/config-hicpro.txt -s merge_persample > $wdir/nohup.sample.log
nohup HiC-Pro -i $wdir/5.sample/hic_results/data -o $wdir/5.sample -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.sample.log )&

#4.连接所有数据到2.mergesample和5.map
fs=(1.sequence 3.replicate 5.sample)
gs=(1.sequence 2.replicate 3.sample)
for i in {0..2}; do f=${fs[i]};g=${gs[i]}
ln -f $wdir/$f/hic_results/data/*/*allValidPairs ~/workspace/8.NT-HiC/2.merge_sample/$g/1.allValidPairs/
ln -f $wdir/$f/hic_results/data/*/*mergestat ~/workspace/8.NT-HiC/2.merge_sample/$g/2.stats/
ln -f $wdir/$f/hic_results/matrix/*/raw/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/1.raw
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/2.iced
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*biases ~/workspace/8.NT-HiC/5.maps/$g/3.biases
done

#5.统计日志到表格
Rmd HiC-Pro_stat.Rmd
#6.排除某些重复；有两个新的sample: cc 12h，没有改变6h和morula
ln -s -f ~/workspace/8.NT-HiC/2.merge_sample/3.sample/1.allValidPairs/[!6m]* ~/workspace/8.NT-HiC/2.merge_sample/4.except/1.allValidPairs/
ln -s -f ~/workspace/8.NT-HiC/5.maps/3.sample/1.raw/[!6m]* ~/workspace/8.NT-HiC/5.maps/4.except/1.raw/
ln -s -f ~/workspace/8.NT-HiC/5.maps/3.sample/2.iced/[!6m]* ~/workspace/8.NT-HiC/5.maps/4.except/2.iced/
ln -s -f ~/workspace/8.NT-HiC/5.maps/3.sample/3.biases/[!6m]* ~/workspace/8.NT-HiC/5.maps/4.except/3.biases/

#####################################################################################################################################################################
##2018-07-04 12h_rep2加测  05h_rep4
#####################################################################################################################################################################
ddir=~/workspace/8.NT-HiC/0.Pre/n.date/6.20180704/2.cutadapt
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/j.20180704
for i in $ddir/*.R1.fastq.gz; do key=$(basename $i .R1.fastq.gz)
mkdir -p $wdir/0.input/$key
ln -s $i $wdir/0.input/$key
ln -s ${i/R1/R2} $wdir/0.input/$key
done
nohup HiC-Pro -i $wdir/0.input -o $wdir/1.sequence -c $ConfigHP/config-hicpro-p32.txt > $wdir/nohup.align.log &

#连接validpairs 和stats到1.align
ln $wdir/1.sequence/hic_results/data/*/*validPairs ~/workspace/8.NT-HiC/1.align/a.all_validPairs/
ln $wdir/1.sequence/bowtie_results/bwt2/*/*bwt2pairs.pairstat ~/workspace/8.NT-HiC/1.align/b.all_stats/
#for this report

#2.合并seq到rep，为了方便，即使只有一个seq，也做一遍
for i in $ddir/*.R1.fastq.gz; do key=${i##*/};key=${key%_*}; 
mkdir -p $wdir/2.replicate_valid/$key; 
ln -s ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/2.replicate_valid/$key; 
done
(nohup HiC-Pro -i $wdir/2.replicate_valid -o $wdir/3.replicate -c $ConfigHP/config-hicpro.txt -s merge_persample > $wdir/nohup.rep.log
nohup HiC-Pro -i $wdir/3.replicate/hic_results/data -o $wdir/3.replicate -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.rep.log)&

#3.合并rep到sample，可以与上一步同步进行
for i in $ddir/*.R1.fastq.gz; do 
key=${i##*/};key=${key%_rep*}; 
mkdir -p $wdir/4.sample_valids/$key; 
ln -s -f ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/4.sample_valids/$key; 
done  #如果有两个同一样本的重复，有File exists错误
(nohup HiC-Pro -i $wdir/4.sample_valids -o $wdir/5.sample -c $ConfigHP/config-hicpro.txt -s merge_persample > $wdir/nohup.sample.log
nohup HiC-Pro -i $wdir/5.sample/hic_results/data -o $wdir/5.sample -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.sample.log )&

#4.连接所有数据到2.mergesample和5.map
fs=(1.sequence 3.replicate 5.sample)
gs=(1.sequence 2.replicate 3.sample)
for i in {0..2}; do f=${fs[i]};g=${gs[i]}
ln -f $wdir/$f/hic_results/data/*/*allValidPairs ~/workspace/8.NT-HiC/2.merge_sample/$g/1.allValidPairs/
ln -f $wdir/$f/hic_results/data/*/*mergestat ~/workspace/8.NT-HiC/2.merge_sample/$g/2.stats/
ln -f $wdir/$f/hic_results/matrix/*/raw/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/1.raw
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/2.iced
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*biases ~/workspace/8.NT-HiC/5.maps/$g/3.biases
done

#5.统计日志到表格
Rmd HiC-Pro_stat.Rmd
#6.排除某些重复；05h_rep1 4cell_rep1 6h_rep2 morula_rep2
ln -s -f ~/workspace/8.NT-HiC/2.merge_sample/3.sample/1.allValidPairs/[!6m]* ~/workspace/8.NT-HiC/2.merge_sample/4.except/1.allValidPairs/
ln -s -f ~/workspace/8.NT-HiC/5.maps/3.sample/1.raw/[!6m]* ~/workspace/8.NT-HiC/5.maps/4.except/1.raw/
ln -s -f ~/workspace/8.NT-HiC/5.maps/3.sample/2.iced/[!6m]* ~/workspace/8.NT-HiC/5.maps/4.except/2.iced/
ln -s -f ~/workspace/8.NT-HiC/5.maps/3.sample/3.biases/[!6m]* ~/workspace/8.NT-HiC/5.maps/4.except/3.biases/

#重新确认了删除的重复05h_rep1 4cell_rep1 6h_rep2 morula_rep2
for i in 05h_rep1 4cell_rep1 6h_rep2 morula_rep2; do key=${i%%_*}
mkdir -p $wdir/6.except_valids/$key
ln -s -f ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/6.except_valids/$key
rm $wdir/6.except_valids/$key/$i*
done
(nohup HiC-Pro -i $wdir/6.except_valids -o $wdir/7.except -c $ConfigHP/config-hicpro.txt -s merge_persample > $wdir/nohup.except.log
nohup HiC-Pro -i $wdir/7.except/hic_results/data -o $wdir/7.except -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.except.log )&



#####################################################################################################################################################################
##2018-07-04 icm_rep1_seq2 icm_rep2_seq2
#####################################################################################################################################################################
ddir=~/workspace/8.NT-HiC/0.Pre/n.date/7.20180719/2.cutadapt
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/k.20180720
for i in $ddir/*.R1.fastq.gz; do key=$(basename $i .R1.fastq.gz)
mkdir -p $wdir/0.input/$key
ln -s $i $wdir/0.input/$key
ln -s ${i/R1/R2} $wdir/0.input/$key
done
nohup HiC-Pro -i $wdir/0.input -o $wdir/1.sequence -c $ConfigHP/config-hicpro-p32.txt > $wdir/nohup.align.log &

#连接validpairs 和stats到1.align
ln $wdir/1.sequence/hic_results/data/*/*validPairs ~/workspace/8.NT-HiC/1.align/a.all_validPairs/
ln $wdir/1.sequence/bowtie_results/bwt2/*/*bwt2pairs.pairstat ~/workspace/8.NT-HiC/1.align/b.all_stats/
#for this report

#2.合并seq到rep，为了方便，即使只有一个seq，也做一遍
for i in $ddir/*.R1.fastq.gz; do key=${i##*/};key=${key%_*}; 
mkdir -p $wdir/2.replicate_valid/$key; 
ln -s ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/2.replicate_valid/$key; 
done
(nohup HiC-Pro -i $wdir/2.replicate_valid -o $wdir/3.replicate -c $ConfigHP/config-hicpro.txt -s merge_persample > $wdir/nohup.rep.log
nohup HiC-Pro -i $wdir/3.replicate/hic_results/data -o $wdir/3.replicate -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.rep.log)&

#3.合并rep到sample，可以与上一步同步进行
for i in $ddir/*.R1.fastq.gz; do 
key=${i##*/};key=${key%_rep*}; 
mkdir -p $wdir/4.sample_valids/$key; 
ln -s -f ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/4.sample_valids/$key; 
done  #如果有两个同一样本的重复，有File exists错误
(nohup HiC-Pro -i $wdir/4.sample_valids -o $wdir/5.sample -c $ConfigHP/config-hicpro.txt -s merge_persample > $wdir/nohup.sample.log
nohup HiC-Pro -i $wdir/5.sample/hic_results/data -o $wdir/5.sample -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.sample.log )&

#4.连接所有数据到2.mergesample和5.map
fs=(1.sequence 3.replicate 5.sample)
gs=(1.sequence 2.replicate 3.sample)
for i in {0..2}; do f=${fs[i]};g=${gs[i]}
ln -f $wdir/$f/hic_results/data/*/*allValidPairs ~/workspace/8.NT-HiC/2.merge_sample/$g/1.allValidPairs/
ln -f $wdir/$f/hic_results/data/*/*mergestat ~/workspace/8.NT-HiC/2.merge_sample/$g/2.stats/
ln -f $wdir/$f/hic_results/matrix/*/raw/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/1.raw
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/2.iced
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*biases ~/workspace/8.NT-HiC/5.maps/$g/3.biases
done


################################################
#去掉一个05h_rep1 4
##############################################
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/l.20180801
key=05h
mkdir -p $wdir/6.except_valids/$key
ln -s -f ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/6.except_valids/$key
rm $wdir/6.except_valids/05h/05h_rep[14]*

(nohup HiC-Pro -i $wdir/6.except_valids -o $wdir/7.except -c $ConfigHP/config-hicpro.txt -s merge_persample > $wdir/nohup.except.log
nohup HiC-Pro -i $wdir/7.except/hic_results/data -o $wdir/7.except -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.except.log )&

ln -f $wdir/7.except/hic_results/data/*/*allValidPairs ~/workspace/8.NT-HiC/2.merge_sample/4.except/1.allValidPairs/
ln -f $wdir/7.except/hic_results/data/*/*mergestat ~/workspace/8.NT-HiC/2.merge_sample/4.except/2.stats/
ln -f $wdir/7.except/hic_results/matrix/*/raw/*/*matrix ~/workspace/8.NT-HiC/5.maps/4.except/1.raw
ln -f $wdir/7.except/hic_results/matrix/*/iced/*/*matrix ~/workspace/8.NT-HiC/5.maps/4.except/2.iced
ln -f $wdir/7.except/hic_results/matrix/*/iced/*/*biases ~/workspace/8.NT-HiC/5.maps/4.except/3.biases

#重新决定只使用6h_rep3 6h_rep5
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/e.20180505
f=7.except
g=4.except
ln -f $wdir/$f/hic_results/data/*/*allValidPairs ~/workspace/8.NT-HiC/2.merge_sample/$g/1.allValidPairs/
ln -f $wdir/$f/hic_results/data/*/*mergestat ~/workspace/8.NT-HiC/2.merge_sample/$g/2.stats/
ln -f $wdir/$f/hic_results/matrix/*/raw/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/1.raw
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/2.iced
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*biases ~/workspace/8.NT-HiC/5.maps/$g/3.biases


##去除1h_rep3 2h_rep3
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/m.20180813
key=2h #1h
mkdir -p $wdir/6.except_valids/$key
ln -s -f ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/6.except_valids/$key
rm $wdir/6.except_valids/*/*_rep3*

(nohup HiC-Pro -i $wdir/6.except_valids -o $wdir/7.except -c $ConfigHP/config-hicpro.txt -s merge_persample > $wdir/nohup.except.log
nohup HiC-Pro -i $wdir/7.except/hic_results/data -o $wdir/7.except -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.except.log )&

ln -f $wdir/7.except/hic_results/data/*/*allValidPairs ~/workspace/8.NT-HiC/2.merge_sample/4.except/1.allValidPairs/
ln -f $wdir/7.except/hic_results/data/*/*mergestat ~/workspace/8.NT-HiC/2.merge_sample/4.except/2.stats/
ln -f $wdir/7.except/hic_results/matrix/*/raw/*/*matrix ~/workspace/8.NT-HiC/5.maps/4.except/1.raw
ln -f $wdir/7.except/hic_results/matrix/*/iced/*/*matrix ~/workspace/8.NT-HiC/5.maps/4.except/2.iced
ln -f $wdir/7.except/hic_results/matrix/*/iced/*/*biases ~/workspace/8.NT-HiC/5.maps/4.except/3.biases


#####################################################################################################################################################################
##2018-08-24 12h_rep3
#####################################################################################################################################################################
ddir=~/workspace/8.NT-HiC/0.Pre/n.date/8.20180823/2.cutadapt
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/n.20180823
for i in $ddir/*.R1.fastq.gz; do key=$(basename $i .R1.fastq.gz)
mkdir -p $wdir/0.input/$key
ln -s $i $wdir/0.input/$key
ln -s ${i/R1/R2} $wdir/0.input/$key
done
nohup HiC-Pro -i $wdir/0.input -o $wdir/1.sequence -c $ConfigHP/config-hicpro-p32.txt > $wdir/nohup.align.log &

#连接validpairs 和stats到1.align
ln $wdir/1.sequence/hic_results/data/*/*validPairs ~/workspace/8.NT-HiC/1.align/a.all_validPairs/
ln $wdir/1.sequence/bowtie_results/bwt2/*/*bwt2pairs.pairstat ~/workspace/8.NT-HiC/1.align/b.all_stats/
#for this report

#2.合并seq到rep，为了方便，即使只有一个seq，也做一遍
for i in $ddir/*.R1.fastq.gz; do key=${i##*/};key=${key%_*}; 
mkdir -p $wdir/2.replicate_valid/$key; 
ln -s ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/2.replicate_valid/$key; 
done
(nohup HiC-Pro -i $wdir/2.replicate_valid -o $wdir/3.replicate -c $ConfigHP/config-hicpro.txt -s merge_persample > $wdir/nohup.rep.log
nohup HiC-Pro -i $wdir/3.replicate/hic_results/data -o $wdir/3.replicate -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.rep.log)&

#3.合并rep到sample，可以与上一步同步进行
for i in $ddir/*.R1.fastq.gz; do 
key=${i##*/};key=${key%_rep*}; 
mkdir -p $wdir/4.sample_valids/$key; 
ln -s -f ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/4.sample_valids/$key; 
done  #如果有两个同一样本的重复，有File exists错误
(nohup HiC-Pro -i $wdir/4.sample_valids -o $wdir/5.sample -c $ConfigHP/config-hicpro.txt -s merge_persample > $wdir/nohup.sample.log
nohup HiC-Pro -i $wdir/5.sample/hic_results/data -o $wdir/5.sample -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.sample.log )&

#4.连接所有数据到2.mergesample和5.map
fs=(1.sequence 3.replicate 5.sample 5.sample)
gs=(1.sequence 2.replicate 3.sample 4.except)
for i in {0..2}; do f=${fs[i]};g=${gs[i]}
ln -f $wdir/$f/hic_results/data/*/*allValidPairs ~/workspace/8.NT-HiC/2.merge_sample/$g/1.allValidPairs/
ln -f $wdir/$f/hic_results/data/*/*mergestat ~/workspace/8.NT-HiC/2.merge_sample/$g/2.stats/
ln -f $wdir/$f/hic_results/matrix/*/raw/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/1.raw
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/2.iced
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*biases ~/workspace/8.NT-HiC/5.maps/$g/3.biases
done


###############################删除的rep
dels=(05h_rep1 05h_rep4 1h_rep3 2h_rep3 3h_rep1 4cell_rep1 6h_rep1 6h_rep2 6h_rep3 6h_rep4 e2cell_rep1 e2cell_rep2 morula_rep2)
for i in ${dels[@]}; do 
mv ~/workspace/8.NT-HiC/2.merge_sample/2.replicate/1.allValidPairs/${i}_allValidPairs ~/workspace/8.NT-HiC/2.merge_sample/2.replicate/3.remove/
done

05h_rep2	05h_rep3	12h_rep1	12h_rep2	12h_rep3	1h_rep1	1h_rep2	2h_rep1	2h_rep2	4cell_rep2	4cell_rep3	6h_rep5	6h_rep6	8cell_rep1	8cell_rep2	8cell_rep3	cc_rep1	cc_rep2	cc_rep3	cc_rep4	e2cell_rep3	e2cell_rep4	icm_rep1	icm_rep2	kdm4d-e2cell_rep1	kdm4d-e2cell_rep2	l2cell_rep1	l2cell_rep2	morula_rep1	morula_rep3	morula_rep4	partheno-1hpa_rep1	partheno-1hpa_rep2	partheno-e2cell_rep1	partheno-e2cell_rep2	sertoli-e2cell_rep1	te_rep1	te_rep2	te_rep3	tsa-e2cell_rep1



#####################################################################################################################################################################
#2019-08-30 revise第一次补测数据
#6h_rep6_seq1/  e2cell_rep4_seq1/  kdm4d-e2cell_rep1_seq1/  kdm4d-e2cell_rep2_seq1/  sertoli-e2cell_rep1_seq1/
#####################################################################################################################################################################
ddir=~/workspace/8.NT-HiC/0.Pre/n.date/9.20190829/2.cutadapt
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/o.20190829
for i in $ddir/*.R1.fastq.gz; do key=$(basename $i .R1.fastq.gz)
mkdir -p $wdir/0.input/$key
ln -s $i $wdir/0.input/$key
ln -s ${i/R1/R2} $wdir/0.input/$key
done
nohup HiC-Pro -i $wdir/0.input -o $wdir/1.sequence -c $ConfigHP/config-hicpro-p64.txt > $wdir/nohup.align.log &

#连接validpairs 和stats到1.align
ln $wdir/1.sequence/hic_results/data/*/*validPairs ~/workspace/8.NT-HiC/1.align/a.all_validPairs/
ln $wdir/1.sequence/bowtie_results/bwt2/*/*bwt2pairs.pairstat ~/workspace/8.NT-HiC/1.align/b.all_stats/
#for this report

#2.合并seq到rep，为了方便，即使只有一个seq，也做一遍
for i in $ddir/*.R1.fastq.gz; do key=${i##*/};key=${key%_*}; 
mkdir -p $wdir/2.replicate_valid/$key; 
ln -s ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/2.replicate_valid/$key; 
done
(nohup HiC-Pro -i $wdir/2.replicate_valid -o $wdir/3.replicate -c $ConfigHP/config-hicpro.txt -s merge_persample > $wdir/nohup.rep.log
nohup HiC-Pro -i $wdir/3.replicate/hic_results/data -o $wdir/3.replicate -c $ConfigHP/config-hicpro.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.rep.log)&
###test
(nohup HiC-Pro -i $wdir/2.replicate_valid -o $wdir/6.test -c $ConfigHP/config-hicpro-p64.txt -s merge_persample > $wdir/nohup.rep.log
nohup HiC-Pro -i $wdir/6.test/hic_results/data -o $wdir/6.test -c $ConfigHP/config-hicpro-p64.txt  -s build_contact_maps -s ice_norm >> $wdir/nohup.rep2.log)&



#3.合并rep到sample，可以与上一步同步进行;
keys=(6h e2cell kdm4d-e2cell sertoli-e2cell)
for key in ${keys[@]}; do 
mkdir -p $wdir/4.sample_valids/$key; 
ln -s -f ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/4.sample_valids/$key; 
done
#去除了不要的replicates
for i in 6h_rep1_seq1 6h_rep2_seq1 6h_rep3_seq1 6h_rep4 e2cell_rep1_seq1 e2cell_rep2_seq1;do
rm $wdir/4.sample_valids/*/$i*validPairs;
done
rm $wdir/4.sample_valids/sertoli-e2cell #去除只有一个重复的sample
(nohup HiC-Pro -i $wdir/4.sample_valids -o $wdir/5.sample -c $ConfigHP/config-hicpro-p64.txt -s merge_persample -s build_contact_maps -s ice_norm> $wdir/nohup.sample.log)&

#4.连接所有数据到2.mergesample和5.map
#不要sample
fs=(1.sequence 3.replicate)
gs=(1.sequence 2.replicate)
for i in {0..2}; do f=${fs[i]};g=${gs[i]}
ln -f $wdir/$f/hic_results/data/*/*allValidPairs ~/workspace/8.NT-HiC/2.merge_sample/$g/1.allValidPairs/
ln -f $wdir/$f/hic_results/data/*/*mergestat ~/workspace/8.NT-HiC/2.merge_sample/$g/2.stats/
ln -f $wdir/$f/hic_results/matrix/*/raw/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/1.raw/
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/2.iced/
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*biases ~/workspace/8.NT-HiC/5.maps/$g/3.biases/
done
f=5.sample;g=4.except
ln -f $wdir/$f/hic_results/data/*/*allValidPairs ~/workspace/8.NT-HiC/2.merge_sample/$g/1.allValidPairs/
ln -f $wdir/$f/hic_results/data/*/*mergestat ~/workspace/8.NT-HiC/2.merge_sample/$g/2.stats/
ln -f $wdir/$f/hic_results/matrix/*/raw/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/1.raw/
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/2.iced/
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*biases ~/workspace/8.NT-HiC/5.maps/$g/3.biases/


###########BLASTn 
ddir=~/workspace/8.NT-HiC/0.Pre/n.date/9.20190829/2.cutadapt
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/o.20190829/a.blastn
step=400000
for i in $ddir/6h*fastq.gz;do o=$(basename $i .fastq.gz)
(zcat  $i| awk -v step=$step 'NR%(4*step)==1{print ">"$0} NR%(4*step)==2{print $0}' > $wdir/$o.pick.fa
nohup blastn -max_target_seqs 1 -num_threads 16 -query $wdir/$o.pick.fa -db /home/share/BlastDB/blastdb_new/nt_new/nt \
-outfmt "7 qacc sacc evalue pident gaps qcovus stitle ssciname" \
-out $wdir/$o.self -evalue 1e-30 -perc_identity 99.33 
echo -e 'count\tratio\tspecies' > $o.statistic.tab
awk 'BEGIN {FS="[\t]+" ; OFS="\t"} $0!~/^#/{organism[$8]+=1 ; total+=1} END {for (item in organism) print organism[item],organism[item]/total*100,item }' $o.self | sort -k2nr,2 >> $o.statistic.tab  )&
done


#####################################################################################################################################################################
#2019-08-30 revise第二次补测数据
#partheno-1hpa_rep1_seq1/  partheno-1hpa_rep2_seq1/  partheno-e2cell_rep1_seq1/  partheno-e2cell_rep2_seq1/  tsa-e2cell_rep1_seq1/
#####################################################################################################################################################################
ddir=~/workspace/8.NT-HiC/0.Pre/n.date/a.20190907/2.cutadapt
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/p.20190907
for i in $ddir/*.R1.fastq.gz; do key=$(basename $i .R1.fastq.gz)
mkdir -p $wdir/0.input/$key
ln -s $i $wdir/0.input/$key
ln -s ${i/R1/R2} $wdir/0.input/$key
done
nohup HiC-Pro -i $wdir/0.input -o $wdir/1.sequence -c $ConfigHP/config-hicpro-p64.txt > $wdir/nohup.align.log &

#连接validpairs 和stats到1.align
ln $wdir/1.sequence/hic_results/data/*/*validPairs ~/workspace/8.NT-HiC/1.align/a.all_validPairs/
ln $wdir/1.sequence/bowtie_results/bwt2/*/*bwt2pairs.pairstat ~/workspace/8.NT-HiC/1.align/b.all_stats/

#2.合并seq到rep，为了方便，即使只有一个seq，也做一遍
for i in $ddir/*.R1.fastq.gz; do key=${i##*/};key=${key%_*}; 
mkdir -p $wdir/2.replicate_valid/$key; 
ln -s ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/2.replicate_valid/$key; 
done
(nohup HiC-Pro -i $wdir/2.replicate_valid -o $wdir/3.replicate -c $ConfigHP/config-hicpro-p64.txt -s merge_persample -s build_contact_maps -s ice_norm > $wdir/nohup.rep.log)&

#3.合并rep到sample，可以与上一步同步进行
for key in  partheno-1hpa partheno-e2cell; do 
mkdir -p $wdir/4.sample_valids/$key; 
ln -s -f ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/4.sample_valids/$key; 
done 
#去除只有一个重复的sample
rm tsa-e2cell*
(nohup HiC-Pro -i $wdir/4.sample_valids -o $wdir/5.sample -c $ConfigHP/config-hicpro-p64.txt -s merge_persample -s build_contact_maps -s ice_norm >> $wdir/nohup.sample.log > $wdir/nohup.sample.log)&

#4.连接所有数据到2.mergesample和5.map
fs=(1.sequence 3.replicate 5.sample 5.sample)
gs=(1.sequence 2.replicate 3.sample 4.except)
for i in {0..3}; do f=${fs[i]};g=${gs[i]}
ln -f $wdir/$f/hic_results/data/*/*allValidPairs ~/workspace/8.NT-HiC/2.merge_sample/$g/1.allValidPairs/
ln -f $wdir/$f/hic_results/data/*/*mergestat ~/workspace/8.NT-HiC/2.merge_sample/$g/2.stats/
ln -f $wdir/$f/hic_results/matrix/*/raw/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/1.raw
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/2.iced
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*biases ~/workspace/8.NT-HiC/5.maps/$g/3.biases
done



#####################################################################################################################################################################
#2019-08-30 revise第三次补测数据
#sertoli-e2cell_rep2_seq1 tsa-e2cell_rep2_seq1
#####################################################################################################################################################################
ddir=~/workspace/8.NT-HiC/0.Pre/n.date/b.20190930/2.cutadapt
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/q.20190930
for i in $ddir/*.R1.fastq.gz; do key=$(basename $i .R1.fastq.gz)
mkdir -p $wdir/0.input/$key
ln -s $i $wdir/0.input/$key
ln -s ${i/R1/R2} $wdir/0.input/$key
done
nohup HiC-Pro -i $wdir/0.input -o $wdir/1.sequence -c $ConfigHP/config-hicpro-p64.txt > $wdir/nohup.align.log &

#连接validpairs 和stats到1.align
ln $wdir/1.sequence/hic_results/data/*/*validPairs ~/workspace/8.NT-HiC/1.align/a.all_validPairs/
ln $wdir/1.sequence/bowtie_results/bwt2/*/*bwt2pairs.pairstat ~/workspace/8.NT-HiC/1.align/b.all_stats/


#2.合并seq到rep，为了方便，即使只有一个seq，也做一遍
for i in $ddir/*.R1.fastq.gz; do key=${i##*/};key=${key%_*}; 
mkdir -p $wdir/2.replicate_valid/$key; 
ln -s ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/2.replicate_valid/$key; 
done
(nohup HiC-Pro -i $wdir/2.replicate_valid -o $wdir/3.replicate -c $ConfigHP/config-hicpro-p64.txt -s merge_persample -s build_contact_maps -s ice_norm > $wdir/nohup.rep.log)&

#3.合并rep到sample，可以与上一步同步进行
for key in sertoli-e2cell tsa-e2cell; do 
mkdir -p $wdir/4.sample_valids/$key; 
ln -s -f ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/4.sample_valids/$key; 
done 

(nohup HiC-Pro -i $wdir/4.sample_valids -o $wdir/5.sample -c $ConfigHP/config-hicpro-p64.txt -s merge_persample -s build_contact_maps -s ice_norm >> $wdir/nohup.sample.log > $wdir/nohup.sample.log)&

fs=(1.sequence 3.replicate 5.sample 5.sample)
gs=(1.sequence 2.replicate 3.sample 4.except)
for i in {0..3}; do f=${fs[i]};g=${gs[i]}
ln -f $wdir/$f/hic_results/data/*/*allValidPairs ~/workspace/8.NT-HiC/2.merge_sample/$g/1.allValidPairs/
ln -f $wdir/$f/hic_results/data/*/*mergestat ~/workspace/8.NT-HiC/2.merge_sample/$g/2.stats/
ln -f $wdir/$f/hic_results/matrix/*/raw/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/1.raw
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/2.iced
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*biases ~/workspace/8.NT-HiC/5.maps/$g/3.biases
done

#####################################################################################################################################################################
#2019-10-18 revise第三次补测数据 补测数据
#sertoli-e2cell_rep2_seq2 tsa-e2cell_rep2_seq2
#####################################################################################################################################################################
ddir=~/workspace/8.NT-HiC/0.Pre/n.date/c.20191018/2.cutadapt
wdir=~/workspace/8.NT-HiC/c.HiC-Pro/r.20191018
for i in $ddir/*.R1.fastq.gz; do key=$(basename $i .R1.fastq.gz)
mkdir -p $wdir/0.input/$key
ln -s $i $wdir/0.input/$key
ln -s ${i/R1/R2} $wdir/0.input/$key
done
nohup HiC-Pro -i $wdir/0.input -o $wdir/1.sequence -c $ConfigHP/config-hicpro-p64.txt > $wdir/nohup.align.log &

#连接validpairs 和stats到1.align
ln $wdir/1.sequence/hic_results/data/*/*validPairs ~/workspace/8.NT-HiC/1.align/a.all_validPairs/
ln $wdir/1.sequence/bowtie_results/bwt2/*/*bwt2pairs.pairstat ~/workspace/8.NT-HiC/1.align/b.all_stats/


#2.合并seq到rep，为了方便，即使只有一个seq，也做一遍
for i in $ddir/*.R1.fastq.gz; do key=${i##*/};key=${key%_*}; 
mkdir -p $wdir/2.replicate_valid/$key; 
ln -s ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/2.replicate_valid/$key; 
done
(nohup HiC-Pro -i $wdir/2.replicate_valid -o $wdir/3.replicate -c $ConfigHP/config-hicpro-p64.txt -s merge_persample -s build_contact_maps -s ice_norm > $wdir/nohup.rep.log)&

#3.合并rep到sample，可以与上一步同步进行
for key in sertoli-e2cell tsa-e2cell; do 
mkdir -p $wdir/4.sample_valids/$key; 
ln -s -f ~/workspace/8.NT-HiC/1.align/a.all_validPairs/$key*._mm10.bwt2pairs.validPairs $wdir/4.sample_valids/$key; 
done 

(nohup HiC-Pro -i $wdir/4.sample_valids -o $wdir/5.sample -c $ConfigHP/config-hicpro-p64.txt -s merge_persample -s build_contact_maps -s ice_norm >> $wdir/nohup.sample.log > $wdir/nohup.sample.log)&

fs=(1.sequence 3.replicate 5.sample 5.sample)
gs=(1.sequence 2.replicate 3.sample 4.except)
for i in {0..3}; do f=${fs[i]};g=${gs[i]}
ln -f $wdir/$f/hic_results/data/*/*allValidPairs ~/workspace/8.NT-HiC/2.merge_sample/$g/1.allValidPairs/
ln -f $wdir/$f/hic_results/data/*/*mergestat ~/workspace/8.NT-HiC/2.merge_sample/$g/2.stats/
ln -f $wdir/$f/hic_results/matrix/*/raw/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/1.raw
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*matrix ~/workspace/8.NT-HiC/5.maps/$g/2.iced
ln -f $wdir/$f/hic_results/matrix/*/iced/*/*biases ~/workspace/8.NT-HiC/5.maps/$g/3.biases
done



