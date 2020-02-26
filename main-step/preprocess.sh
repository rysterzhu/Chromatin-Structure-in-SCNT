 keys=("cc" "05h" "1h" "2h" "6h" "e2cell" "l2cell" "4cell" "8cell" "morula")
####################cutadapt
py35
wdir=~/workspace/8.NT-HiC/0.Pre
for i in $wdir/[c8012E]*/0.data/*R1.fastq.gz; do o=${i/0.data/2.cutadapt}
nohup cutadapt -j 16 --trim-n -q 25,25 -m 75 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log}
done &

wdir=~/workspace/8.NT-HiC/0.Pre
for i in $wdir/[6el4m]*/0.data/*R1.fastq.gz; do o=${i/0.data/2.cutadapt}
nohup cutadapt -j 16 --trim-n -q 25,25 -m 75 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log}
done &


#####################################4cell rep2 二次建库 第一次的加测
ddir=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_Hi_C_20171103/20180501
wdir=~/workspace/8.NT-HiC/0.Pre
ln -s $ddir/CM20180208-3_L3_A008.R1.clean.fastq.gz $wdir/4cell/0.data/4cell_rep2_seq3.R1.fastq.gz
ln -s $ddir/CM20180208-3_L3_A008.R2.clean.fastq.gz $wdir/4cell/0.data/4cell_rep2_seq3.R2.fastq.gz

nohup fastqc -o 1.qc/ 0.data/4cell_rep2_seq3.R* -t 2 &

i=$wdir/4cell/0.data/4cell_rep2_seq3.R1.fastq.gz;o=${i/0.data/2.cutadapt}
nohup cutadapt --trim-n -q 25,25 -m 75 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &


##################################20180505 6h_rep5_seq2 3h_rep1_seq1 12h_rep1_seq1 e2cell_rep3_seq1 l2cell_rep2_seq1
ddir=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_Hi_C_20171103/20180505/1.rawdata
wdir=~/workspace/8.NT-HiC/0.Pre/n.date/1.20180505
mkdir -p $wdir/1.qc $wdir/2.cutadapt $wdir/3.cutadapt_qc

ln -s $ddir/CM-0325-3_TKD180302111/CM-0325-3_TKD180302111_1.fq.gz $wdir/6h_rep5_seq2.R1.fastq.gz
ln -s $ddir/CM0420-1_TKD180401857/CM0420-1_TKD180401857_1.fq.gz $wdir/3h_rep1_seq1.R1.fastq.gz
ln -s $ddir/CM0420-2_TKD180401858/CM0420-2_TKD180401858_1.fq.gz $wdir/12h_rep1_seq1.R1.fastq.gz
ln -s $ddir/CM0420-3_TKD180401859/CM0420-3_TKD180401859_1.fq.gz $wdir/e2cell_rep3_seq1.R1.fastq.gz
ln -s $ddir/CM0420-4_TKD180401860/CM0420-4_TKD180401860_1.fq.gz $wdir/l2cell_rep2_seq1.R1.fastq.gz

ln -s $ddir/CM-0325-3_TKD180302111/CM-0325-3_TKD180302111_2.fq.gz $wdir/6h_rep5_seq2.R2.fastq.gz
ln -s $ddir/CM0420-1_TKD180401857/CM0420-1_TKD180401857_2.fq.gz $wdir/3h_rep1_seq1.R2.fastq.gz
ln -s $ddir/CM0420-2_TKD180401858/CM0420-2_TKD180401858_2.fq.gz $wdir/12h_rep1_seq1.R2.fastq.gz
ln -s $ddir/CM0420-3_TKD180401859/CM0420-3_TKD180401859_2.fq.gz $wdir/e2cell_rep3_seq1.R2.fastq.gz
ln -s $ddir/CM0420-4_TKD180401860/CM0420-4_TKD180401860_2.fq.gz $wdir/l2cell_rep2_seq1.R2.fastq.gz

nohup fastqc -o $wdir/1.qc $wdir/*gz -t 10 &

py35
for i in $wdir/*R1.fastq.gz; do o=$wdir/2.cutadapt/${i##*/};
nohup cutadapt -j 8 --trim-n -q 25,25 -m 75 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
nohup fastqc -o $wdir/3.cutadapt_qc $wdir/2.cutadapt/*gz -t 10 &


##################################2018年5月15日 1h_rep3_seq1
ddir=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_Hi_C_20171103/20180514
wdir=~/workspace/8.NT-HiC/0.Pre/n.date/2.20180514
mkdir -p $wdir/1.qc $wdir/2.cutadapt $wdir/3.cutadapt_qc
ln -s $ddir/CM-0425-1_TKD180402509_1.fq.gz $wdir/2h_rep3_seq1.R1.fastq.gz
ln -s $ddir/CM-0425-1_TKD180402509_2.fq.gz $wdir/2h_rep3_seq1.R2.fastq.gz


#################################2018年5月27日    20180524 + 20180527 两次的数据一起做了
ddir1=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_Hi_C_20171103/20180524/1.rawdata
ddir2=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_Hi_C_20171103/20180527/1.rawdata
wdir=~/workspace/8.NT-HiC/0.Pre/n.date/3.20180527
mkdir -p $wdir/1.qc $wdir/2.cutadapt $wdir/3.cutadapt_qc

ln -s $ddir1/CM-0425-1_TKD180402509_1.fq.gz $wdir/2h_rep3_seq2.R1.fastq.gz
ln -s $ddir1/CM-0425-1_TKD180402509_2.fq.gz $wdir/2h_rep3_seq2.R2.fastq.gz
ln -s $ddir1/CM0420-2_TKD180401858_1.fq.gz $wdir/12h_rep1_seq2.R1.fastq.gz
ln -s $ddir1/CM0420-2_TKD180401858_2.fq.gz $wdir/12h_rep1_seq2.R2.fastq.gz

ln -s $ddir2/CM0512-1_TKD180501003_1.fq.gz $wdir/te_rep1_seq1.R1.fastq.gz
ln -s $ddir2/CM0512-2_TKD180501004_1.fq.gz $wdir/icm_rep1_seq1.R1.fastq.gz
ln -s $ddir2/CM0512-3_TKD180501005_1.fq.gz $wdir/1h_rep3_seq1.R1.fastq.gz
ln -s $ddir2/CM0512-1_TKD180501003_2.fq.gz $wdir/te_rep1_seq1.R2.fastq.gz
ln -s $ddir2/CM0512-2_TKD180501004_2.fq.gz $wdir/icm_rep1_seq1.R2.fastq.gz
ln -s $ddir2/CM0512-3_TKD180501005_2.fq.gz $wdir/1h_rep3_seq1.R2.fastq.gz

nohup fastqc -o $wdir/1.qc $wdir/*gz -t 20 &

py35
for i in $wdir/*R1.fastq.gz; do o=$wdir/2.cutadapt/${i##*/};
nohup cutadapt -j 8 --trim-n -q 25,25 -m 75 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
nohup fastqc -o $wdir/3.cutadapt_qc $wdir/2.cutadapt/*gz -t 10 &

#################################2018年6月6日    CC_rep3 CC_rep4 te_rep2 12h_rep2
ddir=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_Hi_C_20171103/20180606/1.rawdata
wdir=~/workspace/8.NT-HiC/0.Pre/n.date/4.20180607
mkdir -p $wdir/1.qc $wdir/2.cutadapt $wdir/3.cutadapt_qc $wdir/0.raw

ln -s $ddir/CM0524-1_TKD180502220_1.fq.gz $wdir/0.raw/cc_rep3_seq1.R1.fastq.gz
ln -s $ddir/CM0524-1_TKD180502220_2.fq.gz $wdir/0.raw/cc_rep3_seq1.R2.fastq.gz
ln -s $ddir/CM0524-2_TKD180502221_1.fq.gz $wdir/0.raw/cc_rep4_seq1.R1.fastq.gz
ln -s $ddir/CM0524-2_TKD180502221_2.fq.gz $wdir/0.raw/cc_rep4_seq1.R2.fastq.gz
ln -s $ddir/CM0524-3_TKD180502222_1.fq.gz $wdir/0.raw/te_rep2_seq1.R1.fastq.gz
ln -s $ddir/CM0524-3_TKD180502222_2.fq.gz $wdir/0.raw/te_rep2_seq1.R2.fastq.gz
ln -s $ddir/CM0524-4_TKD180502223_1.fq.gz $wdir/0.raw/12h_rep2_seq1.R1.fastq.gz
ln -s $ddir/CM0524-4_TKD180502223_2.fq.gz $wdir/0.raw/12h_rep2_seq1.R2.fastq.gz

nohup fastqc -o $wdir/1.qc $wdir/0.raw/*gz -t 20 &

py35
for i in $wdir/0.raw/*R1.fastq.gz; do o=$wdir/2.cutadapt/${i##*/};
nohup cutadapt -j 8 --trim-n -q 25,25 -m 75 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
nohup fastqc -o $wdir/3.cutadapt_qc $wdir/2.cutadapt/*gz -t 10 &


##########################2018年6月9日 NT-TE-rep1 NT-TE-rep1 加测
###########################2018年6月24日 NT-ICM-rep2 NT-TE-rep3
ddir1=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_Hi_C_20171103/20180609/1.rawdata
ddir2=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_Hi_C_20171103/20180624/1.rawdata
wdir=~/workspace/8.NT-HiC/0.Pre/n.date/5.20180624
mkdir -p $wdir/1.qc $wdir/2.cutadapt $wdir/3.cutadapt_qc $wdir/0.raw

ln -s $ddir1/CM0512-1_TKD180501003_1.fq.gz $wdir/0.raw/te_rep1_seq2.R1.fastq.gz
ln -s $ddir1/CM0512-1_TKD180501003_2.fq.gz $wdir/0.raw/te_rep1_seq2.R2.fastq.gz
ln -s $ddir1/CM0512-3_TKD180501005_1.fq.gz $wdir/0.raw/1h_rep3_seq2.R1.fastq.gz
ln -s $ddir1/CM0512-3_TKD180501005_2.fq.gz $wdir/0.raw/1h_rep3_seq2.R2.fastq.gz

ln -s $ddir2/CM0607-1_TKD180600761_1.fq.gz $wdir/0.raw/icm_rep2_seq1.R1.fastq.gz
ln -s $ddir2/CM0607-1_TKD180600761_2.fq.gz $wdir/0.raw/icm_rep2_seq1.R2.fastq.gz
ln -s $ddir2/CM0607-2_TKD180600762_1.fq.gz $wdir/0.raw/te_rep3_seq1.R1.fastq.gz
ln -s $ddir2/CM0607-2_TKD180600762_2.fq.gz $wdir/0.raw/te_rep3_seq1.R2.fastq.gz

nohup fastqc -o $wdir/1.qc $wdir/0.raw/*gz -t 20 &

py35
for i in $wdir/0.raw/*R1.fastq.gz; do o=$wdir/2.cutadapt/${i##*/};
nohup cutadapt -j 8 --trim-n -q 25,25 -m 75 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
nohup fastqc -o $wdir/3.cutadapt_qc $wdir/2.cutadapt/*gz -t 10 &


#####################2018年7月4日 12h_rep2加测  05h_rep3
ddir=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_Hi_C_20171103/20180630/1.rawdata
wdir=~/workspace/8.NT-HiC/0.Pre/n.date/6.20180704
mkdir -p $wdir/1.qc $wdir/2.cutadapt $wdir/3.cutadapt_qc $wdir/0.raw

ln -s $ddir/CM0524-4_TKD180502223_1.fq.gz $wdir/0.raw/12h_rep2_seq2.R1.fastq.gz
ln -s $ddir/CM0524-4_TKD180502223_2.fq.gz $wdir/0.raw/12h_rep2_seq2.R2.fastq.gz
ln -s $ddir/CM0614-1_TKD180601378_1.fq.gz $wdir/0.raw/05h_rep4_seq1.R1.fastq.gz
ln -s $ddir/CM0614-1_TKD180601378_2.fq.gz $wdir/0.raw/05h_rep4_seq1.R2.fastq.gz

nohup fastqc -o $wdir/1.qc $wdir/0.raw/*gz -t 20 &

py35
for i in $wdir/0.raw/*R1.fastq.gz; do o=$wdir/2.cutadapt/${i##*/};
nohup cutadapt -j 8 --trim-n -q 25,25 -m 75 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
nohup fastqc -o $wdir/3.cutadapt_qc $wdir/2.cutadapt/*gz -t 10 &


######################2018年7月19日 icm_rep1_seq2 icm_rep2_seq2
ddir=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_Hi_C_20171103/20180718/1.rawdata
wdir=~/workspace/8.NT-HiC/0.Pre/n.date/7.20180719
mkdir -p $wdir/1.qc $wdir/2.cutadapt $wdir/3.cutadapt_qc $wdir/0.raw

ln -s $ddir/CM0512-2_TKD180501004_HNCM3CCXY_L2_1.fq.gz $wdir/0.raw/icm_rep1_seq2.R1.fastq.gz
ln -s $ddir/CM0512-2_TKD180501004_HNCM3CCXY_L2_2.fq.gz $wdir/0.raw/icm_rep1_seq2.R2.fastq.gz
ln -s $ddir/CM0607-1_TKD180600761_HNCM3CCXY_L2_1.fq.gz $wdir/0.raw/icm_rep2_seq2.R1.fastq.gz
ln -s $ddir/CM0607-1_TKD180600761_HNCM3CCXY_L2_2.fq.gz $wdir/0.raw/icm_rep2_seq2.R2.fastq.gz

nohup fastqc -o $wdir/1.qc $wdir/0.raw/*gz -t 20 &

py35
for i in $wdir/0.raw/*R1.fastq.gz; do o=$wdir/2.cutadapt/${i##*/};
nohup cutadapt -j 8 --trim-n -q 25,25 -m 75 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait; de
nohup fastqc -o $wdir/3.cutadapt_qc $wdir/2.cutadapt/*gz -t 10 &


######################2018年8月24日 12h_rep3
ddir=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_Hi_C_20171103/20180823
wdir=~/workspace/8.NT-HiC/0.Pre/n.date/8.20180823
mkdir -p $wdir/1.qc $wdir/2.cutadapt $wdir/3.cutadapt_qc $wdir/0.raw

ln -s $ddir/CM-0731-8_TKD180800127_HT2GHCCXY_L2_1.fq.gz $wdir/0.raw/12h_rep3_seq1.R1.fastq.gz
ln -s $ddir/CM-0731-8_TKD180800127_HT2GHCCXY_L2_2.fq.gz $wdir/0.raw/12h_rep3_seq1.R2.fastq.gz

nohup fastqc -o $wdir/1.qc $wdir/0.raw/*gz -t 20 &

py35
for i in $wdir/0.raw/*R1.fastq.gz; do o=$wdir/2.cutadapt/${i##*/};
nohup cutadapt -j 16 --trim-n -q 25,25 -m 75 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait; de
nohup fastqc -o $wdir/3.cutadapt_qc $wdir/2.cutadapt/*gz -t 10 &



######################2018年8月24日 
ddir=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_Hi_C_20171103/20190829
wdir=~/workspace/8.NT-HiC/0.Pre/n.date/9.20190829
mkdir -p $wdir/1.qc $wdir/2.cutadapt $wdir/3.cutadapt_qc $wdir/0.raw

ln -s $ddir/CM0819-1_FKDL190756888-1a_1.fq.gz $wdir/0.raw/kdm4d-e2cell_rep1_seq1.R1.fastq.gz
ln -s $ddir/CM0819-2_FKDL190756889-1a_1.fq.gz $wdir/0.raw/kdm4d-e2cell_rep2_seq1.R1.fastq.gz
ln -s $ddir/CM0819-3_FKDL190756890-1a_1.fq.gz $wdir/0.raw/6h_rep6_seq1.R1.fastq.gz
ln -s $ddir/CM0819-4_FKDL190756891-1a_1.fq.gz $wdir/0.raw/e2cell_rep4_seq1.R1.fastq.gz
ln -s $ddir/CM0819-5_FKDL190756892-1a_1.fq.gz $wdir/0.raw/sertoli-e2cell_rep1_seq1.R1.fastq.gz

ln -s $ddir/CM0819-1_FKDL190756888-1a_2.fq.gz $wdir/0.raw/kdm4d-e2cell_rep1_seq1.R2.fastq.gz
ln -s $ddir/CM0819-2_FKDL190756889-1a_2.fq.gz $wdir/0.raw/kdm4d-e2cell_rep2_seq1.R2.fastq.gz
ln -s $ddir/CM0819-3_FKDL190756890-1a_2.fq.gz $wdir/0.raw/6h_rep6_seq1.R2.fastq.gz
ln -s $ddir/CM0819-4_FKDL190756891-1a_2.fq.gz $wdir/0.raw/e2cell_rep4_seq1.R2.fastq.gz
ln -s $ddir/CM0819-5_FKDL190756892-1a_2.fq.gz $wdir/0.raw/sertoli-e2cell_rep1_seq1.R2.fastq.gz

nohup fastqc -o $wdir/1.qc $wdir/0.raw/*gz -t 20 &

py35
for i in $wdir/0.raw/*R1.fastq.gz; do o=$wdir/2.cutadapt/${i##*/};
nohup cutadapt -j 32 --trim-n -q 25,25 -m 50 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait-m cutadapt; de
nohup fastqc -o $wdir/3.cutadapt_qc $wdir/2.cutadapt/*gz -t 30 &
wait-m qc



#########################################2019年9月7日
linuxnd login -u P101SC18072044-01-F035 -p musjwctf 
linuxnd list oss://CP2018091214825/C101SC18072043/KY_kehu_JK/P101SC18072044-01/P101SC18072044-01-F035
linuxnd cp -d oss://CP2018091214825/C101SC18072043/KY_kehu_JK/P101SC18072044-01/P101SC18072044-01-F035/CM0829-2_FKDL190759435-1a_H327YCCX2_L4 /home/share/DATA/NT-HiC_qszhu/

linuxnd login -u P101SC18072044-01-F034 -p t1d7893x 
linuxnd list oss://CP2018091214825/C101SC18072043/KY_kehu_JK/P101SC18072044-01/P101SC18072044-01-F034/
linuxnd cp -d oss://CP2018091214825/C101SC18072043/KY_kehu_JK/P101SC18072044-01/P101SC18072044-01-F034/1.rawdata/CM0829_1_FKDL190759434-1a_H2YTNCCX2_L6 /home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_Hi_C_20171103/20190907

ddir=/home/share/DATA/clean_reads/GaoSR_mouse_nuclear_transfer_Hi_C_20171103/20190907
wdir=~/workspace/8.NT-HiC/0.Pre/n.date/a.20190907
mkdir -p $wdir/1.qc $wdir/2.cutadapt $wdir/3.cutadapt_qc $wdir/0.raw

ln -s $ddir/CM0823-1*_1.fq.gz $wdir/0.raw/partheno-1hpa_rep1_seq1.R1.fastq.gz
ln -s $ddir/CM0823-2*_1.fq.gz $wdir/0.raw/partheno-e2cell_rep1_seq1.R1.fastq.gz
ln -s $ddir/CM0823-3*_1.fq.gz $wdir/0.raw/partheno-e2cell_rep2_seq1.R1.fastq.gz
ln -s $ddir/CM0829_1*_1.fq.gz $wdir/0.raw/tsa-e2cell_rep1_seq1.R1.fastq.gz
ln -s $ddir/CM0829-2*_1.fq.gz $wdir/0.raw/partheno-1hpa_rep2_seq1.R1.fastq.gz

ln -s $ddir/CM0823-1*_2.fq.gz $wdir/0.raw/partheno-1hpa_rep1_seq1.R2.fastq.gz
ln -s $ddir/CM0823-2*_2.fq.gz $wdir/0.raw/partheno-e2cell_rep1_seq1.R2.fastq.gz
ln -s $ddir/CM0823-3*_2.fq.gz $wdir/0.raw/partheno-e2cell_rep2_seq1.R2.fastq.gz
ln -s $ddir/CM0829_1*_2.fq.gz $wdir/0.raw/tsa-e2cell_rep1_seq1.R2.fastq.gz
ln -s $ddir/CM0829-2*_2.fq.gz $wdir/0.raw/partheno-1hpa_rep2_seq1.R2.fastq.gz

nohup fastqc -o $wdir/1.qc $wdir/0.raw/*gz -t 20 &

conda activate /home/qszhu/Anaconda2/envs/py35
for i in $wdir/0.raw/*R1.fastq.gz; do o=$wdir/2.cutadapt/${i##*/};
nohup cutadapt -j 32 --trim-n -q 25,25 -m 50 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait-m cutadapt; de
nohup fastqc -o $wdir/3.cutadapt_qc $wdir/2.cutadapt/*gz -t 30 &
wait-m qc


#############
linuxnd login -u P101SC18072044-01-F038 -p 8g7znndd 
linuxnd list oss://CP2018091214825/C101SC18072043/KY_kehu_JK/P101SC18072044-01/P101SC18072044-01-F038/ 
linuxnd cp -d oss://CP2018091214825/C101SC18072043/KY_kehu_JK/P101SC18072044-01/P101SC18072044-01-F038/ /home/qszhu/data/GaoSR_mouse_nuclear_transfer_Hi_C_20171103/20190929

linuxnd login -u P101SC18072044-01-F039 -p 4dyaag4x 
linuxnd list oss://CP2018091214825/C101SC18072043/KY_kehu_JK/P101SC18072044-01/P101SC18072044-01-F039/ 
linuxnd cp -d oss://CP2018091214825/C101SC18072043/KY_kehu_JK/P101SC18072044-01/P101SC18072044-01-F039/ /home/qszhu/data/GaoSR_mouse_nuclear_transfer_Hi_C_20171103/20190929



ddir=/home/share/DATA/clean_reads/NT-HiC_qszhu/20190929
ddir=/mnt/pangu/home/share/DATA/clean_reads/NT-HiC_qszhu/20190929
wdir=~/workspace/8.NT-HiC/0.Pre/n.date/b.20190930
mkdir -p $wdir/1.qc $wdir/2.cutadapt $wdir/3.cutadapt_qc $wdir/0.raw

ln -s $ddir/CM0917-1_FKDL190763224-1a_H3TVTCCX2_L8_1.fq.gz $wdir/0.raw/sertoli-e2cell_rep2_seq1.R1.fastq.gz
ln -s $ddir/CM0917-2_FKDL190763225-1a_H3TVTCCX2_L7_1.fq.gz $wdir/0.raw/tsa-e2cell_rep2_seq1.R1.fastq.gz

ln -s $ddir/CM0917-1_FKDL190763224-1a_H3TVTCCX2_L8_2.fq.gz $wdir/0.raw/sertoli-e2cell_rep2_seq1.R2.fastq.gz
ln -s $ddir/CM0917-2_FKDL190763225-1a_H3TVTCCX2_L7_2.fq.gz $wdir/0.raw/tsa-e2cell_rep2_seq1.R2.fastq.gz

nohup fastqc -o $wdir/1.qc $wdir/0.raw/*gz -t 20 &

conda activate /home/qszhu/Anaconda2/envs/py35
for i in $wdir/0.raw/*R1.fastq.gz; do o=$wdir/2.cutadapt/${i##*/};
nohup cutadapt -j 32 --trim-n -q 25,25 -m 50 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait-m cutadapt;
nohup fastqc -o $wdir/3.cutadapt_qc $wdir/2.cutadapt/*gz -t 30 &
wait-m qc

######################################补测tsa sertoli rep2 seq2
ddir=~/data/GaoSR_mouse_nuclear_transfer_Hi_C_20171103/20191018
wdir=~/workspace/8.NT-HiC/0.Pre/n.date/c.20191018
mkdir -p $wdir/1.qc $wdir/2.cutadapt $wdir/3.cutadapt_qc $wdir/0.raw

ln -s $ddir/CM0917_1*_1.fq.gz $wdir/0.raw/sertoli-e2cell_rep2_seq2.R1.fastq.gz
ln -s $ddir/CM0917_2*_1.fq.gz $wdir/0.raw/tsa-e2cell_rep2_seq2.R1.fastq.gz

ln -s $ddir/CM0917_1*_2.fq.gz $wdir/0.raw/sertoli-e2cell_rep2_seq2.R2.fastq.gz
ln -s $ddir/CM0917_2*_2.fq.gz $wdir/0.raw/tsa-e2cell_rep2_seq2.R2.fastq.gz

nohup fastqc -o $wdir/1.qc $wdir/0.raw/*gz -t 20 &

conda activate /home/qszhu/Anaconda2/envs/py35
for i in $wdir/0.raw/*R1.fastq.gz; do o=$wdir/2.cutadapt/${i##*/};
nohup cutadapt -j 32 --trim-n -q 25,25 -m 50 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $o -p ${o/R1/R2} $i ${i/R1/R2} > ${o/R1.fastq.gz/log} &
done
wait-m cutadapt;
nohup fastqc -o $wdir/3.cutadapt_qc $wdir/2.cutadapt/*gz -t 30 &
wait-m qc




