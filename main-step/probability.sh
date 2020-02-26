################################################all sample probability
ddir=~/workspace/8.NT-HiC/d.HiC-Pro_analysis/1.merge_samples/3.matrix
wdir=~/workspace/8.NT-HiC/a.probability/1.all_samples_probability
#取100k分辨率的matrix，分出interchrom和intrachrom
res=100000
for i in $ddir/*${res}_iced.matrix; do o=${i##*/};#o=${o/_100000_iced.matrix}
awk -v out=$o -v wdir=$wdir 'NR==FNR{if($4>a[$1]){a[$1]=$4}} NR>FNR{for(i in a){if($1<=a[i]&&$2>a[i]){print $0 > wdir"/inter."out;next}};
print $0 > wdir"/intra."out}' $ConfigHP/${res}_mm10.bed $i &
done
#将所有的等距的互作数加和，距离以100k为单位记录
for i in $wdir/intra*matrix;do o=${i##*intra.}; o=${o/matrix/tab};
awk '{a[$2-$1]+=$3} END{for(i in a)print i,a[i]}' $i | sort -k1n,1 > $wdir/$o & done


 
ddir=~/workspace/8.NT-HiC/d.HiC-Pro_analysis/a.merge_samples_except/3.matrix
wdir=~/workspace/8.NT-HiC/a.probability/2.all_except_probability


##################################################chr1 cc probability
ddir=~/workspace/8.NT-HiC/d.HiC-Pro_analysis/a.merge_samples_except/3.matrix
wdir=~/workspace/8.NT-HiC/a.probability/4.all_except_20180323/chrs
res=40000
#chr=chr2
key=6h
for chr in chr2 chr14 chr19 chrX; do
awk -v out=$o -v wdir=$wdir -v chr=$chr 'NR==FNR&&$1==chr{a[$4]} NR>FNR{if(($1 in a) && ($2 in a)){print $0}}' $ConfigHP/${res}_mm10.bed $ddir/${key}_${res}_iced.matrix | \
awk '{a[$2-$1]+=$3} END{for(i in a)print i,a[i]}' | sort -k1n,1 > $wdir/${key}_$chr.tab &
done



###############################################2018年3月23日   
ddir=~/workspace/8.NT-HiC/d.HiC-Pro_analysis/a.merge_samples_except/3.matrix
wdir=~/workspace/8.NT-HiC/a.probability/4.all_except_20180323
#取100k分辨率的matrix，分出interchrom和intrachrom
res=40000
for i in $ddir/*${res}_iced.matrix; do o=${i##*/};#o=${o/_100000_iced.matrix}
awk -v out=$o -v wdir=$wdir 'NR==FNR{if($4>a[$1]){a[$1]=$4}} NR>FNR{for(i in a){if($1<=a[i]&&$2>a[i]){print $0 > wdir"/inter."out;next}};
print $0 > wdir"/intra."out}' $ConfigHP/${res}_mm10.bed $i &
done
wait 
#将所有的等距的互作数加和，距离以100k为单位记录
for i in $wdir/intra*matrix;do o=${i##*intra.}; o=${o/matrix/tab};
awk '{a[$2-$1]+=$3} END{for(i in a)print i,a[i]}' $i | sort -k1n,1 > $wdir/$o & done


###################################2018年3月29日 测试是否可以用Ps~s-1 normalize intrachrom的matrix
wdir=~/workspace/8.NT-HiC/a.probability/4.all_except_20180323/matrix
for i in ../intra*_100000_iced.matrix; do
awk '{$3=$3*($2-$1+1)/(4000000/(100000));print $0}' $i > ${i#*intra.} &
done

i=chr19
res=100000
ddir=~/workspace/8.NT-HiC/a.probability/4.all_except_20180323/matrix/1.intra_test
wdir=~/workspace/8.NT-HiC/a.probability/4.all_except_20180323/matrix/1.intra_test
python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix \
$ddir/05h_${res}_iced.matrix \
$ddir/1h_${res}_iced.matrix \
$ddir/2h_${res}_iced.matrix \
$ddir/6h_${res}_iced.matrix \
$ddir/e2cell_${res}_iced.matrix \
$ddir/l2cell_${res}_iced.matrix \
$ddir/4cell_${res}_iced.matrix \
$ddir/8cell_${res}_iced.matrix \
$ddir/morula_${res}_iced.matrix \
-tri 1 -bed $ConfigHP/${res}_mm10.bed -n ${keys[@]} -chr $i -o $wdir/$i -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1


#####################################2018年4月1日 normalize for depth
###############################################2018年3月23日   
ddir=~/workspace/8.NT-HiC/d.HiC-Pro_analysis/a.merge_samples_except/4.norm_depth
wdir=~/workspace/8.NT-HiC/a.probability/5.norm_for_depth
#取40k分辨率的matrix，分出interchrom和intrachrom
res=40000
for i in $ddir/*${res}_iced.matrix; do o=${i##*/}; #o=${o/_100000_iced.matrix}
awk -v out=$o -v wdir=$wdir 'NR==FNR{if($4>a[$1]){a[$1]=$4}} NR>FNR{for(i in a){if($1<=a[i]&&$2>a[i]){print $0 > wdir"/inter."out;next}};
print $0 > wdir"/intra."out}' $ConfigHP/${res}_mm10.bed $i &
done
wait 
#将所有的等距的互作数加和，距离以100k为单位记录
for i in $wdir/intra*_${res}_iced.matrix;do o=${i##*intra.}; o=${o/matrix/tab};
awk '{a[$2-$1]+=$3} END{for(i in a)print i,a[i]}' $i | sort -k1n,1 > $wdir/$o & done

########es500
ddir=~/workspace/8.NT-HiC/6.hic-pro_20171123/f.merge_es500_20180404/2.norm_to_depth


######################################2018年4月12日 + 6h 4cell morula 
bash ~/codes/hic_probability.sh -d ~/workspace/8.NT-HiC/d.HiC-Pro_analysis/1.merge_samples_except_20180411/4.norm_to_depth -w ~/workspace/8.NT-HiC/a.probability/6.all_except_20180412 


#########################################2018年5月13日 + 6h_rep5 3h 12h e2cell l2cell
bash ~/codes/hic_probability.sh -d ~/workspace/8.NT-HiC/5.maps/4.except/6.norm_to_mean -w ~/workspace/8.NT-HiC/a.probability/7.all_except_20180513 
##########+ 2h_rep3 同上



#########################################2018年5月13日 + 6h_rep5 3h 12h e2cell l2cell
bash ~/codes/hic_probability.sh -d ~/workspace/8.NT-HiC/5.maps/4.except/6.norm_to_mean -w ~/workspace/8.NT-HiC/a.probability/7.all_except_20180513 

#########################################2018年5月13日 + 6h_rep5 3h 12h e2cell l2cell
bash ~/codes/hic_probability.sh -d ~/workspace/8.NT-HiC/5.maps/4.except/5.norm_to_mean -w ~/workspace/8.NT-HiC/a.probability/8.all_except_20180531 -r 40000
#try 100k
bash ~/codes/hic_probability.sh -d ~/workspace/8.NT-HiC/5.maps/4.except/5.norm_to_mean -w ~/workspace/8.NT-HiC/a.probability/8.all_except_20180531 -r 100000


###########################2018年6月8日 测试10kb的raw data
bash ~/codes/hic_probability.sh -d ~/workspace/8.NT-HiC/5.maps/4.except/1.raw -w ~/workspace/8.NT-HiC/a.probability/9.test_raw_20180608 -r 10000
#raw data 是可以的
bash ~/codes/hic_probability.sh -d ~/workspace/8.NT-HiC/5.maps/4.except/1.raw -w ~/workspace/8.NT-HiC/a.probability/a.except_raw_20180608 -r 10000

######################2018年6月11日 allreplicate
bash ~/codes/hic_probability.sh -d ~/workspace/8.NT-HiC/5.maps/2.replicate/1.raw -w ~/workspace/8.NT-HiC/a.probability/b.allRep_raw_20180608 -r 10000


############### ~ 1M ~ 10M ~ 互作数比例
ddir=~/workspace/8.NT-HiC/5.maps/2.replicate/1.raw
wdir=~/workspace/8.NT-HiC/a.probability/c.cc_1h_chrs_20180621
mkdir -p $wdir
for i in $ddir/cc*_10000.matrix $ddir/05h*_10000.matrix $ddir/1h*_10000.matrix; do o=$wdir/$(basename $i _10000.matrix).tab
python ~/codes/probability_chrs.py -i $i -o $o &
done

####分染色体
ddir=~/workspace/8.NT-HiC/5.maps/2.replicate/1.raw
wdir=~/workspace/8.NT-HiC/a.probability/c.cc_1h_chrs_20180621
bash ~/codes/hic_probability2.sh -d $ddir -w $wdir -r 10000
#结果不好


######################2018年6月29日 allreplicate
#2018年8月3日
bash ~/codes/hic_probability.sh -d ~/workspace/8.NT-HiC/5.maps/2.replicate/1.raw -w ~/workspace/8.NT-HiC/a.probability/d.allRep_raw_20180629 -r 10000


#2018年8月3日
bash ~/codes/hic_probability.sh -d ~/workspace/8.NT-HiC/5.maps/4.except/1.raw -w ~/workspace/8.NT-HiC/a.probability/a.except_raw_20180608 -r 10000
rm 3h_10000.*

#2018年9月5日
bash ~/codes/hic_probability.sh -d ~/workspace/8.NT-HiC/5.maps/2.replicate/1.raw -w ~/workspace/8.NT-HiC/a.probability/g.allRep_raw_20180905 -r 10000


#####################2018年8月26日 + 12h_rep3
bash ~/codes/hic_probability.sh -d ~/workspace/8.NT-HiC/5.maps/4.except/1.raw -w ~/workspace/8.NT-HiC/a.probability/f.except_raw_20180826 -r 10000


##################2018年11月28日 all rep for XW
bash ~/codes/hic_probability.sh -d ~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/7.replicate/3.raw -w ~/workspace/8.NT-HiC/a.probability/h.xw_allRep_20181128 -r 10000

bash ~/codes/hic_probability.sh -d ~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/3.raw -w ~/workspace/8.NT-HiC/a.probability/i.xw_20181128 -r 10000


bash ~/codes/hic_probability.sh -d ~/workspace/8.NT-HiC/5.maps/2.replicate/2.iced -w ~/workspace/8.NT-HiC/a.probability/j.allRep_iced_20181128 -r 10000

bash ~/codes/hic_probability.sh -d ~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/7.replicate/4.iced -w ~/workspace/8.NT-HiC/a.probability/k.xw_allRep_iced_20181128 -r 10000

#################################pairsqc test
wdir=~/workspace/8.NT-HiC/a.probability/l.pairsqc
ddir=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/1.allValidPairs

pairix -f -s1 -b3 -e4 -u8 -T /n/data1/hms/dbmi/park/sl325/4dn/pairs/GM12878_insitu_DpnII_merged_nodups.tab.bsorted.txt.gz
awk '{print $1,$2,$3,$5,$6,$4,$7}' $ddir/early_2cell_allValidPairs | sort -k2,2 -k4,4 -k3n,3 > early_2cell_allValidPairs
bgzip early_2cell_allValidPairs
pairix early_2cell_allValidPairs.gz
python ~/1.HiC-software/l.pairsqc/pairsqc.py -p $ddir/early_2cell_allValidPairs -c $ConfigHP/chrom_mm10.sizes -O $wdir/NF -s e2cell -M 8.2 &


python pairsqc.py -p test_samples/merged_nodup.tab.chrblock_sorted.txt.gz -c test_samples/GRCh37.chrom.sizes.mainonly.female -t M
Rscript plot.r
open report/pairsqc_report.html
zip report.zip report



##############################for allelic in zygote
#-r C57BL_6N female G1
#-a PWK_PhJ male G2
wdir=~/workspace/8.NT-HiC/a.probability/m.allelic_zygote_20181202

ln -s ~/workspace/8.NT-HiC/5.maps/4.except/1.raw/12h_10000.matrix $wdir/0.data/
ln -s ~/workspace/8.NT-HiC/5.maps/4.except/1.raw/6h_10000.matrix $wdir/0.data/
ln -s ~/workspace/8.NT-HiC/3.public_data/e.Allelic_xw/3.raw/PN3_zygote_10000_G1.matrix $wdir/0.data/PN3_mat_10000.matrix
ln -s ~/workspace/8.NT-HiC/3.public_data/e.Allelic_xw/3.raw/PN3_zygote_10000_G2.matrix $wdir/0.data/PN3_pat_10000.matrix
ln -s ~/workspace/8.NT-HiC/3.public_data/e.Allelic_xw/3.raw/PN5_zygote_10000_G1.matrix $wdir/0.data/PN5_mat_10000.matrix
ln -s ~/workspace/8.NT-HiC/3.public_data/e.Allelic_xw/3.raw/PN5_zygote_10000_G2.matrix $wdir/0.data/PN5_pat_10000.matrix
ln -s ~/workspace/8.NT-HiC/3.public_data/j.Flyamer/3.raw/1.mat_10000.matrix $wdir/0.data/Flyamer_mat_10000.matrix
ln -s ~/workspace/8.NT-HiC/3.public_data/j.Flyamer/3.raw/2.pat_10000.matrix $wdir/0.data/Flyamer_pat_10000.matrix

bash ~/codes/hic_probability.sh -d $wdir/0.data -w $wdir -r 10000



########################revise 1hpa
ddir=~/workspace/8.NT-HiC/5.maps/4.except/1.raw
wdir=~/workspace/8.NT-HiC/a.probability/n.partheno_20191006
bash ~/codes/hic_probability.sh -d $ddir -w $wdir -r 10000



ddir1=~/workspace/8.NT-HiC/5.maps/4.except/1.raw
ddir2=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/3.raw
wdir=~/workspace/8.NT-HiC/a.probability/p.e2cell-except_20191011
mkd $wdir
ln -s $ddir1/*e2cell*10000.matrix $wdir
ln -s $ddir2/early_2cell*10000.matrix $wdir
bash ~/codes/hic_probability.sh -d $wdir -w $wdir -r 10000


