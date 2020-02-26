odir=~/workspace/8.NT-HiC/5.maps
ddir=~/workspace/8.NT-HiC/2.merge_sample/1.sequence/1.allValidPairs

for i in $ddir/*; do k=${i##*/};k=${k%_*}
bash ~/codes/hic_pro_buildMaps.sh -i $i -w $odir/temp/$k -k $k -o $odir/1.sequence &
done


odir=~/workspace/8.NT-HiC/5.maps
ddir=~/workspace/8.NT-HiC/2.merge_sample/2.replicate/1.allValidPairs

for i in $ddir/*; do k=${i##*/};k=${k%_*}
bash ~/codes/hic_pro_buildMaps.sh -i $i -w $odir/temp/$k -k $k -o $odir/2.replicate &
done

odir=~/workspace/8.NT-HiC/5.maps
ddir=~/workspace/8.NT-HiC/2.merge_sample/3.sample/1.allValidPairs

for i in $ddir/*; do k=${i##*/};k=${k%_*}
bash ~/codes/hic_pro_buildMaps.sh -i $i -w $odir/temp/$k -k $k -o $odir/3.sample &
done
cp -r $odir/3.sample/* $odir/4.except
rm $odir/4.except/*/6h*
rm $odir/4.except/*/morula*

for i in ~/workspace/8.NT-HiC/2.merge_sample/4.except/1.allValidPairs/*; do k=${i##*/};k=${k%_*}
bash ~/codes/hic_pro_buildMaps.sh -i $i -w $odir/temp/$k -k $k -o $odir/4.except &
done
####################
##########################################这样不好，太麻烦了


################################################
###########normalize
################################################
ddir=~/workspace/8.NT-HiC/5.maps/4.except/2.iced
wdir=~/workspace/8.NT-HiC/5.maps/4.except

mkdir -p $wdir
for i in $ddir/*matrix; do 
awk -v OFS="\t" -v key=${i##*/} '$1!=$2{depth+=$3} END{print key,depth}' $i >> $wdir/depth.tab &
echo "calc depth: "$i
done

res=40000
f=1e8   #1e7 太低了，看不清
wdir=~/workspace/8.NT-HiC/5.maps/4.except/5.norm_to_1e8
mkdir -p $wdir
for i in $ddir/*$res*matrix; do 
awk -v OFS="\t" -v key=${i##*/} -v f=$f 'NR==FNR{if(key==$1)depth=$2} NR!=FNR{print $1,$2,$3*f/depth}' $wdir/../depth.tab $i > $wdir/${i##*/} &
echo "normalize: "$i
done

ddir=~/workspace/8.NT-HiC/5.maps/4.except/2.iced
wdir=~/workspace/8.NT-HiC/5.maps/4.except/6.norm_to_mean
res=500000
mkdir -p $wdir
f=`awk -v res=_${res}_ '$1~res{a+=$2;b+=1} END{print a/b}' $wdir/../depth.tab`  #mean depth of res
for i in $ddir/*_${res}_*matrix; do 
awk -v OFS="\t" -v key=${i##*/} -v f=$f 'NR==FNR{if(key==$1)depth=$2} NR!=FNR{print $1,$2,$3*f/depth}' $wdir/../depth.tab $i > $wdir/${i##*/} &
echo "normalize: "$i
done

##depth may be error
for i in $ddir/*matrix; do 
awk -v OFS="\t" -v key=${i##*/} '{depth+=$3} END{print key,depth}' $i >> $wdir/depth.tab &
echo "calc depth: "$i
done

bash ~/codes/hic_norm.sh -w ~/workspace/8.NT-HiC/5.maps/4.except -r 40000 -e 1e8 &
#没有什么问题

###################normalize 2018年5月31日
#1.计算每个的depth，先删除原来的，用raw计算，不同的resolution应该一样
ddir=~/workspace/8.NT-HiC/5.maps/4.except/1.raw
wdir=~/workspace/8.NT-HiC/5.maps/4.except
rm $wdir/depth.tab
for i in $ddir/*matrix; do 
awk -v OFS="\t" -v key=${i##*/} '{depth+=$3} END{print key,depth}' $i >> $wdir/depth.tab &
echo "calc depth: "$i
done
wait
sort -k1,1 $wdir/depth.tab > $wdir/depth.sorted.tab

#1.还是只用一个resolution计算depth就好了,resolution大计算快
ddir=~/workspace/8.NT-HiC/5.maps/4.except/1.raw
wdir=~/workspace/8.NT-HiC/5.maps/4.except
rm $wdir/depth.raw.tab
for i in $ddir/*_2000000.matrix; do k=$(basename $i _2000000.matrix)
awk -v OFS="\t" -v key=$k '{depth+=$3} END{print key,depth}' $i >> $wdir/depth.raw.tab &
done;wait
sort -k1,1 $wdir/depth.raw.tab | awk '{total+=$2;n+=1;print} END{print "total",total;printf("%s\t%d\n","mean",total/n)}' > $wdir/depth.sorted.tab
#平均 59343525 (59.3M pairs；最大e2cell 107M pairs； 可尝试norm到100M pairs
05h	18725549
12h	51366048
1h	45639641
2h	76086333
3h	17834870
4cell	98172467
6h	60195645
8cell	50251849
cc	51247669
e2cell	107577939
icm	30297331
l2cell	89933675
morula	87310712
te	46169629
total	830809357
mean	59343525

#2.norm to mean 
ddir=~/workspace/8.NT-HiC/5.maps/4.except/2.iced
odir=$wdir/4.norm_to_100M
mkdir -p $wdir
f=100000000
for res in 10000 25000 100000 250000 500000 1000000 2000000; do 
for i in $ddir/*_${res}_iced.matrix; do key=${i##*/};key=${key%%_*}
depth=`awk -v key=$key '{if(key==$1)print $2}' $wdir/depth.sorted.tab`
awk -v f=$f -v depth=$depth '{print $1,$2,$3*f/depth}' $i > $odir/${i##*/} &
echo "normalize: "$i
done;wait
done

#3.norm to mean
ddir=~/workspace/8.NT-HiC/5.maps/4.except/2.iced
odir=$wdir/5.norm_to_mean
mkdir -p $wdir
f=59343525 
for res in 10000 25000 40000 100000 250000 500000 1000000 2000000; do 
for i in $ddir/*_${res}_iced.matrix; do key=${i##*/};key=${key%%_*}
depth=`awk -v key=$key '{if(key==$1)print $2}' $wdir/depth.sorted.tab`
awk -v f=$f -v depth=$depth '{print $1,$2,$3*f/depth}' $i > $odir/${i##*/} &
echo "normalize: "$i
done;wait
done

#########################################2018年6月9日  终结版
###运行记录： 2018年6月29日； 
###2018年7月13日 换不同的rep选择；排除05h_rep1 4cell_rep1 6h_rep2 morula_rep2
###2018年7月27日 icm seq2
###2018年8月2日 去除 05h_rep4
###2018年8月13日 去除 1h_rep3 2h_rep3
###2018年8月25日  12h_rep3
################################################################################
#1.只用一个resolution计算depth就好了,resolution大计算快, 用raw计算depth更合理，并且与iced的一致
ddir=~/workspace/8.NT-HiC/5.maps/4.except/1.raw
wdir=~/workspace/8.NT-HiC/5.maps/4.except
mv $wdir/depth.sorted.tab $wdir/OLD_depth/depth.20180907.tab
rm $wdir/depth.raw.tab
for i in $wdir/1.raw/*_2000000.matrix; do k=$(basename $i _2000000.matrix)
awk -v OFS="\t" -v key=$k '{depth+=$3} END{print key,depth}' $i >> $wdir/depth.raw.tab &
done;wait
sort -k1,1 $wdir/depth.raw.tab | awk '{total+=$2;n+=1;print} END{print "total",total;printf("%s\t%d\n","mean",total/n)}' > $wdir/depth.sorted.tab
##平均 89247975 (59.3M pairs；最大e2cell 107M pairs； 可尝试norm到100M pairs
05h	13811561	13811561
12h	123499780	159254490
1h	66368273	11319110
2h	76086333	8923725
3h	17834870	17834870
4cell	100194900	100194900
6h	60195645	60195645
8cell	50251849	50251849
cc	137002285	137002285
e2cell	107577939	107577939
icm	123118472	123118472
l2cell	89933675	89933675
morula	87310712	87310712
te	165825537	165825537
total	1219011831	1132554770
mean	87072273	80896769
#2.norm to 100M
f=100000000; odir=4.norm_to_100M  #2018年7月27日 使用这个更好
#f=87803753; odir=5.norm_to_mean #2018年7月27日 
for res in 40000 100000 500000 10000 25000 250000 1000000 2000000; do 
for i in $wdir/2.iced/*_${res}_iced.matrix; do key=${i##*/};key=${key%%_*}
depth=`awk -v key=$key '{if(key==$1)print $2}' $wdir/depth.sorted.tab`
awk -v f=$f -v depth=$depth '{print $1,$2,$3*f/depth}' $i > $wdir/$odir/${i##*/} &
echo "normalize: "$i
done;wait
done

#3.norm raw data to mean
f=81747285; odir=6.norm_raw_to_mean
#mkdir -p $wdir/$odir
for res in 40000 100000 500000 10000 25000 250000 1000000 2000000; do 
for i in $wdir/1.raw/*_${res}.matrix; do key=${i##*/};key=${key%%_*}
depth=`awk -v key=$key '{if(key==$1)print $2}' $wdir/depth.sorted.tab`
awk -v f=$f -v depth=$depth '{print $1,$2,$3*f/depth}' $i > $wdir/$odir/${i##*/} &
echo "normalize: "$i
done;wait
done

######################################################
###分rep
ddir=~/workspace/8.NT-HiC/5.maps/2.replicate/1.raw
wdir=~/workspace/8.NT-HiC/5.maps/2.replicate

f=27877021; odir=5.norm_to_mean
for res in 40000 100000 500000; do 
for i in $wdir/2.iced/*_${res}_iced.matrix; do key=$(basename $i _${res}_iced.matrix)
depth=`awk -v key=$key '{if(key==$1)print $2}' $wdir/depth.sorted.tab`
awk -v f=$f -v depth=$depth '{print $1,$2,$3*f/depth}' $i > $wdir/$odir/${i##*/} &
echo "normalize: "$i
done;wait
done

f=27877021; odir=6.raw_to_mean
mkdir -p $wdir/$odir
for res in 100000 500000; do 
for i in $ddir/*_${res}.matrix; do key=$(basename $i _${res}.matrix)
depth=`awk -v key=$key '{if(key==$1)print $2}' $wdir/depth.sorted.tab`
awk -v f=$f -v depth=$depth '{print $1,$2,$3*f/depth}' $i > $wdir/$odir/${i##*/} &
echo "normalize: "$i
done;wait
done

###############2018年8月27日
#1.只用一个resolution计算depth就好了,resolution大计算快, 用raw计算depth更合理，并且与iced的一致
ddir=~/workspace/8.NT-HiC/5.maps/2.replicate/1.raw
wdir=~/workspace/8.NT-HiC/5.maps/2.replicate
mv $wdir/depth.sorted.tab $wdir/OLD_depth/depth.20180825.tab
rm $wdir/depth.raw.tab
for i in $wdir/1.raw/*_2000000.matrix; do k=$(basename $i _2000000.matrix)
awk -v OFS="\t" -v key=$k '{depth+=$3} END{print key,depth}' $i >> $wdir/depth.raw.tab &
done;wait
sort -k1,1 $wdir/depth.raw.tab | awk '{total+=$2;n+=1;print} END{print "total",total;printf("%s\t%d\n","mean",total/n)}' > $wdir/depth.sorted.tab
#2.norm to 100M
f=100000000; odir=4.norm_to_100M  #2018年7月27日 使用这个更好
mkdir -p $wdir/$odir
for res in 40000 100000 500000 10000 25000 250000 1000000 2000000; do 
for i in $wdir/2.iced/*_${res}_iced.matrix; do key=$(basename $i _${res}_iced.matrix)
depth=`awk -v key=$key '{if(key==$1)print $2}' $wdir/depth.sorted.tab`
awk -v f=$f -v depth=$depth '{print $1,$2,$3*f/depth}' $i > $wdir/$odir/${i##*/} &
echo "normalize: "$i
done;wait
done

###########################################
###########random
#########################################
wdir=~/workspace/8.NT-HiC/5.maps/4.except/5.random
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
mkdir -p $wdir
res=40000

awk 'NR==1{chr=$1;s=$4;e=$4} NR>1&&chr==$1{e=$4} NR>1&&chr!=$1{print chr,s,e;chr=$1;s=$4;e=$4} END{print chr,s,e}' $ConfigHP/${res}_mm10.bed > $wdir/incl_$res.bed

awk '{print $1,1,$2}' $ConfigHP/chrom_mm10.sizes | bedtools subtract -a - -b $wdir/incl_$res.bed > $wdir/excl_$res.bed

awk 'NR==FNR{s[$1]=$2;e[$1]=$3} NR>FNR{for(i in s){if(s[i]<=$1&&$2<=e[i]){print i,$1,$2,$3;next}}}' $wdir/incl_$res.bed $ddir/icm_${res}_iced.matrix > $wdir/icm_$res.temp  
bedtools shuffle -i $wdir/icm_$res.temp -excl $wdir/excl_$res.bed -chrom -seed 100 -g $ConfigHP/chrom_mm10.sizes -f 1E-9 > $wdir/icm_$res.random

sort -k2n,2 -k3n,3 icm_40000.random > icm_40000.temp2
cut -f 2-4 icm_40000.temp2 >icm_40000.matrix
##以上方法有问题，bedtools shuffle的incl参数不正确

awk '{srand($1);a=int(rand()*(68148-$2+$1));print a,a+$2-$1,$3}' ~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M/icm_40000_iced.matrix >temp
sort -k1n,1 -k2n,2 -S 5% temp > ~/workspace/8.NT-HiC/5.maps/4.except/5.random/2.DI/random_40000_iced.matrix
#这种方法对DI不行

awk 'NR==FNR{a[$1]=$2;b[$1]=$3} NR>FNR{for(i in a){if($1>=a[i]&&$2<=b[i]){
srand($1);c=int(rand()*(b[i]-a[i]-$2+$1)) + a[i];print c,c+$2-$1,$3}}}' $wdir/incl_40000.bed ~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M/icm_40000_iced.matrix > temp2 &
sort -k1n,1 -k2n,2 -S 5% temp2 > ~/workspace/8.NT-HiC/5.maps/4.except/5.random/random_40000_iced.matrix

#Insulation score
mkdir -p $wdir/1.IS/
for j in $(seq 1 19) X; do chr=chr$j
(python ~/codes/convert_3col_to_matrix_for_cworld.py -i random_40000_iced.matrix -I $ConfigHP/40000_mm10.bed -c $chr -o random.$chr.mat
perl -I ~/1.HiC-software/1.cworld/lib/ ~/1.HiC-software/1.cworld/perl/matrix2insulation.pl --is 1000000 --ids 200000 --nt 0.25 --bmoe 0 -i random.$chr.mat
)&
done
cat random.chr*boundaries.bed | awk -v OFS="\t" '$1~/^chr/{print}' | sort -k1V,1 -k2n,2 > random.boundary
cat random.chr*bedGraph | awk -v OFS="\t" '$1~/^chr/{if($4=="NA")$4=0;print}' | sort -k1V,1 -k2n,2 > random.insulation.bedGraph
rm *chr*
#Directional index
python ~/software/hic-pro/HiC-Pro_2.9.0/bin/utils/sparseToDense.py -o random -b $ConfigHP/40000_mm10.bed random_40000_iced.matrix --perchr -d
for i in $(seq 1 19) X; do 
perl ~/1.HiC-software/8.DI/1DI_from_matrix.pl chr${i}_random 40000 2000000 $ConfigHP/chrom_mm10.sizes > chr${i}.di
done
cat chr*.di | awk -v OFS="\t" '$1~/[0-9]+/{$0=$0;print} $1~"X"{$1=20;print} $1~"Y"{$1=21;print} ' | sort -k1n,1 -k2n,2 > all.di
awk -v OFS="\t" '$1==20{$1="X"} $1==21||$1==""||$2==""||$3==""||$4==""{next} {$1="chr"$1;print}' all.di | sort -k1V,1 -k2n,2 > random.bedGraph 

#####################################
#重新ICE
####################################
ddir=~/workspace/8.NT-HiC/5.maps/4.except/1.raw
wdir=~/workspace/8.NT-HiC/5.maps/4.except/6.ICED
for i in ${NT[@]}; do k=${i%%:*}
python /home/qszhu/software/hic-pro/HiC-Pro_2.9.0/scripts/ice --results_filename $wdir/${k}_40000_iced.matrix --filter_low_counts_perc 0.05 --filter_high_counts_perc 0 --max_iter 100 --eps 0.0001 --remove-all-zeros-loci --output-bias 1 --verbose 1 $ddir/${k}_40000.matrix >> $wdir/logs/$k.log &
done

wdir=~/workspace/8.NT-HiC/5.maps/4.except
f=100000000; odir=7.ICED_norm  #2018年7月27日 使用这个更好
for res in 40000; do 
for i in $wdir/6.ICED/*_${res}_iced.matrix; do key=${i##*/};key=${key%%_*}
depth=`awk -v key=$key '{if(key==$1)print $2}' $wdir/depth.sorted.tab`
awk -v f=$f -v depth=$depth '{print $1,$2,$3*f/depth}' $i > $wdir/$odir/${i##*/} &
echo "normalize: "$i
done;wait
done


###xw
ddir=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/3.raw
wdir=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/9.ICED
for i in ${NF[@]}; do k=${i%%:*}
python /home/qszhu/software/hic-pro/HiC-Pro_2.9.0/scripts/ice --results_filename $wdir/${k}_40000_iced.matrix --filter_low_counts_perc 0.05 --filter_high_counts_perc 0 --max_iter 100 --eps 0.0001 --remove-all-zeros-loci --output-bias 1 --verbose 1 $ddir/${k}_40000.matrix >> $wdir/logs/$k.log &
done

wdir=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw
f=100000000; odir=10.ICED_norm  #2018年7月27日 使用这个更好
for res in 40000; do 
for i in $wdir/9.ICED/*_${res}_iced.matrix; do key=${i##*/};key=${key%%_40000*}
depth=`awk -v key=$key '{if(key==$1)print $2}' $wdir/depth.sorted.tab`
awk -v f=$f -v depth=$depth '{print $1,$2,$3*f/depth}' $i > $wdir/$odir/${i##*/} &
echo "normalize: "$i
done;wait
done


#####################2019年8月19日 重新buildMap
ddir=~/workspace/8.NT-HiC/2.merge_sample/4.except/1.allValidPairs
wdir=~/workspace/8.NT-HiC/5.maps/4.except/a.rebuild
#mkd $wdir/logs
#for i in $ddir/*_allValidPairs; do k=$(basename $i _allValidPairs);mkd $wdir/data/$k;ln -s $i $wdir/data/$k;done
rm 6h* e2cell* #这两个sample补了重复
nohup HiC-Pro -i $wdir/data -o $wdir -c $ConfigHP/config-hicpro-p64.txt -s build_contact_maps -s ice_norm > $wdir/logs/build.log &

ln -f $wdir/hic_results/matrix/*/raw/*/*matrix ~/workspace/8.NT-HiC/5.maps/4.except/1.raw
ln -f $wdir/hic_results/matrix/*/iced/*/*matrix ~/workspace/8.NT-HiC/5.maps/4.except/2.iced
ln -f $wdir/hic_results/matrix/*/iced/*/*biases ~/workspace/8.NT-HiC/5.maps/4.except/3.biases


ddir=~/workspace/8.NT-HiC/5.maps/4.except/1.raw
wdir=~/workspace/8.NT-HiC/5.maps/4.except
mv $wdir/depth.sorted.tab $wdir/OLD_depth/depth.$(date +%Y-%m-%d).tab
rm $wdir/depth.raw.tab
for i in $wdir/1.raw/*_2000000.matrix; do k=$(basename $i _2000000.matrix)
awk -v OFS="\t" -v key=$k '{depth+=$3} END{print key,depth}' $i >> $wdir/depth.raw.tab &
done;wait
sort -k1,1 $wdir/depth.raw.tab | awk '{total+=$2;n+=1;print} END{print "total",total;printf("%s\t%d\n","mean",total/n)}' > $wdir/depth.sorted.tab

wdir=~/workspace/8.NT-HiC/5.maps/4.except
f=100000000; odir=4.norm_to_100M
mkd $wdir/$odir
for res in 500000; do 
for i in $wdir/2.iced/*_${res}_iced.matrix; do key=$(basename $i _${res}_iced.matrix)
depth=`awk -v key=$key '{if(key==$1)print $2}' $wdir/depth.sorted.tab`
awk -v f=$f -v depth=$depth '{print $1,$2,$3*f/depth}' $i > $wdir/$odir/${i##*/} &
echo "normalize: "$i
done;wait
done



###以下弃用
#######################Pangu
#####################################
#重新ICE Replicates 2019年9月12日
####################################
ddir=~/workspace/8.NT-HiC/5.maps/2.replicate/1.raw
wdir=~/workspace/8.NT-HiC/5.maps/2.replicate
mv $wdir/depth.sorted.tab $wdir/OLD_depth/depth.20190912.tab
rm $wdir/depth.raw.tab
for i in $wdir/1.raw/*_2000000.matrix; do k=$(basename $i _2000000.matrix)
awk -v OFS="\t" -v key=$k '{depth+=$3} END{print key,depth}' $i >> $wdir/depth.raw.tab &
done;wait
sort -k1,1 $wdir/depth.raw.tab | awk '{total+=$2;n+=1;print} END{print "total",total;printf("%s\t%d\n","mean",total/n)}' > $wdir/depth.sorted.tab


ddir=~/workspace/8.NT-HiC/5.maps/2.replicate/1.raw
wdir=~/workspace/8.NT-HiC/5.maps/2.replicate/6.ICED
mkd $wdir/logs
for res in 40000; do  #一次不能跑太多
for i in $ddir/[a-z]*_${res}.matrix; do k=$(basename $i _${res}.matrix)
python /usr/local/software/HiC-Pro_2.11.1/scripts/ice --results_filename $wdir/${k}_${res}_iced.matrix --filter_low_counts_perc 0.02 --filter_high_counts_perc 0 --max_iter 100 --eps 0.1 --remove-all-zeros-loci --output-bias 1 --verbose 1 $ddir/${k}_${res}.matrix >> $wdir/logs/${k}_${res}.log &
done;done

wdir=~/workspace/8.NT-HiC/5.maps/2.replicate
f=100000000; odir=7.ICED_norm  
mkd $wdir/$odir
for res in 40000 100000; do 
for i in $wdir/6.ICED/*_${res}_iced.matrix; do key=$(basename $i _${res}_iced.matrix)
depth=`awk -v key=$key '{if(key==$1)print $2}' $wdir/depth.sorted.tab`
awk -v f=$f -v depth=$depth '{print $1,$2,$3*f/depth}' $i > $wdir/$odir/${i##*/} &
echo "normalize: "$i
done;
done

###########################2019年9月13日 e2cell-rep3 rep4 Kdm4d-rep1 rep2
wdir=~/workspace/8.NT-HiC/5.maps/4.except/b.e2cell
ddir=~/workspace/8.NT-HiC/2.merge_sample/1.sequence/1.allValidPairs
for k in e2cell kdm4d-e2cell; do 
mkd $wdir/0.data/$k;done
cc_rep1_seq1._mm10.bwt2pairs.validPairs
ln -s $ddir/e2cell_rep3_seq1_allValidPairs $wdir/0.data/e2cell/e2cell_rep3._mm10.bwt2pairs.validPairs
ln -s $ddir/e2cell_rep4_seq1_allValidPairs $wdir/0.data/e2cell/e2cell_rep4._mm10.bwt2pairs.validPairs
ln -s $ddir/kdm4d-e2cell_rep1_seq1_allValidPairs $wdir/0.data/kdm4d-e2cell/kdm4d-e2cell_rep1._mm10.bwt2pairs.validPairs
ln -s $ddir/kdm4d-e2cell_rep2_seq1_allValidPairs $wdir/0.data/kdm4d-e2cell/kdm4d-e2cell_rep2._mm10.bwt2pairs.validPairs
mkd $wdir/logs
nohup HiC-Pro -i $wdir/0.data -o $wdir -c $ConfigHP/config-hicpro-p64.txt -s merge_persample -s build_contact_maps -s ice_norm > $wdir/logs/build.log &

ln -f ~/workspace/8.NT-HiC/5.maps/4.except/b.e2cell/hic_results/matrix/*/raw/*/*matrix 1.raw/
ln -f ~/workspace/8.NT-HiC/5.maps/4.except/b.e2cell/hic_results/matrix/*/iced/*/*matrix 2.iced/
ln -f ~/workspace/8.NT-HiC/5.maps/4.except/b.e2cell/hic_results/matrix/*/iced/*/*biases 3.biases/

ddir=~/workspace/8.NT-HiC/5.maps/4.except/1.raw
wdir=~/workspace/8.NT-HiC/5.maps/4.except
mv $wdir/depth.sorted.tab $wdir/OLD_depth/depth.$(date +%Y-%m-%d).tab
rm $wdir/depth.raw.tab
for i in $wdir/1.raw/*_2000000.matrix; do k=$(basename $i _2000000.matrix)
awk -v OFS="\t" -v key=$k '{depth+=$3} END{print key,depth}' $i >> $wdir/depth.raw.tab &
done;wait
sort -k1,1 $wdir/depth.raw.tab | awk '{total+=$2;n+=1;print} END{print "total",total;printf("%s\t%d\n","mean",total/n)}' > $wdir/depth.sorted.tab

wdir=~/workspace/8.NT-HiC/5.maps/4.except
f=100000000; odir=4.norm_to_100M
mkd $wdir/$odir
for res in 40000 100000; do 
for i in $wdir/2.iced/*_${res}_iced.matrix; do key=$(basename $i _${res}_iced.matrix)
depth=`awk -v key=$key '{if(key==$1)print $2}' $wdir/depth.sorted.tab`
awk -v f=$f -v depth=$depth '{print $1,$2,$3*f/depth}' $i > $wdir/$odir/${i##*/} &
echo "normalize: "$i
done;wait
done




