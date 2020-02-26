wdir=~/workspace/8.NT-HiC/k.fithic/1.all_except_10k_20180614
ddir=~/workspace/8.NT-HiC/5.maps/4.except/6.norm_raw_to_mean
bdir=~/workspace/8.NT-HiC/5.maps/4.except/3.biases
res=10000
for key in ${keys[@]}; do
(mkdir -p $wdir/$key/data $wdir/$key/output $wdir/$key/logs
nohup python ~/1.HiC-software/3.HiC-Pro_untils/hicpro2fithic.py -i $ddir/${key}_${res}.matrix -s $bdir/${key}_${res}_iced.matrix.biases \
-b $ConfigHP/${res}_mm10.bed -o $wdir/$key/data -r ${res} > $wdir/$key/logs/hicpro2fithic.log )&
done
##fithic占用内存很大，10k的fithic做不了
#应该normalize


wdir=~/workspace/8.NT-HiC/k.fithic/2.all_except_100k_20180614
ddir=~/workspace/8.NT-HiC/5.maps/4.except/1.raw
bdir=~/workspace/8.NT-HiC/5.maps/4.except/3.biases
res=100000
for key in ${keys[@]}; do
(mkdir -p $wdir/$key/data $wdir/$key/output $wdir/$key/logs
nohup python ~/1.HiC-software/3.HiC-Pro_untils/hicpro2fithic.py -i $ddir/${key}_${res}.matrix -s $bdir/${key}_${res}_iced.matrix.biases \
-b $ConfigHP/${res}_mm10.bed -o $wdir/$key/data -r ${res} > $wdir/$key/logs/hicpro2fithic.log 
nohup fithic -f $wdir/$key/data/fithic.fragmentMappability.gz -i $wdir/$key/data/fithic.interactionCounts.gz -t $wdir/$key/data/fithic.biases.gz \
-l $key -o $wdir/$key/output -v > $wdir/$key/logs/fithic.log )&
done


wdir=~/workspace/8.NT-HiC/k.fithic/3.all_except_40k_20180614
ddir=~/workspace/8.NT-HiC/5.maps/4.except/1.raw
bdir=~/workspace/8.NT-HiC/5.maps/4.except/3.biases
res=40000
for key in ${keys[@]}; do
mkdir -p $wdir/$key/data $wdir/$key/output $wdir/$key/logs
nohup python ~/1.HiC-software/3.HiC-Pro_untils/hicpro2fithic.py -i $ddir/${key}_${res}.matrix -s $bdir/${key}_${res}_iced.matrix.biases \
-b $ConfigHP/${res}_mm10.bed -o $wdir/$key/data -r ${res} > $wdir/$key/logs/hicpro2fithic.log 
nohup fithic -f $wdir/$key/data/fithic.fragmentMappability.gz -i $wdir/$key/data/fithic.interactionCounts.gz -t $wdir/$key/data/fithic.biases.gz \
-l $key -o $wdir/$key/output -v > $wdir/$key/logs/fithic.log
done

#####
zcat fithic.interactionCounts.gz | awk '$1=="chr1"&&$2==9305000{n+=$5} $3=="chr1"&&$4==9305000{n+=$5} END{print n}'
705
is same as fithic.fragmentMappability.gz
chr1    9300000 9305000 705     1


#########################to HiCPro
wdir=~/workspace/8.NT-HiC/k.fithic/2.all_except_100k_20180614/
for i in $wdir/1.output/*spline_pass1.significances.txt.gz; do k=$(basename $i .spline_pass1.significances.txt.gz)
python ~/codes/fithic2hicpro.py -i $i -o $wdir/3.toHicpro/$k.matrix -r 100000 &
done

for i in $wdir/3.toHicpro/*matrix; do 
cut -f 1-3 $i > ${i/.matrix/_100000_iced.matrix} &
done

key=3.qvalue_0.05
mkdir -p $wdir/3.toHicpro/$key
for i in $wdir/3.toHicpro/*matrix; do k=$(basename $i .matrix)
awk '$5<0.05{print $1,$2,$3}' $i > $wdir/3.toHicpro/$key/${k}_100000_iced.matrix &
done
wait
bash ~/codes/hic_plotter.sh -d $wdir/3.toHicpro/$key -w $wdir/3.toHicpro/$key -c chr1 -s 100 -e 500 -r 100000 -ext png &
bash ~/codes/hic_plotter.sh -d $wdir/3.toHicpro/$key -w $wdir/3.toHicpro/$key -c chr7 -r 100000 -ext png &

####################不要TAD plot is
bash ~/codes/hic_plotter.sh -d $wdir/3.toHicpro/1.all -w $wdir/3.toHicpro/1.all -c chr1 -s 100 -e 500 -r 100000 -ext png &
bash ~/codes/hic_plotter.sh -d $wdir/3.toHicpro/1.all -w $wdir/3.toHicpro/1.all -c chr7 -r 100000 -ext png &


#############################################################
#XieWei +NT
###################################################################
wdir=~/workspace/8.NT-HiC/k.fithic/3.all_except_40k_20181127
ddir=~/workspace/8.NT-HiC/5.maps/4.except/1.raw
bdir=~/workspace/8.NT-HiC/5.maps/4.except/3.biases
res=40000
key=e2cell
(mkdir -p $wdir/$key/data $wdir/$key/output $wdir/$key/logs
nohup python ~/1.HiC-software/3.HiC-Pro_untils/hicpro2fithic.py -i $ddir/${key}_${res}.matrix -s $bdir/${key}_${res}_iced.matrix.biases \
-b $ConfigHP/${res}_mm10.bed -o $wdir/$key/data -r ${res} > $wdir/$key/logs/hicpro2fithic.log 
nohup fithic -f $wdir/$key/data/fithic.fragmentMappability.gz -i $wdir/$key/data/fithic.interactionCounts.gz -t $wdir/$key/data/fithic.biases.gz -p 3 \
-l $key -o $wdir/$key/output -v > $wdir/$key/logs/fithic.log )&

wdir=~/workspace/8.NT-HiC/k.fithic/4.xw_40k_20181127
ddir=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/3.raw
bdir=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/5.biases
res=40000
key=8cell
(mkdir -p $wdir/$key/data $wdir/$key/output $wdir/$key/logs
nohup python ~/1.HiC-software/3.HiC-Pro_untils/hicpro2fithic.py -i $ddir/${key}_${res}.matrix -s $bdir/${key}_${res}_iced.matrix.biases \
-b $ConfigHP/${res}_mm10.bed -o $wdir/$key/data -r ${res} > $wdir/$key/logs/hicpro2fithic.log 
nohup fithic -f $wdir/$key/data/fithic.fragmentMappability.gz -i $wdir/$key/data/fithic.interactionCounts.gz -t $wdir/$key/data/fithic.biases.gz -p 3 \
-l $key -o $wdir/$key/output -v > $wdir/$key/logs/fithic.log )&
#do early2cell ICM icm e2cell PN5 late_2cell 8cell

##############500k
wdir=~/workspace/8.NT-HiC/k.fithic/5.all_except_500k_20181127
ddir=~/workspace/8.NT-HiC/5.maps/4.except/1.raw
bdir=~/workspace/8.NT-HiC/5.maps/4.except/3.biases
res=500000
key=e2cell
(mkdir -p $wdir/$key/data $wdir/$key/output $wdir/$key/logs
nohup python ~/1.HiC-software/3.HiC-Pro_untils/hicpro2fithic.py -i $ddir/${key}_${res}.matrix -s $bdir/${key}_${res}_iced.matrix.biases \
-b $ConfigHP/${res}_mm10.bed -o $wdir/$key/data -r ${res} > $wdir/$key/logs/hicpro2fithic.log 
nohup fithic -f $wdir/$key/data/fithic.fragmentMappability.gz -i $wdir/$key/data/fithic.interactionCounts.gz -t $wdir/$key/data/fithic.biases.gz -p 3 \
-l $key -o $wdir/$key/output -v > $wdir/$key/logs/fithic.log )&

wdir=~/workspace/8.NT-HiC/k.fithic/6.xw_500k_20181127
ddir=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/3.raw
bdir=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/5.biases
res=500000
key=early_2cell
(mkdir -p $wdir/$key/data $wdir/$key/output $wdir/$key/logs
nohup python ~/1.HiC-software/3.HiC-Pro_untils/hicpro2fithic.py -i $ddir/${key}_${res}.matrix -s $bdir/${key}_${res}_iced.matrix.biases \
-b $ConfigHP/${res}_mm10.bed -o $wdir/$key/data -r ${res} > $wdir/$key/logs/hicpro2fithic.log 
nohup fithic -f $wdir/$key/data/fithic.fragmentMappability.gz -i $wdir/$key/data/fithic.interactionCounts.gz -t $wdir/$key/data/fithic.biases.gz -p 3 \
-l $key -o $wdir/$key/output -v > $wdir/$key/logs/fithic.log )&

##################################################################################################
##############100k
##################################################################################################
wdir=~/workspace/8.NT-HiC/k.fithic/1.all_except_100k_20181202
ddir=~/workspace/8.NT-HiC/5.maps/4.except/1.raw
bdir=~/workspace/8.NT-HiC/5.maps/4.except/3.biases
res=100000
key=kdm4d-e2cell
mkdir -p $wdir/$key/data $wdir/$key/output $wdir/$key/logs $wdir/inter
(nohup python ~/1.HiC-software/3.HiC-Pro_untils/hicpro2fithic.py -i $ddir/${key}_${res}.matrix -s $bdir/${key}_${res}_iced.matrix.biases \
-b $ConfigHP/${res}_mm10.bed -o $wdir/$key/data -r ${res} > $wdir/$key/logs/hicpro2fithic.log 
nohup fithic -f $wdir/$key/data/fithic.fragmentMappability.gz -i $wdir/$key/data/fithic.interactionCounts.gz -t $wdir/$key/data/fithic.biases.gz -p 3 \
-l $key -o $wdir/$key/output -v > $wdir/$key/logs/fithic.log) &
#mkdir -p $wdir/inter
#nohup fithic -f $wdir/$key/data/fithic.fragmentMappability.gz -i $wdir/$key/data/fithic.interactionCounts.gz -t $wdir/$key/data/fithic.biases.gz -p 3 \
#-l $key -o $wdir/$key/inter -x interOnly > $wdir/$key/logs/inter.log &

wdir=~/workspace/8.NT-HiC/k.fithic/2.xw_100k_20181202
ddir=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/3.raw
bdir=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/5.biases
res=100000
key=MII_oocyte
(mkdir -p $wdir/$key/data $wdir/$key/output $wdir/$key/logs
nohup python ~/1.HiC-software/3.HiC-Pro_untils/hicpro2fithic.py -i $ddir/${key}_${res}.matrix -s $bdir/${key}_${res}_iced.matrix.biases \
-b $ConfigHP/${res}_mm10.bed -o $wdir/$key/data -r ${res} > $wdir/$key/logs/hicpro2fithic.log 
nohup fithic -f $wdir/$key/data/fithic.fragmentMappability.gz -i $wdir/$key/data/fithic.interactionCounts.gz -t $wdir/$key/data/fithic.biases.gz -p 3 -l $key -o $wdir/$key/output -v > $wdir/$key/logs/fithic.log )&
mkdir -p $wdir/inter
nohup fithic -f $wdir/$key/data/fithic.fragmentMappability.gz -i $wdir/$key/data/fithic.interactionCounts.gz -t $wdir/$key/data/fithic.biases.gz -p 3 -l $key -o $wdir/$key/inter -x interOnly > $wdir/$key/logs/inter.log &

for i in ~/workspace/8.NT-HiC/k.fithic/[12]*/*/output/*3.significances.txt.gz; do 
echo "zcat $i > ${i/.gz}"
done




###################10kb
wdir=~/workspace/8.NT-HiC/k.fithic/7.xw_10k_20181213
ddir=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/3.raw
bdir=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/5.biases
res=10000
key=ICM
mkdir -p $wdir/$key/data $wdir/$key/output $wdir/$key/logs
nohup python ~/1.HiC-software/3.HiC-Pro_untils/hicpro2fithic.py -i $ddir/${key}_${res}.matrix -s $bdir/${key}_${res}_iced.matrix.biases \
-b $ConfigHP/${res}_mm10.bed -o $wdir/$key/data -r ${res} > $wdir/$key/logs/hicpro2fithic.log &
zcat fithic.biases.gz | awk '{print > "biases."$1}' &
zcat fithic.fragmentMappability.gz | awk '{print > "fragmentMappability."$1}' &
zcat fithic.interactionCounts.gz | awk '$1==$3{print > "interactionCounts."$1}' &
for i in $(seq 5 19) X; do chr=chr$i; 
nohup fithic -f $wdir/$key/data/fragmentMappability.$chr.gz -i $wdir/$key/data/interactionCounts.$chr.gz -t $wdir/$key/data/biases.$chr.gz -p 3 -l $chr -o $wdir/$key/output -v > $wdir/$key/logs/fithic.$chr.log ; done &

wdir=~/workspace/8.NT-HiC/k.fithic/8.except_10k_20181213
ddir=~/workspace/8.NT-HiC/5.maps/4.except/1.raw
bdir=~/workspace/8.NT-HiC/5.maps/4.except/3.biases
res=10000
key=icm
mkdir -p $wdir/$key/data $wdir/$key/output $wdir/$key/logs
nohup python ~/1.HiC-software/3.HiC-Pro_untils/hicpro2fithic.py -i $ddir/${key}_${res}.matrix -s $bdir/${key}_${res}_iced.matrix.biases \
-b $ConfigHP/${res}_mm10.bed -o $wdir/$key/data -r ${res} > $wdir/$key/logs/hicpro2fithic.log &
zcat fithic.biases.gz | awk '{print > "biases."$1}' &
zcat fithic.fragmentMappability.gz | awk '{print > "fragmentMappability."$1}' &
zcat fithic.interactionCounts.gz | awk '$1==$3{print > "interactionCounts."$1}' &
for i in $(seq 1 19) X; do chr=chr$i; 
nohup fithic -f $wdir/$key/data/fragmentMappability.$chr.gz -i $wdir/$key/data/interactionCounts.$chr.gz -t $wdir/$key/data/biases.$chr.gz -p 3 -l $chr -o $wdir/$key/output -v > $wdir/$key/logs/fithic.$chr.log ; done &

#############################
#合并normalized interaction frequency 到fithic output
ddir=~/workspace/8.NT-HiC/5.maps/4.except/8.chroms
wdir=~/workspace/8.NT-HiC/k.fithic/f.merge_100k
ln -s ~/workspace/8.NT-HiC/k.fithic/1*/e*/output/*pass3.significances.txt NT-e2cell.fithic
ln -s ~/workspace/8.NT-HiC/k.fithic/1*/icm/output/*pass3.significances.txt NT-ICM.fithic
ln -s ~/workspace/8.NT-HiC/k.fithic/2*/e*/output/*pass3.significances.txt NF-e2cell.fithic
ln -s ~/workspace/8.NT-HiC/k.fithic/2*/ICM/output/*pass3.significances.txt NF-ICM.fithic
ln -s ~/workspace/8.NT-HiC/k.fithic/1*/l2cell/output/*pass3.significances.txt NT-l2cell.fithic
ln -s ~/workspace/8.NT-HiC/k.fithic/1*/8cell/output/*pass3.significances.txt NT-8cell.fithic
ln -s ~/workspace/8.NT-HiC/k.fithic/2*/late_2cell/output/*pass3.significances.txt NF-l2cell.fithic
ln -s ~/workspace/8.NT-HiC/k.fithic/2*/8cell/output/*pass3.significances.txt NF-8cell.fithic
ln -s ~/workspace/8.NT-HiC/k.fithic/1*/cc/output/*pass3.significances.txt CC.fithic
ln -s ~/workspace/8.NT-HiC/k.fithic/2*/Sperm/output/*pass3.significances.txt Sperm.fithic
ln -s ~/workspace/8.NT-HiC/k.fithic/2*/MII*/output/*pass3.significances.txt MII-oocyte.fithic ##

res=100000
for j in Sperm; do 
for i in $(seq 1 19) X; do chr=chr$i
awk -v res=$res -v chr=$chr 'ARGIND==1&&$1==res&&$2==chr{start=$3} 
ARGIND==2{print chr,($1-start+0.5)*res,chr,($2-start+0.5)*res,$3}' $ConfigHP/HP_chrom_mm10.tab $ddir/${j}_${res}_iced.$chr.matrix > $wdir/${j}.$chr.mat &
done
done

for j in Sperm; do 
cat $wdir/${j}*mat | awk 'NR==FNR{a[$1"\t"$2"\t"$3"\t"$4]=$5} NR>FNR&&FNR>1{print $0,a[$1"\t"$2"\t"$3"\t"$4]}' - $wdir/${j}.fithic > $wdir/${j}.pass3.txt &
done





##########################################################
#先测试pass2和pass3的qvalue<0.01 的NT和NF的交集数量
#在非交集的peaks中，分别取等数量最显著的peaks
#以这些peaks来作为各自的显著peaks
#########################################
wdir=~/workspace/8.NT-HiC/k.fithic/c.pass_overlap_test/2.pass3
res=100000


zcat ~/workspace/8.NT-HiC/k.fithic/2*/ICM/output/*pass3.significances.txt.gz | awk -v res=$res '$7<0.0001&&($4-$2)>res' > NF-ICM.peaks &
zcat ~/workspace/8.NT-HiC/k.fithic/1*/icm/output/*pass3.significances.txt.gz | awk -v res=$res '$7<0.0001&&($4-$2)>res' > NT-ICM.peaks &
zcat ~/workspace/8.NT-HiC/k.fithic/2*/e*/output/*pass3.significances.txt.gz | awk -v OFS="\t" -v res=$res '$7<0.0001&&($4-$2)>res' > NF-e2cell.peaks &
zcat ~/workspace/8.NT-HiC/k.fithic/1*/e*/output/*pass3.significances.txt.gz | awk -v OFS="\t" -v res=$res '$7<0.0001&&($4-$2)>res' > NT-e2cell.peaks &
zcat ~/workspace/8.NT-HiC/k.fithic/1*/l2*/output/*pass3.significances.txt.gz | awk -v OFS="\t" -v res=$res '$7<0.0001&&($4-$2)>res' > NT-l2cell.peaks &
zcat ~/workspace/8.NT-HiC/k.fithic/2*/late*/output/*pass3.significances.txt.gz | awk -v OFS="\t" -v res=$res '$7<0.0001&&($4-$2)>res' > NF-l2cell.peaks &
zcat ~/workspace/8.NT-HiC/k.fithic/2*/8*/output/*pass3.significances.txt.gz | awk -v res=$res '$7<0.0001&&($4-$2)>res' > NF-8cell.peaks &
zcat ~/workspace/8.NT-HiC/k.fithic/1*/8*/output/*pass3.significances.txt.gz | awk -v res=$res '$7<0.0001&&($4-$2)>res' > NT-8cell.peaks &
zcat ~/workspace/8.NT-HiC/k.fithic/2*/PN3*/output/*pass3.significances.txt.gz | awk -v res=$res '$7<0.0001&&($4-$2)>res' > NF-PN3.peaks &
zcat ~/workspace/8.NT-HiC/k.fithic/1*/6h*/output/*pass3.significances.txt.gz | awk -v res=$res '$7<0.0001&&($4-$2)>res' > NT-PN3.peaks &
zcat ~/workspace/8.NT-HiC/k.fithic/1*/cc*/output/*pass3.significances.txt.gz | awk -v res=$res '$7<0.0001&&($4-$2)>res' > CC.peaks &


for i in e2cell l2cell 8cell ICM; do 
a=`cat NF-$i.peaks | wc -l`
b=`cat NT-$i.peaks | wc -l`
c=`cat *$i.peaks | cut -f 1-4 | sort | uniq -c | awk '$1==2' | wc -l`
awk ''
echo -e $i"\t"$a"\t"$b"\t"$c
done

for i in e2cell l2cell 8cell ICM; do 
awk -v samefile=NT-$i.same.peaks -v difffile=NT-$i.diff.temp 'NR==FNR{a[$1"\t"$2"\t"$3"\t"$4]} NR>FNR{if($1"\t"$2"\t"$3"\t"$4 in a){print $0 > samefile}else{print $0 > difffile}}' NF-$i.peaks NT-$i.peaks &
awk -v samefile=NF-$i.same.peaks -v difffile=NF-$i.diff.temp 'NR==FNR{a[$1"\t"$2"\t"$3"\t"$4]} NR>FNR{if($1"\t"$2"\t"$3"\t"$4 in a){print $0 > samefile}else{print $0 > difffile}}' NT-$i.peaks NF-$i.peaks &
done

#固定数量不太好，取到最低的那个NT或NF的总数
for i in e2cell:17831 l2cell:23700 8cell:8902 ICM:86354; do 
for j in NF NT; do k=$j-${i%%:*};d=${i##*:}
(sort -k7g,7 $k.diff.temp | head -$d > $k.diff.peaks 
cat $k.same.peaks $k.diff.peaks > $k.signi.peaks)&
done; done 

#inter
for i in e2cell:48839 l2cell:90702 8cell:91025 ICM:101089; do 
for j in NF NT; do k=$j-${i%%:*};d=${i##*:}
sort -k7g,7 $k.diff.temp | head -$d | cat $k.same.peaks - > $k.signi.peaks &
done; done 

####signi.peaks 在TAD中的数量
for j in NF NT; do 
for i in e2cell l2cell 8cell ICM; do k=${j}-$i.signi.peaks
a=`awk '{print $1,$2+50000,$4-50000}' $k | intersectBed -a - -b ~/workspace/8.NT-HiC/g.DI_ALL/h.delta_RTI_reps/3.bed -f 1 -u | wc -l`
b=`awk '{print $1,$2+50000,$4-50000}' $k | intersectBed -a - -b ~/workspace/8.NT-HiC/g.DI_ALL/h.delta_RTI_reps/4.bed -f 1 -u | wc -l`
c=`awk '{print $1,$2+50000,$4-50000}' $k | intersectBed -a - -b ~/workspace/8.NT-HiC/g.DI_ALL/h.delta_RTI_reps/all.RTI.bed -f 1 -u | wc -l`
d=`awk '{print $1,$2+50000,$4-50000}' $k | wc -l `
echo -e ${j}"\t"$i"\t"$a"\t"$b"\t"$c"\t"$d
done;done

############################
#
wdir=~/workspace/8.NT-HiC/k.fithic/c.pass_overlap_test/6.normalized
ddir=~/workspace/8.NT-HiC/k.fithic/f.merge_100k
res=100000
for i in NT NF; do for j in e2cell l2cell 8cell ICM; do 
awk -v res=$res '$7<0.01&&($4-$2)>res&&$8>15' $ddir/$i-$j.pass3.txt > $wdir/$i-$j.peaks &
done;done




####
peaks=~/workspace/8.NT-HiC/k.fithic/c.pass_overlap_test/6.normalized
wdir=~/workspace/8.NT-HiC/k.fithic/c.pass_overlap_test/6.normalized/test_in_SE
k27ac=~/workspace/9.NT-ChIP/2.public/a.Renbin/5.homer/1.bed/2cell_H3K27ac.super5.bed
f=NT-e2cell
awk -v f=$f '{n+=1;print $1,$2-50000,$2+50000,f"-"n > f".left.sequence";print $3,$4-50000,$4+50000,f"-"n > f".right.sequence"}' $peaks/$f.signi.peaks
cat $f.left.sequence $f.right.sequence | intersectBed -a - -b $k27ac -wa -u | cut -f 4 | sort | uniq -c | awk '{print $1}' | sort | uniq -c

 
##############make like-fithic file from Enhancer and promoter
edir=~/workspace/9.NT-ChIP/2.public/a.Renbin/5.homer/1.bed
promoter=~/ann/mm10_promoter_1000-1000.bed
wdir=~/workspace/8.NT-HiC/k.fithic/h.EP_fithic-like-file2
ddir=~/workspace/8.NT-HiC/2.merge_sample/4.except/3.allValidPairs_byChr
mkdir -p $wdir
res=100000
for e in $edir/*super5.bed; do k=$(basename $e _H3K27ac.super5.bed)
awk -v res=$res '{a=int($2/res);b=int($3/res)} NR==FNR{for(i=a;i<=b;i++)en[$1"\t"i]} NR>FNR{for(i=a;i<=b;i++){for(j in en){print $1,i,j}}}' $e $promoter | \
awk -v res=$res '{if($2>$4){i=$2;$2=$4;$4=i}} $1==$3{print $1,($2+0.5)*res,$3,($4+0.5)*res,1}' | \
sort -k1,1 -k2n,2 -k4n,4 | uniq > $wdir/$k.SEP.peaks &
done
wait

for f in e2cell:2cell l2cell:2cell 8cell:8cell ICM:ES; do k=${f%%:*};p=${f##*:};mkdir -p $wdir/$k
for i in $(seq 1 19) X; do chr=chr$i;for j in NF NT; do 
python ~/codes/contacts_in_loops_from_validPairs.py -i $ddir/$j-$k.$chr.Pairs -f $wdir/$p.SEP.peaks -r $res -c $chr -o $wdir/$k/$j.$chr.Pairs &
done;done;wait
done

for f in e2cell:2cell l2cell:2cell 8cell:8cell ICM:ES; do k=${f%%:*};p=${f##*:};cd $wdir/$k
k27ac=~/workspace/9.NT-ChIP/2.public/a.Renbin/5.homer/1.bed/${p}_H3K27ac.super5.bed
for j in NF NT; do
awk 'FNR>1{if($3-$2<500){if($4=="+"){$2=$3-500}else{$3=$2+500}} if($7-$6<500){if($8=="+"){$6=$7-500}else{$7=$6+500}};print $1,$2,$3,$5,$6,$7,$9;}' $j.chr*Pairs > $j.all.Pairs 
awk '{print $1,$2,$3,$7;print $4,$5,$6,$7}' $j.all.Pairs >$j.all.sequences 
intersectBed -a $j.all.sequences -b $promoter -r -e -f 0.5 -wa -u > temp$j &
intersectBed -a $j.all.sequences -b $k27ac -r -e -f 0.5 -wa -u > temp3$j &
done
wait
awk 'FILENAME=="tempNF"{a[$0]} FILENAME=="temp3NF"{b[$0]} FILENAME=="NF.all.Pairs"{left=$1"\t"$2"\t"$3"\t"$7;right=$4"\t"$5"\t"$6"\t"$7;
if((left in a)&&(right in b)){print left,right > "NF.ep.pairs"};
if((left in b)&&(right in a)){print right,left > "NF.ep.pairs"}
if((left in a)&&(right in a)){print left,right > "NF.pp.pairs";print right,left > "NF.pp.pairs"};
if((left in a)&&!(right in b)&&!(right in a)){print left,right > "NF.np.pairs"};
if((right in a)&&!(left in b)&&!(left in a)){print right,left > "NF.np.pairs"}
if(!(left in a)&&!(right in b)&&!(left in b)&&!(right in a)){print left,right > "NF.null.pairs"}
if(!((left in a)&&(right in b))&&!((left in b)&&(right in a))){print left,right > "NF.other.pairs"}
}' tempNF temp3NF NF.all.Pairs &
awk 'FILENAME=="tempNT"{a[$0]} FILENAME=="temp3NT"{b[$0]} FILENAME=="NT.all.Pairs"{left=$1"\t"$2"\t"$3"\t"$7;right=$4"\t"$5"\t"$6"\t"$7;
if((left in a)&&(right in b)){print left,right > "NT.ep.pairs"};
if((left in b)&&(right in a)){print right,left > "NT.ep.pairs"}
if((left in a)&&(right in a)){print left,right > "NT.pp.pairs";print right,left > "NT.pp.pairs"};
if((left in a)&&!(right in b)&&!(right in a)){print left,right > "NT.np.pairs"};
if((right in a)&&!(left in b)&&!(left in a)){print right,left > "NT.np.pairs"}
if(!(left in a)&&!(right in b)&&!(left in b)&&!(right in a)){print left,right > "NT.null.pairs"}
if(!((left in a)&&(right in b))&&!((left in b)&&(right in a))){print left,right > "NT.other.pairs"}
}' tempNT temp3NT NT.all.Pairs &
wait
done
mkdir -p $wdir/figures

#########random Enhancer
for e in ~/workspace/9.NT-ChIP/2.public/a.Renbin/5.homer/1.bed/*super5.bed; do k=$(basename $e)
shuffleBed -i $e -g $ConfigHP/chrom_mm10.sizes -seed 100 > ~/workspace/9.NT-ChIP/2.public/a.Renbin/5.homer/1.bed/random/$k &
done
edir=~/workspace/9.NT-ChIP/2.public/a.Renbin/5.homer/1.bed/random
promoter=~/ann/mm10_promoter_1000-1000.bed
wdir=~/workspace/8.NT-HiC/k.fithic/i.EP_fithic-like-file_random
ddir=~/workspace/8.NT-HiC/2.merge_sample/4.except/3.allValidPairs_byChr
mkdir -p $wdir
res=100000


########################################
#
wdir=~/workspace/8.NT-HiC/k.fithic/g.contacts_distance_density
ddir=~/workspace/8.NT-HiC/2.merge_sample/4.except/3.allValidPairs_byChr
for f in e2cell:2cell l2cell:2cell 8cell:8cell ICM:ES; do k=${f%%:*};p=${f##*:};mkdir -p $wdir/$k
for i in $(seq 1 19) X; do chr=chr$i;for j in NF NT; do 
python ~/codes/contacts_from_validPairs.py -i $ddir/$j-$k.$chr.Pairs -f 500 -r $wdir/$k/$j.$chr.dist -c $chr -o $wdir/$k/$j.$chr.Pairs &
done;done;wait
done

k=1cell;mkdir -p $wdir/$k
for j in T-1h A-12h A-6h A-1h PN5-zygote PN3-zygote MII-oocyte; do
for i in $(seq 1 19) X; do chr=chr$i;
python ~/codes/contacts_from_validPairs.py -i $ddir/$j.$chr.Pairs -f 500 -r $wdir/$k/$j.$chr.dist -c $chr -o $wdir/$k/$j.$chr.Pairs &
done;done;wait

k=CC;mkdir -p $wdir/$k
for j in T-30min CC; do
for i in $(seq 1 19) X; do chr=chr$i;
python ~/codes/contacts_from_validPairs.py -i $ddir/$j.$chr.Pairs -f 500 -r $wdir/$k/$j.$chr.dist -c $chr -o $wdir/$k/$j.$chr.Pairs &
done;done;wait

