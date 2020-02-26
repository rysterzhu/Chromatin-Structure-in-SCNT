#按lj的方法，将rearrangment (retain2) 拆出了lost
keys=(cc 05h 1h 2h 3h 6h 12h e2cell l2cell 4cell 8cell morula icm te)
ddir1=~/workspace/8.NT-HiC/g.DI_ALL/9.except_norm_mean_raw_20180713/2.unoverlap_tads
wdir=~/workspace/8.NT-HiC/g.DI_ALL/b.tad_shift_analysis/1.NT_all
A=~/workspace/8.NT-HiC/g.DI_ALL/9.except_norm_mean_raw_20180713/2.unoverlap_tads/cc.tad
(
echo -e " \tstable\tmerge\tsplit\trearrangment\tgain\tlost\ttotal"
for i in $(seq 1 13); do B=${keys[i]}.tad;
intersectBed -b $A -a $ddir1/$B -F 1.0 -wa | sort | uniq -c | awk '$1>1{print $2,$3,$4}' >  $wdir/${B/tad/temp} #
intersectBed -b $A -a $ddir1/$B -F 1.0 -wo | awk 'NR==FNR{a[$0]} NR!=FNR{if($1"\t"$2"\t"$3 in a){print $4,$5,$6}}' $wdir/${B/tad/temp} - > $wdir/${B/tad/merge}
intersectBed -a $A -b $ddir1/$B -F 1.0 -wa | sort | uniq -c | awk '$1>1{print $2,$3,$4}' >  $wdir/${B/tad/split}
cat $wdir/${B/tad/merge} $wdir/${B/tad/split} | awk 'NR==FNR{a[$0]} NR!=FNR{if(!($0 in a)){print}}' - $A > $wdir/${B/tad/retain}
intersectBed -a $wdir/${B/tad/retain} -b $ddir1/$B -r -f 0.75 -wa > $wdir/${B/tad/stable}
intersectBed -a $wdir/${B/tad/retain} -b $ddir1/$B -r -f 0.75 -v -wa > $wdir/${B/tad/retain2}
intersectBed -a $wdir/${B/tad/retain2} -b $ddir1/$B -wa -u > $wdir/${B/tad/rearrangment}
#intersectBed -a $wdir/${B/tad/retain2} -b $ddir1/$B -v -wa > $wdir/${B/tad/lost}
intersectBed -a $A -b $ddir1/$B -v -wa > $wdir/${B/tad/lost}
intersectBed -a $ddir1/$B -b $A -v -wa > $wdir/${B/tad/gain}
stable=`cat $wdir/${B/tad/stable} | wc -l `
merge=`cat $wdir/${B/tad/merge} | wc -l `
split=`cat $wdir/${B/tad/split} | wc -l `
rearrangment=`cat $wdir/${B/tad/rearrangment} | wc -l `
gain=`cat $wdir/${B/tad/gain} | wc -l `
lost=`cat $wdir/${B/tad/lost} | wc -l `
total=`cat $A | wc -l`
echo -e ${keys[i]}"\t"$stable"\t"$merge"\t"$split"\t"$rearrangment"\t"$gain"\t"$lost"\t"$total
done
)

for i in $wdir/*; do flag=${i##*.}
sed -i "s/$/\t${flag}/" $i;
done

for i in $(seq 1 13); do B=${keys[i]}
cat $B.stable $B.merge $B.split $B.retain2 | sort -k1V,1 -k2n,2 > $B.all
done

for i in $(seq 1 13); do B=${keys[i]}
cat $B.stable $B.merge $B.split $B.gain $B.lost $B.rearrangment | sort -k1V,1 -k2n,2 > $B.all2
done
#################################################
keys=(12h e2cell l2cell 8cell icm)
ddir1=~/workspace/8.NT-HiC/g.DI_ALL/9.except_norm_mean_raw_20180713/2.unoverlap_tads
wdir=~/workspace/8.NT-HiC/g.DI_ALL/b.tad_shift_analysis/2.NT-NF/NT
A=~/workspace/8.NT-HiC/g.DI_ALL/9.except_norm_mean_raw_20180713/2.unoverlap_tads/cc.tad
mkdir -p $wdir

keys=(PN5_zygote early_2cell late_2cell 8cell ICM)
ddir1=~/workspace/8.NT-HiC/3.public_data/f.XW_analysis/5.DI/2.unoverlap_tads
wdir=~/workspace/8.NT-HiC/g.DI_ALL/b.tad_shift_analysis/2.NT-NF/NF
A=~/workspace/8.NT-HiC/g.DI_ALL/9.except_norm_mean_raw_20180713/2.unoverlap_tads/cc.tad
mkdir -p $wdir

(
echo -e " \tstable\tmerge\tsplit\trearrangment\tgain\tlost\ttotal"
for i in $(seq 0 4); do B=${keys[i]}.tad;
intersectBed -b $A -a $ddir1/$B -F 1.0 -wa | sort | uniq -c | awk '$1>1{print $2,$3,$4}' >  $wdir/${B/tad/temp} #
intersectBed -b $A -a $ddir1/$B -F 1.0 -wo | awk 'NR==FNR{a[$0]} NR!=FNR{if($1"\t"$2"\t"$3 in a){print $4,$5,$6}}' $wdir/${B/tad/temp} - > $wdir/${B/tad/merge}
intersectBed -a $A -b $ddir1/$B -F 1.0 -wa | sort | uniq -c | awk '$1>1{print $2,$3,$4}' >  $wdir/${B/tad/split}
cat $wdir/${B/tad/merge} $wdir/${B/tad/split} | awk 'NR==FNR{a[$0]} NR!=FNR{if(!($0 in a)){print}}' - $A > $wdir/${B/tad/retain}
intersectBed -a $wdir/${B/tad/retain} -b $ddir1/$B -r -f 0.75 -wa > $wdir/${B/tad/stable}
intersectBed -a $wdir/${B/tad/retain} -b $ddir1/$B -r -f 0.75 -v -wa > $wdir/${B/tad/retain2}
intersectBed -a $wdir/${B/tad/retain2} -b $ddir1/$B -wa -u > $wdir/${B/tad/rearrangment}
#intersectBed -a $wdir/${B/tad/retain2} -b $ddir1/$B -v -wa > $wdir/${B/tad/lost}
intersectBed -a $A -b $ddir1/$B -v -wa > $wdir/${B/tad/lost}
intersectBed -a $ddir1/$B -b $A -v -wa > $wdir/${B/tad/gain}
stable=`cat $wdir/${B/tad/stable} | wc -l `
merge=`cat $wdir/${B/tad/merge} | wc -l `
split=`cat $wdir/${B/tad/split} | wc -l `
rearrangment=`cat $wdir/${B/tad/rearrangment} | wc -l `
gain=`cat $wdir/${B/tad/gain} | wc -l `
lost=`cat $wdir/${B/tad/lost} | wc -l `
total=`cat $A | wc -l`
echo -e ${keys[i]}"\t"$stable"\t"$merge"\t"$split"\t"$rearrangment"\t"$gain"\t"$lost"\t"$total
done
)
#在各种类型的bed的第四列添加flag
for i in $wdir/*; do flag=${i##*.}
sed -i "s/$/\t${flag}/" $i;
done

# for i in $(seq 0 4); do B=$wdir/${keys[i]}
# cat $B.stable $B.merge $B.split $B.retain2 | sort -k1V,1 -k2n,2 > $B.all
# done
#合并各种类型到一个tad文件
for i in $(seq 0 4); do B=$wdir/${keys[i]}
cat $B.stable $B.merge $B.split $B.gain $B.lost $B.rearrangment | sort -k1V,1 -k2n,2 > $B.all2
done
#

#################################################
wdir=~/workspace/8.NT-HiC/g.DI_ALL/b.tad_shift_analysis/3.cmm
ddir=~/workspace/8.NT-HiC/g.DI_ALL/2.except_norm_100M_20180825/2.unoverlap_tads
cp $ddir/cc.tad $wdir/
cp $ddir/icm.tad $wdir/
cp ~/workspace/8.NT-HiC/3.public_data/f.XW_analysis/5.DI/2.unoverlap_tads/ICM.tad $wdir/

keys=(cc ICM)
ddir1=$wdir
A=$wdir/icm.tad
ratio1=0.75
ratio2=0.9
(
echo -e " \tstable\tmerge\tsplit\tshift\tgain\tlost\ttotal"
for i in ${keys[@]}; do B=$i.tad; #A到B的变化
intersectBed -b $A -a $ddir1/$B -F $ratio2 -wa | sort | uniq -c | awk '$1>1{print $2,$3,$4}' >  $wdir/${B/tad/temp} #把至少两个A在一个B里的B筛选出来得到temp
intersectBed -b $A -a $ddir1/$B -F $ratio2 -wo | awk 'NR==FNR{a[$0]} NR!=FNR{if($1"\t"$2"\t"$3 in a){print $4,$5,$6}}' $wdir/${B/tad/temp} - > $wdir/${B/tad/merge} #temp中存在的那些A即为merge
intersectBed -a $A -b $ddir1/$B -F $ratio2 -wa | sort | uniq -c | awk '$1>1{print $2,$3,$4}' >  $wdir/${B/tad/split} #把至少两个B在一个A里的A即为split
cat $wdir/${B/tad/merge} $wdir/${B/tad/split} | awk 'NR==FNR{a[$0]} NR!=FNR{if(!($0 in a)){print}}' - $A > $wdir/${B/tad/retain} #所有的A中去除merge和split剩余的得到retain
intersectBed -a $wdir/${B/tad/retain} -b $ddir1/$B -r -f $ratio1 -wa > $wdir/${B/tad/stable} #retain中与B交集超过ratio1的为stable
intersectBed -a $wdir/${B/tad/retain} -b $ddir1/$B -r -f $ratio1 -v -wa > $wdir/${B/tad/retain2} #剩余的得到retain2
intersectBed -a $wdir/${B/tad/retain2} -b $ddir1/$B -wa -u > $wdir/${B/tad/shift} #retain2中与B有交集的为shift（交集小于ratio1）
#intersectBed -a $wdir/${B/tad/retain2} -b $ddir1/$B -v -wa > $wdir/${B/tad/lost} #retain2中与B没有交集的A为lost
intersectBed -a $A -b $ddir1/$B -v -wa > $wdir/${B/tad/lost} #与B没有交集的A为lost
intersectBed -a $ddir1/$B -b $A -v -wa > $wdir/${B/tad/gain} #与A没有交集的B为gain
stable=`cat $wdir/${B/tad/stable} | wc -l `
merge=`cat $wdir/${B/tad/merge} | wc -l `
split=`cat $wdir/${B/tad/split} | wc -l `
shift=`cat $wdir/${B/tad/shift} | wc -l `
gain=`cat $wdir/${B/tad/gain} | wc -l `
lost=`cat $wdir/${B/tad/lost} | wc -l `
total=`cat $A | wc -l`
echo -e $i"\t"$stable"\t"$merge"\t"$split"\t"$shift"\t"$gain"\t"$lost"\t"$total
done
)
 	stable	merge	split	shift	gain	lost	total
cc	890	216	270	408	43	3	1787
ICM	1108	184	260	234	17	1	1787

flags=(stable merge split shift gain lost)
colors=("27,158,119" "217,95,2" "117,112,179" "231,41,138" "102,166,30" "230,171,2")
for i in ${keys[@]}; do rm $wdir/$i.all
for j in $(seq 0 5); do 
awk -v k=${flags[j]} -v c=${colors[j]} '{print $1,$2,$3,k,0,".",$2,$3,c}' $wdir/$i.${flags[j]} >> $wdir/$i.all
done
sort -k1V,1 -k2n,2 $wdir/$i.all > $wdir/$i.color.bed
cut -f 1-4 $wdir/$i.color.bed > $wdir/$i.all2
done

awk '{print $0,"stable"}' icm.tad > icm2.all2


flags=(stable merge split retain2)
colors=("255,0,0" "0,0,255" "0,255,0" "0,0,0")
for i in ${keys[@]}; do rm $wdir/$i.all
for j in $(seq 0 3); do 
awk -v k=${flags[j]} -v c=${colors[j]} '{print $1,$2,$3,k,0,".",$2,$3,c}' $wdir/$i.${flags[j]} >> $wdir/$i.retain2.all
done
sort -k1V,1 -k2n,2 $wdir/$i.retain2.all > $wdir/$i.color.retain2.bed
cut -f 1-4 $wdir/$i.color.retain2.bed > $wdir/$i.retain2.all2
done