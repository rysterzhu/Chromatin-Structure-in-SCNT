#!/bin/bash
res=40000

while [ -n "$1" ]; do
case "$1" in
    -h) help; shift 1;;
    -w) wdir=$2; shift 2;;
    -m) matrix=$2; shift 2;;
	-r) res=$2; shift 2;;
    -*) echo "error: no such option $1. -h for help";exit 1;;
     *) break;
esac
done

mkdir -p $wdir/0.raw $wdir/3.plot_tads $wdir/1.tads $wdir/4.directional_index $wdir/2.unoverlap_tads $wdir/logs

key=${matrix##*/}; key=${key/_${res}_iced.matrix}
###1.合并replicates（见HiCPro）
###2.转为di input格式
mkdir -p $wdir/0.raw/$key
cd $wdir/0.raw/$key
python ~/software/hic-pro/HiC-Pro_2.9.0/bin/utils/sparseToDense.py -o $key.mat -b $ConfigHP/${res}_mm10.bed $matrix --perchr -d  #必须得拆开做，不然太大了
#3.
for i in *$key.mat; do
perl ~/1.HiC-software/8.DI/1DI_from_matrix.pl $i ${res} $(( res * 50 )) $ConfigHP/chrom_mm10.sizes > ${i/mat/di}
done
#4.cat
cat *$key.di | awk -v OFS="\t" '$1~/[0-9]+/{$0=$0;  print} $1~"X"{$1=20;print} $1~"Y"{$1=21;print} ' | sort -k1n,1 -k2n,2 > $wdir/0.raw/${key}.di 
echo "calc di done: "$i
#5.matlab
cd $wdir
matlab -nodesktop -nosplash -r "inputfile='$wdir/0.raw/${key}.di',outputfile='$wdir/0.raw/${key}.hmm';HMM_calls;quit" > $wdir/logs/${key}.matlab.log

##6.
perl ~/1.HiC-software/8.DI/3file_ends_cleaner.pl $wdir/0.raw/${key}.hmm $wdir/0.raw/${key}.di | perl ~/1.HiC-software/8.DI/4converter_7col.pl > $wdir/0.raw/${key}.7col
#7.
perl ~/1.HiC-software/8.DI/5hmm_probablity_correcter.pl $wdir/0.raw/${key}.7col 2 0.99 $res | \
perl ~/1.HiC-software/8.DI/6hmm-state_caller.pl $ConfigHP/chrom_mm10.sizes chr | \
perl ~/1.HiC-software/8.DI/7hmm-state_domains.pl > $wdir/0.raw/${key}.raw.tad

#8.others
# 去除一些错误的TAD
awk -v OFS="\t" '{chrs[$1]}
$2>$3{t=$2;$2=$3;$3=t;} #将TAD的方向换成正的
$1=="chr20"{$1="chrX"}  #将TAD的chr20换成X
$1=="chr21"||$1==""||$2==""||$3==""{next}   #不要chrY,空行
{$0=$0;print}
' $wdir/0.raw/${key}.raw.tad | sort -k1V,1 -k2n,2 > $wdir/1.tads/${key}.tad

# 去除重叠的TAD中的大的，便于计算
awk -v OFS="\t" 'NR>1&&(a!=$1||c<=$2){print a,b,c} {a=$1;b=$2;c=$3} END{print a,b,c}' $wdir/1.tads/${key}.tad >  $wdir/2.unoverlap_tads/${key}.tad

# 为没有TAD的染色体添加一个假的TAD，用于画图
awk -v OFS="\t" -v res=$res '{chrs[$1];print}
END{for(i=1;i<=19;i++){
if(!("chr"i in chrs)){print "chr"i,0,res}};
if(!("chrX" in chrs)){print "chrX",0,res}
}
' $wdir/2.unoverlap_tads/${key}.tad | sort -k1V,1 -k2n,2 > $wdir/3.plot_tads/${key}.tad

# 将di转成bedGraph
awk -v OFS="\t" '
$1==20{$1="X"}
$1==21||$1==""||$2==""||$3==""||$4==""{next}   #不要chrY,空行
{$1="chr"$1;print}
' $wdir/0.raw/${key}.di | sort -k1V,1 -k2n,2 > $wdir/4.directional_index/${key}.bedGraph 


