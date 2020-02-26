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

mkdir -p $wdir/0.raw $wdir/3.plot_tads $wdir/1.tads $wdir/4.directional_index $wdir/2.unoverlap_tads $wdir/logs $wdir/6.JuiceBox_tad
key=${matrix##*/}; key=${key/_${res}_iced.matrix}
##6.
echo -en $key\\t
perl ~/1.HiC-software/8.DI/3file_ends_cleaner.pl $wdir/0.raw/${key}.hmm $wdir/0.raw/${key}.di | perl ~/1.HiC-software/8.DI/4converter_7col.pl > $wdir/0.raw/${key}.7col
echo -en 7col\\t
##6.5 split 7col to chrs
rm -f $wdir/0.raw/$key/*7col
awk -v dir=$wdir/0.raw/$key '{print >> dir"/"$1".7col"}' $wdir/0.raw/${key}.7col
#7.
for i in $wdir/0.raw/$key/*7col; do chr=$(basename $i .7col)
echo $chr
perl ~/1.HiC-software/8.DI/5hmm_probablity_correcter.pl $i 2 0.99 $res | \
perl ~/1.HiC-software/8.DI/6hmm-state_caller.pl $ConfigHP/chrom_mm10.sizes $chr | \
perl ~/1.HiC-software/8.DI/7hmm-state_domains.pl > ${i/7col/raw.tad}
done
echo -en raw.tad\\t

#8.others
# 去除一些错误的TAD
cat $wdir/0.raw/$key/*raw.tad | \
awk -v OFS="\t" '{chrs[$1]}
$2>$3{t=$2;$2=$3;$3=t;} #将TAD的方向换成正的
$1=="chr20"{$1="chrX"}  #将TAD的chr20换成X
$1=="chr21"||$1==""||$2==""||$3==""{next}   #不要chrY,空行
{$0=$0;print}
' | sort -k1V,1 -k2n,2 > $wdir/1.tads/${key}.tad

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

# 转Juicebox格式
awk -v OFS="\t" 'BEGIN{OFS="\t";print "chr1\tx1\tx2\tchr2\ty1\ty2\tcolor\tcomment";n=1} 
{print $1,$2,$3,$1,$2,$3,"0,0,0","Tad"n;n+=1}' $wdir/2.unoverlap_tads/${key}.tad > $wdir/6.JuiceBox_tad/${key}.tad
echo done.