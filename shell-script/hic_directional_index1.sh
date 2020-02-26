#!/bin/bash
res=40000
length=2000000
while [ -n "$1" ]; do
case "$1" in
    -h) help; shift 1;;
    -w) wdir=$2; shift 2;;
    -m) matrix=$2; shift 2;;
	-r) res=$2; shift 2;;
	-l) length=$2; shift 2;;
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
##3.
for i in *$key.mat; do
perl ~/1.HiC-software/8.DI/1DI_from_matrix.pl $i ${res} $length $ConfigHP/chrom_mm10.sizes > ${i/mat/di}
done
##4.cat
cat *$key.di | awk -v OFS="\t" '$1~/[0-9]+/{$0=$0;print} $1~"X"{$1=20;print} $1~"Y"{$1=21;print} ' | sort -k1n,1 -k2n,2 > $wdir/0.raw/${key}.di 
echo "calc di done: "$i

