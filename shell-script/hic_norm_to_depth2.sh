#!/bin/bash
res=100000

while [ -n "$1" ]; do
case "$1" in
    -h) help; shift 1;;
    -w) wdir=$2; shift 2;;
    -d) ddir=$2; shift 2;;
	-r) res=$2; shift 2;;
	-f) f=$2; shift 2;;
    -*) echo "error: no such option $1. -h for help";exit 1;;
     *) break;
esac
done

mkdir -p $wdir
for i in $ddir/*$res*matrix; do 
awk -v OFS="\t" -v key=${i##*/} '$1!=$2{depth+=$3} END{print key,depth}' $i >> $wdir/depth.tab &
echo "calc depth: "$i
done
wait 

for i in $ddir/*$res*matrix; do 
awk -v OFS="\t" -v key=${i##*/} -v f=$f 'NR==FNR{if(key==$1)depth=$2} NR!=FNR{print $1,$2,$3*f/depth}' $wdir/depth.tab $i > $wdir/${i##*/} &
echo "normalize: "$i
done

