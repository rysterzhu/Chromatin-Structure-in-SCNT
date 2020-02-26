#!/bin/bash --login
res=100000
bs=50
bc=0
while [ -n "$1" ]; do
case "$1" in
    -h) help; shift 1;;
    -w) wdir=$2; shift 2;;
    -d) ddir=$2; shift 2;;
	-r) res=$2; shift 2;;
	-b) bs=$2; shift 2;;
	-c) bc=$2; shift 2;;
    -*) echo "error: no such option $1. -h for help";exit 1;;
     *) break;
esac
done

 
mkdir -p $wdir
#取100k分辨率的matrix，分出interchrom和intrachrom
for i in $ddir/*${res}[._]*matrix; do o=${i##*/};
awk -v OFS="\t" -v intra=$wdir/${o/matrix/intra} -v inter=$wdir/${o/matrix/inter} 'NR==FNR{if($4>a[$1]){a[$1]=$4}} 
NR>FNR{for(i in a){if($1<=a[i]&&$2>a[i]){print $0 > inter;next}}; print $0 > intra}' $ConfigHP/${res}_mm10.bed $i &
echo "split chromosome: "$i
done
wait

for i in $ddir/*${res}[._]*matrix; do o=${i##*/};
awk -v OFS="\t" -v bs=$bs -v bc=$bc '$2-$1<=bs&&$2-$1>=bc{print $1"-"$2,$3}' $wdir/${o/matrix/intra} > $wdir/${o/matrix/txt} &
echo "screen "$bs" bins: "$wdir/${o/matrix/intra}
done
wait

#
for i in $ddir/*${res}[._]*matrix; do o=${i##*/};
awk -v OFS="\t" 'NR==FNR{a+=$3} NR!=FNR{$2=$2*1e8/a;print }' $i $wdir/${o/matrix/txt} > $wdir/${o/matrix/norm} &
echo "normalize: "$wdir/${o/matrix/txt}
done
wait

#rm *inter *intra *txt
