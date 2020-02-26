#!/bin/bash
res=100000

while [ -n "$1" ]; do
case "$1" in
    -h) help; shift 1;;
    -w) wdir=$2; shift 2;;
    -d) ddir=$2; shift 2;;
	-r) res=$2; shift 2;;
	-e) ratio=$2; shift 2;;
    -*) echo "error: no such option $1. -h for help";exit 1;;
     *) break;
esac
done


#ln -s $wdir/hic_results/matrix/*/iced/$res/*_iced.matrix $wdir/3.matrix


# for i in $wdir/3.matrix/*${res}_iced.matrix; do 
# awk -v OFS="\t" -v key=${i##*/} '{depth+=$3} END{print key,depth}' $i >> $wdir/4.norm_to_depth/depth.tab &
# echo "calc depth: "$i
# done
# wait 
# if [ "$ratio" == "" ]; then 
# ratio=`grep _${res}_ $wdir/depth.tab|awk '{a+=$2;b+=1} END{print a/b}'` 
# fi

mkdir -p $wdir

for i in $ddir/*${res}_iced.matrix; do key=${i##*/};key=${key%%_*}
awk -v OFS="\t" -v key=${i##*/} -v ratio=$ratio 'NR==FNR{if(key==$1)depth=$2} NR!=FNR{print $1,$2,$3*ratio/depth}' $wdir/depth.tab $i > $wdir/4.norm_to_depth/${i##*/} &
echo "normalize: "$i
done
wait

# for i in $wdir/4.norm_to_depth/*${res}_iced.matrix; do o=${i##*/}
# awk -v res=$res '{$3=$3*($2-$1+1)*4e4/res;print $0}' $i > $wdir/5.norm_to_ps/${o} &
# done

wait
