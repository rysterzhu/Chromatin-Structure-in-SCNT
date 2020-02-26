#!/bin/bash
#res=40000
while [ -n "$1" ]; do
case "$1" in
    -h) help; shift 1;;
    -w) wdir=$2; shift 2;;
#    -d) ddir=$2; shift 2;;
#	-r) res=$2; shift 2;;
    -*) echo "error: no such option $1. -h for help";exit 1;;
     *) break;
esac
done

mkdir $wdir/stats
cd $wdir
cp $wdir/*/hic_results/data/*/*stat $wdir/stats/
cp $wdir/*/bowtie_results/bwt2/*/*stat $wdir/stats/
cd stats
echo -e "Sample\t\
Total_pairs_processed\t\
Reported_pairs\t%map_pairs\t\
valid_interaction\t%valid_interaction\t\
valid_interaction_rmdup\t%Duplicates\t\
%final_valid\t\
trans_interaction\t%trans\t\
cis_interaction\t%cis\t\
cis_shortRange\tcis_longRange" > $wdir/stats/mergestats

for i in $wdir/stats/*mpairstat; do 
k=${i##*/}
k=${k%%.*}
j=$wdir/stats/${k}_allValidPairs.mergestat
awk -v k=$k 'BEGIN{}
$1=="Total_pairs_processed"{a=($2/1e6)}
$1=="Reported_pairs"{b=($2/1e6)}
$1=="valid_interaction"{c=($2/1e6)}
$1=="valid_interaction_rmdup"{d=($2/1e6)}
$1=="trans_interaction"{e=($2/1e6)}
$1=="cis_interaction"{f=($2/1e6)}
$1=="cis_shortRange"{g=($2/1e6)}
$1=="cis_longRange"{h=($2/1e6)}
END{printf("%s\t%.2f\t%.2f\t%.2f%\t%.2f\t%.2f%\t%.2f\t%.2f%\t%.2f%\t%.2f\t%.2f%\t%.2f\t%.2f%\t%.2f\t%.2f\n",k,a,b,100*b/a,c,100*c/b,d,(c-d)*100/c,100*d/a, e,100*e/d,f,100*f/d,g,h)}
' $i $j >> $wdir/stats/mergestats
done

