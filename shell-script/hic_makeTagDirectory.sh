#!/bin/bash
threads=4
while [ -n "$1" ]; do
case "$1" in
    -h) help; shift 1;;
    -w) wdir=$2; shift 2;;
    -d) ddir=$2; shift 2;;
	-t) threads=$2; shift 2;;
    -*) echo "error: no such option $1. -h for help";exit 1;;
     *) break;
esac
done

trap "exec 1000>&-;exec 1000<&-;exit 0" 2 
mkfifo temp.fifo 
exec 1000<>temp.fifo 
for((i=0;i<$threads;i++));do echo >&1000; done

# ddir=~/workspace/8.NT-HiC/d.HiC-Pro_analysis/a.merge_samples_except/2.allValidPairs
# wdir=~/workspace/8.NT-HiC/d.HiC-Pro_analysis/d.homer_samples_expcept

mkdir -p $wdir/logs $wdir/1.tags
for i in $ddir/*_allValidPairs; do 
read -u 1000 
{
	o=${i##*/}; o=${o/_allValidPairs/.tag}
	nohup makeTagDirectory $wdir/1.tags/$o $i -format HiCsummary > $wdir/logs/${o/.tag/.log} 2>&1
	echo >&1000 
	echo $i" done"
}&   
done &

wait
echo "All threads done."
rm -rf temp.fifo