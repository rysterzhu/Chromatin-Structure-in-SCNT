#!/bin/bash
res=100000
chr=chr2
s=0
e=""
g=""
ext=png

while [ -n "$1" ]; do
case "$1" in
    -h) help; shift 1;;
    -w) wdir=$2; shift 2;;
    -d) ddir=$2; shift 2;;
	-r) res=$2; shift 2;;
	-c) chr=$2; shift 2;;
	-s) s=$2; shift 2;;
	-e) e=$2; shift 2;;
	-g) g=$2; shift 2;;
	-ext) ext=$2; shift 2;;
    -*) echo "error: no such option $1. -h for help";exit 1;;
     *) break;
esac
done
mkdir -p $wdir
k1=("cc" "05h" "1h" "2h" "3h" "6h" "12h")
k2=("e2cell" "l2cell" "4cell" "8cell" "morula" "icm" "te")

if [ "$e" == "" ]; then 
ee=""; 
else 
ee="-e $e"
fi

python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir/05h_${res}_iced.matrix $ddir/1h_${res}_iced.matrix $ddir/2h_${res}_iced.matrix $ddir/3h_${res}_iced.matrix $ddir/6h_${res}_iced.matrix $ddir/12h_${res}_iced.matrix \
-tri 1 -bed $ConfigHP/${res}_mm10.bed -n ${k1[@]} -chr $chr -s $s $ee -o $wdir/1cell -ext $ext -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 &

python ~/codes/HiCPlotter.py -f $ddir/e2cell_${res}_iced.matrix $ddir/l2cell_${res}_iced.matrix $ddir/4cell_${res}_iced.matrix $ddir/8cell_${res}_iced.matrix $ddir/morula_${res}_iced.matrix $ddir/icm_${res}_iced.matrix $ddir/te_${res}_iced.matrix \
-tri 1 -bed $ConfigHP/${res}_mm10.bed -n ${k2[@]} -chr $chr -s $s $ee -o $wdir/2cell -ext $ext -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 &

wait

