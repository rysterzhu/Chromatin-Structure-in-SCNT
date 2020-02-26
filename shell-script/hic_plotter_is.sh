#!/bin/bash
res=100000
chr=chr2
s=0
e=""
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
	-tad) tad=$2; shift 2;;
	-di) di=$2; shift 2;;
	-is) is=$2; shift 2;;
	-ext) ext=$2; shift 2;;
    -*) echo "error: no such option $1. -h for help";exit 1;;
     *) break;
esac
done


k1=("cc" "05h" "1h" "2h" "3h" "6h")
k2=("12h" "e2cell" "l2cell" "4cell" "8cell" "morula")

k1=("cc" "05h" "1h" "2h" "3h" "6h" "12h")
k2=("e2cell" "l2cell" "4cell" "8cell" "morula" "icm" "te")

if [ "$e" == "" ]; then 
ee=""; 
else 
ee="-e $e"
fi
mkdir -p $wdir

python ~/codes/HiCPlotter.py -n ${k1[@]} -tri 1 -bed $ConfigHP/${res}_mm10.bed -r ${res} \
-f $ddir/cc_${res}_iced.matrix $ddir/05h_${res}_iced.matrix $ddir/1h_${res}_iced.matrix $ddir/2h_${res}_iced.matrix $ddir/3h_${res}_iced.matrix $ddir/6h_${res}_iced.matrix $ddir/12h_${res}_iced.matrix \
-hist $di/cc.bedGraph,$is/cc.bedGraph $di/05h.bedGraph,$is/05h.bedGraph $di/1h.bedGraph,$is/1h.bedGraph $di/2h.bedGraph,$is/2h.bedGraph $di/3h.bedGraph,$is/3h.bedGraph $di/6h.bedGraph,$is/6h.bedGraph $di/12h.bedGraph,$is/12h.bedGraph \
-hl Directional_Index,Insulation_Score Directional_Index,Insulation_Score Directional_Index,Insulation_Score Directional_Index,Insulation_Score Directional_Index,Insulation_Score Directional_Index,Insulation_Score Directional_Index,Insulation_Score \
-chr $chr -s $s $ee -o $wdir/1cell -ext $ext -dpi 720 -spi 1 -hmc 1 -ptr 1 &

python ~/codes/HiCPlotter.py -n ${k2[@]} -tri 1 -bed $ConfigHP/${res}_mm10.bed -r ${res} \
-f $ddir/e2cell_${res}_iced.matrix $ddir/l2cell_${res}_iced.matrix $ddir/4cell_${res}_iced.matrix $ddir/8cell_${res}_iced.matrix $ddir/morula_${res}_iced.matrix $ddir/icm_${res}_iced.matrix $ddir/te_${res}_iced.matrix \
-hist $di/e2cell.bedGraph,$is/e2cell.bedGraph $di/l2cell.bedGraph,$is/l2cell.bedGraph $di/4cell.bedGraph,$is/4cell.bedGraph $di/8cell.bedGraph,$is/8cell.bedGraph $di/morula.bedGraph,$is/morula.bedGraph $di/icm.bedGraph,$is/icm.bedGraph $di/te.bedGraph,$is/te.bedGraph \
-hl Directional_Index,Insulation_Score Directional_Index,Insulation_Score Directional_Index,Insulation_Score Directional_Index,Insulation_Score Directional_Index,Insulation_Score Directional_Index,Insulation_Score Directional_Index,Insulation_Score \
-chr $chr -s $s $ee -o $wdir/2cell -ext $ext -dpi 720 -spi 1 -hmc 1 -ptr 1 &

wait

