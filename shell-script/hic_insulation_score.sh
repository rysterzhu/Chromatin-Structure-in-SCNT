#!/bin/bash
res=40000

while [ -n "$1" ]; do
case "$1" in
    -h) help; shift 1;;
    -w) wdir=$2; shift 2;;
    -d) ddir=$2; shift 2;;
	-r) res=$2; shift 2;;
    -*) echo "error: no such option $1. -h for help";exit 1;;
     *) break;
esac
done
# --is       [500000]   optional, insulation square size, size of the insulation square in bp
	# --ss       [0]        optional, insulation smooth size, size of insulation vector smoothing vector in bp
	# --ids      [250000]   optional, insulation delta span (window size of insulationd delta - blue line)
	# --im       [mean]     optional, insulation mode (how to aggregrate signal within insulation square), mean,sum,median
	# --nt       [0.1]      optional, noise threshold, minimum depth of valley
	# --bmoe     [3]        optional, boundary margin of error (specified in number of BINS), added to each side of the boundary
	# --yb       [auto]     optional, -yBound - +yBound of insulation plot
	# --bg       []         FLAG, use transparent background of insulation plot
	# --minDist  []         minimum allowed interaction distance, exclude < N distance (in BP)
	# --maxDist  []         maximum allowed interaction distance, exclude > N distance (in BP)
	# --ez       []         FLAG, ignore 0s in all calculations
	
		# --mbs      []         min boundary strength (0-inf)
	# --mts      []         min tad strength (0-inf)

#xw res:40k is:1M ids:200k --mbs:0.25
#lj res:40k is:480k
	
#keys=("cc" "05h" "1h" "2h" "3h" "6h" "12h" "e2cell" "l2cell" "4cell" "8cell" "morula" "icm" "te")
#wdir=~/workspace/8.NT-HiC/f.IS_ALL/*
mkdir -p "$wdir/boundaries"         #boundaries split by chr sample
mkdir -p "$wdir/insulation"         #insulation score  split by chr sample
mkdir -p "$wdir/logs"
mkdir -p "$wdir/matrix"             #matrix convert from HiC-Pro
mkdir -p "$wdir/5.plotProfile"       #dir for 5.plotProfile
mkdir -p "$wdir/4.bw"  
mkdir -p "$wdir/picture"            #insulation plot 
mkdir -p "$wdir/1.cat_boundary"   #filt 0.25 boundary strength
mkdir -p "$wdir/2.cat_insulation"      #cat insulation to sample
mkdir -p "$wdir/tad"   
mkdir -p "$wdir/3.cat_tad"   

cworld=/usr/local/software/1.HiC-software/1.cworld
cworld=~/1.HiC-software/1.cworld

cd $wdir
for j in $(seq 1 19) X; do chr=chr$j
for i in $ddir/*_${res}_iced.matrix; do o=${i##*/};o=${o/_${res}_iced.matrix}
#python ~/codes/convert_3col_to_matrix_for_cworld.py -i $i -I $ConfigHP/${res}_mm10.bed -c $chr -o $o.$chr.mat
#perl -I ~/1.HiC-software/1.cworld/lib/ ~/1.HiC-software/1.cworld/perl/matrix2insulation.pl --is 1000000 --ids 200000 --nt 0 --nt 0.1 --bmoe 3 -i $o.$chr.mat
perl -I $cworld/lib/ $cworld/perl/matrix2insulation.pl --is 1000000 --ids 200000 --nt 0.1 --bmoe 3 -i $o.$chr.mat #0.5?
#perl -I ~/1.HiC-software/1.cworld/lib/ ~/1.HiC-software/1.cworld/perl/matrix2insulation.pl --is 2000000 --ids 400000 --nt 0.1 --bmoe 1 -i $o.$chr.mat #test14
#perl -I ~/1.HiC-software/1.cworld/lib/ ~/1.HiC-software/1.cworld/perl/matrix2insulation.pl --is 480000 --ids 120000 --nt 0.1 --bmoe 3 --im iqrMean -i $o.$chr.mat
perl -I $cworld/lib/ $cworld/perl/insulation2tads.pl -i $o.$chr--*insulation -b $o.$chr--*boundaries -o $o.$chr --mbs 0.25 
echo "calc IS: "$o.$chr
done &
done 

wait


mv *log logs/
mv *.boundaries* boundaries/
mv *mat matrix/
mv *png picture/
mv *pdf picture/
mv *.insulation* insulation/
mv *.tad* tad/

for i in $ddir/*_${res}_iced.matrix; do 
o=${i##*/};o=${o/_${res}_iced.matrix}
cat $wdir/boundaries/$o.chr*boundaries.bed | awk -v OFS="\t" '$1~/^chr/{print}' | sort -k1V,1 -k2n,2 > $wdir/1.cat_boundary/$o.bed
#awk -v OFS="\t" '$5>0.25{print}' $wdir/1.cat_boundary/$o.bed > $wdir/1.cat_boundary/filted/$o.filted.bed
cat $wdir/insulation/$o.chr*bedGraph | awk -v OFS="\t" '$1~/^chr/{if($4=="NA")$4=0;print}' | sort -k1V,1 -k2n,2 > $wdir/2.cat_insulation/$o.bedGraph
cat $wdir/tad/$o.chr*bedGraph | awk -v OFS="\t" '$1~/^chr/{if($4=="NA")$4=0;print}' | sort -k1V,1 -k2n,2 > $wdir/3.cat_tad/$o.bed
done

for i in $wdir/1.cat_boundary/*bed; do k=$(basename $i .bed)
awk -v i=$k '{print i,$5}' $i >> $wdir/1.cat_boundary/all.ins
done

# for i in ${keys[@]}; do 
# awk -v i=${i} '{print i,$5}' $wdir/1.cat_boundary/filted/$i.filted.bed >> $wdir/1.cat_boundary/filted/all.ins
# done

for i in $wdir/2.cat_insulation/*bedGraph; do k=$(basename $i .bedGraph)
bedSort $i $wdir/4.bw/$k.bedGraph
bedGraphToBigWig $wdir/4.bw/$k.bedGraph ~/ann/mm10.chromSizes $wdir/4.bw/$k.bw
done

