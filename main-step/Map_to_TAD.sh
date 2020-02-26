tdir=~/workspace/8.NT-HiC/g.DI_ALL/9.except_norm_mean_raw_20180713/2.unoverlap_tads
ddir=~/workspace/8.NT-HiC/5.maps/3.sample/2.iced
wdir=~/workspace/8.NT-HiC/g.DI_ALL/c.map_to_tad

for i in ${keys[@]}; do 
python ~/codes/map_to_TAD.py -i $ddir/${i}_40000_iced.matrix -I $tdir/icm.tad -b $ConfigHP/40000_mm10.bed -o $wdir/${i}_icm_40000.mat
done

python ~/codes/HiCPlotter.py -f $wdir/${keys[0]}_icm_40000.mat $wdir/${keys[1]}_icm_40000.mat $wdir/${keys[2]}_icm_40000.mat $wdir/${keys[3]}_icm_40000.mat $wdir/${keys[4]}_icm_40000.mat $wdir/${keys[5]}_icm_40000.mat $wdir/${keys[6]}_icm_40000.mat \
-fh 0 -n cc 05h 1h 2h 3h 6h 12h -o $wdir/all1 -ext png -dpi 360 -chr chr0 &
python ~/codes/HiCPlotter.py -f $wdir/${keys[7]}_icm_40000.mat $wdir/${keys[8]}_icm_40000.mat $wdir/${keys[9]}_icm_40000.mat $wdir/${keys[10]}_icm_40000.mat $wdir/${keys[11]}_icm_40000.mat $wdir/${keys[12]}_icm_40000.mat $wdir/${keys[13]}_icm_40000.mat \
-fh 0 -n e2cell l2cell 4cell 8cell morula icm te -o $wdir/all2 -ext png -dpi 360 -chr chr0 &

python ~/codes/HiCPlotter.py -f $wdir/${keys[2]}_icm_40000.mat $wdir/${keys[3]}_icm_40000.mat $wdir/${keys[4]}_icm_40000.mat $wdir/${keys[5]}_icm_40000.mat $wdir/${keys[6]}_icm_40000.mat $wdir/${keys[7]}_icm_40000.mat $wdir/${keys[8]}_icm_40000.mat \
-fh 0 -n 1h 2h 3h 6h 12h e2cell l2cell -o $wdir/1cell_test2 -ext png -dpi 360 -chr chr0 -hmc 1 -mm 12

wdir=~/workspace/8.NT-HiC/g.DI_ALL/c.map_to_tad/1.diff_between_stages

python ~/codes/HiCPlotter.py -f $wdir/${keys[1]}_icm_40000.mat $wdir/${keys[2]}_icm_40000.mat $wdir/${keys[3]}_icm_40000.mat $wdir/${keys[4]}_icm_40000.mat $wdir/${keys[5]}_icm_40000.mat $wdir/${keys[6]}_icm_40000.mat \
-fh 0 -n 05h 1h 2h 3h 6h 12h -o $wdir/all1 -ext png -dpi 360 -chr chr0 -hmc 5 &
python ~/codes/HiCPlotter.py -f $wdir/${keys[7]}_icm_40000.mat $wdir/${keys[8]}_icm_40000.mat $wdir/${keys[9]}_icm_40000.mat $wdir/${keys[10]}_icm_40000.mat $wdir/${keys[11]}_icm_40000.mat $wdir/${keys[12]}_icm_40000.mat $wdir/${keys[13]}_icm_40000.mat \
-fh 0 -n e2cell l2cell 4cell 8cell morula icm te -o $wdir/all2 -ext png -dpi 360 -chr chr0 &

#####################HiC Plotter 不好用，全部用R画


##################xw
tdir=~/workspace/8.NT-HiC/3.public_data/f.XW_analysis/5.DI/2.unoverlap_tads
ddir=~/workspace/8.NT-HiC/5.maps/3.sample/2.iced
wdir=~/workspace/8.NT-HiC/3.public_data/f.XW_analysis/8.map_to_tad

for i in ${keys[@]}; do 
python ~/codes/map_to_TAD.py -i $ddir/${i}_40000_iced.matrix -I $tdir/ICM.tad -b $ConfigHP/40000_mm10.bed -o $wdir/${i}_icm_40000.mat
done


###########################2018年8月18日 
tdir=~/workspace/8.NT-HiC/g.DI_ALL/11.except_norm_1e8_20180813/2.unoverlap_tads
ddir=~/workspace/8.NT-HiC/5.maps/4.except/2.iced
wdir=~/workspace/8.NT-HiC/g.DI_ALL/c.map_to_tad/2.NT_ICM_20180818

for i in ${keys[@]}; do 
python ~/codes/map_to_TAD.py -i $ddir/${i}_40000_iced.matrix -I $tdir/icm.tad -b $ConfigHP/40000_mm10.bed -o $wdir/${i}_icm_40000.mat &
done

###########################2018年8月18日 deDoc
ddir=~/workspace/8.NT-HiC/5.maps/4.except/2.iced
wdir=~/workspace/8.NT-HiC/g.DI_ALL/c.map_to_tad/3.NT_deDocE_20180818
mkdir -p $wdir
for i in ${keys[@]}; do 
python ~/codes/map_to_TAD.py -i $ddir/${i}_40000_iced.matrix -I ~/workspace/8.NT-HiC/8.softeware_test/g.deDoc/5.NT/allM_sorted.tad -b $ConfigHP/40000_mm10.bed -o $wdir/${i}_icm_40000.mat &
done