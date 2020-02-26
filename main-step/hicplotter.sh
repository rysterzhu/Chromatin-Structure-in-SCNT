python ~/codes/HiCPlotter.py  -f ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/xw_rep1_40000_iced.matrix -tri 1 -bed ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/xw_rep1_40000_abs.bed -n xw_rep1 -chr chr19 -o xw_rep1 -ext png -dpi 200 -r 40000 
21944  python ~/codes/HiCPlotter.py  -f ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/xw_rep1_40000_iced.matrix -tri 1 -bed ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/xw_rep1_40000_abs.bed -n xw_rep1 -chr chr19 -o xw_rep1 -ext png -dpi 200 -r 100000 -spi 1 -hmc 1 
21946  python ~/codes/HiCPlotter.py  -f ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/xw_rep1_100000_iced.matrix -tri 1 -bed ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/xw_rep1_100000_abs.bed -n xw_rep1 -chr chr19 -s 0 -e 61 -o xw_rep1 -ext png -dpi 360 -r 100000 -spi 1 -hmc 1 
21948  python ~/codes/HiCPlotter.py  -f ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/xw_rep1_40000_iced.matrix -tri 1 -bed ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/xw_rep1_100000_abs.bed -n xw_rep1 -chr chr19 -s 0 -e 61 -o xw_rep1 -ext png -dpi 360 -r 40000 -spi 1 -hmc 1 
21949  python ~/codes/HiCPlotter.py  -f ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/xw_rep1_40000_iced.matrix -tri 1 -bed ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/xw_rep1_40000_abs.bed -n xw_rep1 -chr chr19 -s 0 -e 61 -o xw_rep1 -ext png -dpi 360 -r 100000 -spi 1 -hmc 1 
21950  python ~/codes/HiCPlotter.py  -f ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/xw_rep1_40000_iced.matrix -tri 1 -bed ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/xw_rep1_40000_abs.bed -n xw_rep1 -chr chr19 -s 0 -e 61e6 -o xw_rep1 -ext png -dpi 360 -r 100000 -spi 1 -hmc 1 
21951  python ~/codes/HiCPlotter.py  -f ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/xw_rep1_40000_iced.matrix -tri 1 -bed ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/xw_rep1_40000_abs.bed -n xw_rep1 -chr chr19 -s 0 -e 61000000 -o xw_rep1 -ext png -dpi 360 -r 100000 -spi 1 -hmc 1 
21952  python ~/codes/HiCPlotter.py  -f ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/xw_rep1_40000_iced.matrix -tri 1 -bed ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/xw_rep1_40000_abs.bed -n xw_rep1 -chr chr19 -s 35000000 -e 45000000 -o xw_rep1 -ext png -dpi 360 -r 40000 -spi 1 -hmc 1 
21954  python ~/codes/HiCPlotter.py  -f ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/xw_rep1_40000_iced.matrix -tri 1 -bed ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/xw_rep1_40000_abs.bed -n xw_rep1 -chr chr19 -s 875 -e 1125 -o xw_rep1 -ext png -dpi 360 -r 40000 -spi 1 -hmc 1 
21956  python ~/codes/HiCPlotter.py  -f ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/xw_rep1_${res}_iced.matrix -tri 1 -bed ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/xw_rep1_${res}_abs.bed -n xw_rep1 -chr chr19 -o xw_rep1 -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 
21957  python ~/codes/HiCPlotter.py -h
21960  python ~/codes/HiCPlotter.py  -f $ddir/xw_rep1_${res}_iced.matrix $ddir/xw_rep2_${res}_iced.matrix -tri 1 -bed $ddir/xw_rep1_${res}_abs.bed $ddir/xw_rep2_${res}_abs.bed -n xw -chr chr19 -o xw -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 
21963  python ~/codes/HiCPlotter.py  -f $ddir/xw_rep1_${res}_iced.matrix $ddir/xw_rep2_${res}_iced.matrix -tri 1 -bed $ddir/xw_rep1_${res}_abs.bed -n xw -chr chr19 -o xw -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 
21966  python ~/codes/HiCPlotter.py  -f $ddir/xw_rep1_${res}_iced.matrix $ddir/xw_rep2_${res}_iced.matrix -tri 1 -bed $ddir/xw_rep1_${res}_abs.bed -n xw_rep1 xw_rep2 -chr chr19 -o xw -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 
21975  python ~/codes/HiCPlotter.py  -f $ddir/cc1_${res}_iced.matrix $ddir/cc2_${res}_iced.matrix -tri 1 -bed $ddir/cc1_${res}_abs.bed -n cc_rep1 cc_rep2 -chr chr19 -o cc -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 
21978  python ~/codes/HiCPlotter.py  -f $ddir/morula1_${res}_iced.matrix $ddir/morula2_${res}_iced.matrix -tri 1 -bed $ddir/morula1_${res}_abs.bed -n morula_rep1 morula_rep2 -chr chr19 -o morula -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 
21993  python ~/codes/HiCPlotter.py  -f $ddir/xw_rep1.chr19.simpleNorm.40kby40k.txt $ddir/xw_rep1.chr19.simpleNorm.40kby40k.txt -n xw_rep1 xw_rep2 -chr chr19 -o xw -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 
22059  for i in *40k.txt; do cut -f 1,2 --complement > ~/workspace/8.NT-HiC/4.HiCPlotter/1.simpleNorm/$i & done
22060  less ~/workspace/8.NT-HiC/4.HiCPlotter/1.simpleNorm/$i
22067  for i in *40k.txt; do cut -f 1,2 --complement $i > ~/workspace/8.NT-HiC/4.HiCPlotter/1.simpleNorm/$i & done
22072  ddir=~/workspace/8.NT-HiC/4.HiCPlotter/1.simpleNorm
22073  python ~/codes/HiCPlotter.py  -f $ddir/xw_rep1.chr19.simpleNorm.40kby40k.txt $ddir/xw_rep1.chr19.simpleNorm.40kby40k.txt -n xw_rep1 xw_rep2 -chr chr19 -o xw -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 
22076  ddir=~/workspace/8.NT-HiC/4.HiCPlotter/1.simpleNorm
22077  python ~/codes/HiCPlotter.py  -f $ddir/xw_rep1.chr19.simpleNorm.40kby40k.txt $ddir/xw_rep1.chr19.simpleNorm.40kby40k.txt -n xw_rep1 xw_rep2 -chr chr19 -fh 1 -o xw -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 
22079  ddir=~/workspace/8.NT-HiC/4.HiCPlotter/1.simpleNorm
22080  python ~/codes/HiCPlotter.py  -f $ddir/cc.1.chr19.simpleNorm.40kby40k.txt $ddir/cc.2.chr19.simpleNorm.40kby40k.txt -n cc_rep1 cc_rep2 -chr chr19 -fh 1 -o cc -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 
22082  for i in *4M.txt; do cut -f 1,2 --complement $i > ~/workspace/8.NT-HiC/4.HiCPlotter/1.simpleNorm/$i & done
22087  ddir=~/workspace/8.NT-HiC/4.HiCPlotter/1.simpleNorm
22088  python ~/codes/HiCPlotter.py  -f $ddir/cc.1.all.simpleNorm.4Mby4M.txt $ddir/cc.2.all.simpleNorm.4Mby4M.txt -n cc_rep1 cc_rep2 -chr chr19 -fh 1 -o cc_4M -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 
22090  ddir=~/workspace/8.NT-HiC/4.HiCPlotter/1.simpleNorm
22091  python ~/codes/HiCPlotter.py  -f $ddir/cc.1.chr19.simpleNorm.40kby40k.txt $ddir/cc.2.chr19.simpleNorm.40kby40k.txt -n cc_rep1 cc_rep2 -chr chr19 -fh 1 -o cc -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 -mm 5
22193  cd 4.HiCPlotter/
22196  python ~/codes/HiCPlotter.py  -f ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/cc1_${res}_iced.matrix -tri 1 -bed ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/cc1_${res}_abs.bed -n cc_rep1 cc_rep2 -chr chr19 -s 875 -e 1125 -o cc -ext png -dpi 360 -r ${res} \
22198  python ~/codes/HiCPlotter.py  -f ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/cc1_${res}_iced.matrix ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/cc2_${res}_iced.matrix -tri 1 -bed ~/workspace/8.NT-HiC/1.hicup/5.hic-pro/matrix/cc1_${res}_abs.bed -n cc_rep1 cc_rep2 -chr chr19 -s 875 -e 1125 -o cc -ext png -dpi 360 -r ${res} -spi 1 -hmc 1

##########2017年12月6日 es大测数据
ddir=~/workspace/8.NT-HiC/6.hic-pro_20171123/es/1.results
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/2.es_20171206
rdir=~/workspace/8.NT-HiC/6.hic-pro_20171123/1.regions
res=100000
key=es500
python ~/codes/HiCPlotter.py -f $ddir/es100_${res}_iced.matrix $ddir/es200_${res}_iced.matrix $ddir/es500_1_${res}_iced.matrix $ddir/es500_2_${res}_iced.matrix \
-tri 1 -bed $rdir/${res}_mm10.bed -n es100 es200 es500_1 es500_2 -chr chr19 -o $wdir/es_all -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 &
python ~/codes/HiCPlotter.py -f $ddir/es100_${res}_iced.matrix $ddir/es200_${res}_iced.matrix $ddir/es500_1_${res}_iced.matrix $ddir/es500_2_${res}_iced.matrix \
-tri 1 -bed $rdir/${res}_mm10.bed -n es100 es200 es500_1 es500_2 -chr chr19 -o $wdir/es_ptr -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 &  ##good
python ~/codes/HiCPlotter.py -f $ddir/es100_${res}_iced.matrix $ddir/es200_${res}_iced.matrix $ddir/es500_1_${res}_iced.matrix $ddir/es500_2_${res}_iced.matrix \
-tri 1 -bed $rdir/${res}_mm10.bed -n es100 es200 es500_1 es500_2 -chr chr19 -o $wdir/es_ptd -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 -ptd 1 &  ##error

##include xw rb ES
pdir=~/workspace/8.NT-HiC/6.hic-pro_20171123/xw/1.results
rbdir=~/workspace/8.NT-HiC/6.hic-pro_20171123/rb/1.results
python ~/codes/HiCPlotter.py -f $ddir/es100_${res}_iced.matrix $ddir/es200_${res}_iced.matrix $ddir/es500_1_${res}_iced.matrix $ddir/es500_2_${res}_iced.matrix \
$pdir/xw_rep1_${res}_iced.matrix $pdir/xw_rep2_${res}_iced.matrix $rbdir/rb_rep1_${res}_iced.matrix $rbdir/rb_rep2_${res}_iced.matrix \
-tri 1 -bed $rdir/${res}_mm10.bed -n es100 es200 es500_1 es500_2 xw_rep1 xw_rep2 rb_rep1 rb_rep2 -chr chr19 -o $wdir/es_pub_ptr -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 &


###########2017年12月7日 test DI tad domain
pdir=~/workspace/8.NT-HiC/6.hic-pro_20171123/xw/1.results
rbdir=~/workspace/8.NT-HiC/6.hic-pro_20171123/rb/1.results
python ~/codes/HiCPlotter.py -f $ddir/es100_${res}_iced.matrix $ddir/es200_${res}_iced.matrix $ddir/es500_1_${res}_iced.matrix $ddir/es500_2_${res}_iced.matrix \
$pdir/xw_rep1_${res}_iced.matrix $pdir/xw_rep2_${res}_iced.matrix $pdir/mII_rep1_${res}_iced.matrix $pdir/mII_rep2_${res}_iced.matrix \
-tri 1 -bed $rdir/${res}_mm10.bed -n es100 es200 es500_1 es500_2 xw_rep1 xw_rep2 mII_rep1 mII_rep2 -chr chr1 -o $wdir/es_pub_pcd \
-ext png -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 -pcd 1 -pcdf ~/workspace/8.NT-HiC/8.softeware_test/1.DI/chr1.tad ~/workspace/8.NT-HiC/8.softeware_test/1.DI/chr1.tad ~/workspace/8.NT-HiC/8.softeware_test/1.DI/chr1.tad ~/workspace/8.NT-HiC/8.softeware_test/1.DI/chr1.tad ~/workspace/8.NT-HiC/8.softeware_test/1.DI/chr1.tad ~/workspace/8.NT-HiC/8.softeware_test/1.DI/chr1.tad ~/workspace/8.NT-HiC/8.softeware_test/1.DI/chr1.tad ~/workspace/8.NT-HiC/8.softeware_test/1.DI/chr1.tad &

python ~/codes/HiCPlotter.py -f $pdir/xw_rep1_${res}_iced.matrix \
-tri 1 -bed $rdir/${res}_mm10.bed -n xw_rep1 -chr chr19 -o $wdir/xw_rep1_pcd2 \
-ext png -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 -pcd 1 -pcdf ~/workspace/8.NT-HiC/8.softeware_test/1.DI/tad_files/chr19.filt.tad &




###############2017年12月31日 
ddir=~/workspace/8.NT-HiC/6.hic-pro_20171123/morula_8cell_20171229/1.results
ddir2=~/workspace/8.NT-HiC/6.hic-pro_20171123/merge_8cell_20171231/1.results
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/3.8cell
rdir=~/workspace/8.NT-HiC/6.hic-pro_20171123/1.regions
res=100000
key=8cell
python ~/codes/HiCPlotter.py -f $ddir2/8cell_${res}_iced.matrix $ddir/8cell_rep1_${res}_iced.matrix $ddir/8cell_rep2_${res}_iced.matrix $ddir/8cell_rep3_${res}_iced.matrix \
-tri 1 -bed $rdir/${res}_mm10.bed -n merge 8cell_rep1 8cell_rep2 8cell_rep3 -chr chr19 -o $wdir/8cell -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 &

ddir=~/workspace/8.NT-HiC/6.hic-pro_20171123/merge_8cell_20171231/1.results
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/3.8cell
rdir=~/workspace/8.NT-HiC/6.hic-pro_20171123/1.regions
res=100000
key=8cell
python ~/codes/HiCPlotter.py -f $ddir/8cell_rep1_${res}_iced.matrix $ddir/8cell_rep2_${res}_iced.matrix $ddir/8cell_rep3_${res}_iced.matrix \
-tri 1 -bed $rdir/${res}_mm10.bed -n 8cell_rep1 8cell_rep2 8cell_rep3 -chr chr19 -o $wdir/8cell -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 &





ddir=~/workspace/8.NT-HiC/6.hic-pro_20171123/merge_morula1_20171231/1.results
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/4.morula
rdir=~/workspace/8.NT-HiC/6.hic-pro_20171123/1.regions
res=100000
key=morula1
python ~/codes/HiCPlotter.py -f $ddir/${key}_${res}_iced.matrix \
-tri 1 -bed $rdir/${res}_mm10.bed -n morula1 -chr chr19 -o $wdir/morula1_20171231 -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 &

ddir=~/workspace/8.NT-HiC/6.hic-pro_20171123/morula_8cell_20171229/1.results
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/4.morula
rdir=~/workspace/8.NT-HiC/6.hic-pro_20171123/1.regions
res=100000
key=Morula3
python ~/codes/HiCPlotter.py -f $ddir/${key}_${res}_iced.matrix \
-tri 1 -bed $rdir/${res}_mm10.bed -n morula1 -chr chr19 -o $wdir/morula3_20171231 -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 &



#############################################################
##all sample
#############################################################
ddir=~/workspace/8.NT-HiC/d.HiC-Pro_analysis/1.merge_samples/3.matrix
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/7.all_sample_100k_20180127
res=500000
chrs=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chrX")
keys=("cc" "05h" "1h" "2h" "6h" "e2cell" "l2cell" "4cell" "8cell" "morula")
mkdir -p $wdir
for i in ${chrs[@]}; do 
python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix \
$ddir/05h_${res}_iced.matrix \
$ddir/1h_${res}_iced.matrix \
$ddir/2h_${res}_iced.matrix \
$ddir/6h_${res}_iced.matrix \
$ddir/e2cell_${res}_iced.matrix \
$ddir/l2cell_${res}_iced.matrix \
$ddir/4cell_${res}_iced.matrix \
$ddir/8cell_${res}_iced.matrix \
$ddir/morula_${res}_iced.matrix \
-tri 1 -bed $ConfigHP/${res}_mm10.bed -n ${keys[@]} -chr $i -o $wdir/$i -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 &
done

res=4000000
python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix \
$ddir/05h_${res}_iced.matrix \
$ddir/1h_${res}_iced.matrix \
$ddir/2h_${res}_iced.matrix \
$ddir/6h_${res}_iced.matrix \
$ddir/e2cell_${res}_iced.matrix \
$ddir/4cell_${res}_iced.matrix \
$ddir/8cell_${res}_iced.matrix \
$ddir/morula_${res}_iced.matrix \
-tri 1 -bed $ConfigHP/${res}_mm10.bed -n ${keys[@]} -wg 1 -chr chr19 -o $wdir/wg -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 &
#not work: Whole genome can be plotted only as matrix - this feature will be improved in future releases

ddir=~/workspace/8.NT-HiC/d.HiC-Pro_analysis/a.merge_samples_except/3.matrix
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/9.except_6hrep2_20180323
res=100000
python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix \
$ddir/05h_${res}_iced.matrix \
$ddir/1h_${res}_iced.matrix \
-tri 1 -bed $ConfigHP/${res}_mm10.bed -n cc 05h 1h -chr chr2 -s 0 -e 250 -o $wdir/cc-1h -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 &

python ~/codes/HiCPlotter.py -f $ddir/2h_${res}_iced.matrix \
$ddir/6h_${res}_iced.matrix \
$ddir/e2cell_${res}_iced.matrix \
$ddir/l2cell_${res}_iced.matrix \
-tri 1 -bed $ConfigHP/${res}_mm10.bed -n 2h 6h e2cell l2cell -chr chr2 -s 0 -e 250 -o $wdir/2h-e2cell -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 &

#highlight test
tad=~/workspace/8.NT-HiC/d.HiC-Pro_analysis/8.IS_samples_except/cat_boundary/tad
python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix \
$ddir/05h_${res}_iced.matrix \
$ddir/1h_${res}_iced.matrix \
-tri 1 -bed $ConfigHP/${res}_mm10.bed -n cc 05h 1h -chr chr2 -s 0 -e 250 -o $wdir/cc-1h_tad -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 \
-pcd 1 -pcdf $tad/cc.all.tad.bed $tad/05h.all.tad.bed $tad/1h.all.tad.bed &

res=500000
python ~/codes/HiCPlotter.py -f $ddir/2h_${res}_iced.matrix \
$ddir/6h_${res}_iced.matrix \
$ddir/e2cell_${res}_iced.matrix \
$ddir/l2cell_${res}_iced.matrix \
-tri 1 -bed $ConfigHP/${res}_mm10.bed -n 2h 6h e2cell l2cell -chr chr2 -o $wdir/2h-e2cell_tad -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 \
-pcd 1 -pcdf $tad/2h.all.tad.bed $tad/6h.all.tad.bed $tad/e2cell.all.tad.bed $tad/l2cell.all.tad.bed &
#有些tad不存在报错


############insultion score
python ~/codes/HiCPlotter.py -tri 1 -bed $ConfigHP/${res}_mm10.bed -f $ddir/cc_${res}_iced.matrix \
-n cc -chr chr2 -r 40000 -o $wdir/cc_is -s 100 -e 700 \
-hist ~/workspace/8.NT-HiC/d.HiC-Pro_analysis/i.IS_norm_to_depth/insulation/cc.insulation.bedGraph -hl insulation 


#############################################################
##xw development samples
#############################################################
ddir=~/workspace/8.NT-HiC/3.public_data/4.HiC-Pro_xw/1.results
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/8.xw_development_100k_20180127
res=100000
chrs=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chrX")
mkdir $wdir
for i in ${chrs[@]}; do 
python ~/codes/HiCPlotter.py -f $ddir/sperm_${res}_iced.matrix \
$ddir/MII_${res}_iced.matrix \
$ddir/PN3_${res}_iced.matrix \
$ddir/PN5_${res}_iced.matrix \
$ddir/e2cell_${res}_iced.matrix \
$ddir/l2cell_${res}_iced.matrix \
$ddir/8cell_${res}_iced.matrix \
$ddir/ICM_${res}_iced.matrix \
$ddir/mESC500_${res}_iced.matrix \
-tri 1 -bed $ConfigHP/${res}_mm10.bed -n sperm MII PN3 PN5 e2cell l2cell 8cell ICM mESC500 \
-chr $i -o $wdir/$i -ext png -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 &
done




#################################################
bash ~/codes/hic_plotter.sh -w ~/workspace/8.NT-HiC/4.HiCPlotter/10.20180412 \
-d ~/workspace/8.NT-HiC/d.HiC-Pro_analysis/1.merge_samples_except_20180411/4.norm_to_depth \
-c chr19 -r 500000 -ext pdf 

bash ~/codes/hic_plotter.sh -w ~/workspace/8.NT-HiC/4.HiCPlotter/10.20180412 \
-d ~/workspace/8.NT-HiC/d.HiC-Pro_analysis/1.merge_samples_except_20180411/4.norm_to_depth \
-c chr2 -s 100 -e 700 -r 40000 -ext pdf 

####################20180413 tad di is:
alias plot="bash ~/codes/hic_plotter_tad.sh -w ~/workspace/8.NT-HiC/4.HiCPlotter/11.04_tad \
-d ~/workspace/8.NT-HiC/d.HiC-Pro_analysis/1.merge_samples_except_20180411/4.norm_to_depth \
-tad ~/workspace/8.NT-HiC/g.DI_ALL/1.all_except_norm_depth_20180412/plot_tads \
-di ~/workspace/8.NT-HiC/g.DI_ALL/1.all_except_norm_depth_20180412/directional_index \
-is ~/workspace/8.NT-HiC/f.IS_ALL/1.all_except_norm_depth_20180412/cat_insulation"
plot -c chr2 -s 100 -e 700 -r 100000 &
plot -c chr2 -s 100 -e 700 -r 40000 &



##########xiewei
alias plot="bash ~/codes/hic_plotter_xw_tad.sh -w ~/workspace/8.NT-HiC/4.HiCPlotter/11.04_tad \
-d ~/workspace/8.NT-HiC/3.public_data/4.HiC-Pro_xw/4.norm_to_depth \
-tad ~/workspace/8.NT-HiC/3.public_data/9.DI_xw/plot_tads \
-di ~/workspace/8.NT-HiC/3.public_data/9.DI_xw/directional_index \
-is ~/workspace/8.NT-HiC/3.public_data/7.insulation_score/cat_insulation"
plot -c chr19 -r 500000 &
plot -c chr2 -s 100 -e 700 -r 100000 &
plot -c chr2 -s 100 -e 700 -r 40000 &
 
###########norm test
alias plot="bash ~/codes/hic_plotter_tad.sh -w ~/workspace/8.NT-HiC/4.HiCPlotter/12.norm_test \
-d ~/workspace/8.NT-HiC/d.HiC-Pro_analysis/1.merge_samples_except_20180411/8.norm_test \
-tad ~/workspace/8.NT-HiC/g.DI_ALL/1.all_except_norm_depth_20180412/plot_tads \
-di ~/workspace/8.NT-HiC/g.DI_ALL/1.all_except_norm_depth_20180412/directional_index \
-is ~/workspace/8.NT-HiC/f.IS_ALL/1.all_except_norm_depth_20180412/cat_insulation"
plot -c chr2 -s 100 -e 700 -r 40000 &
plot -c chr2 -s 100 -e 700 -r 100000 &
plot -c chr2 -r 500000 &


####################################################
########2018年5月13日
ddir=~/workspace/8.NT-HiC/5.maps/4.except/6.norm_to_mean
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/13.20180513/1.norm_to_mean

bash ~/codes/hic_plotter.sh -w $wdir -d $ddir \
-c chr7 -s 0 -e 400 -r 40000 -ext png &
bash ~/codes/hic_plotter.sh -w $wdir -d $ddir \
-c chr2 -s 100 -e 600 -r 40000 -ext png &

ddir=~/workspace/8.NT-HiC/5.maps/4.except/6.norm_to_mean
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/13.20180513/2.norm_to_1e8

bash ~/codes/hic_plotter.sh -w $wdir -d $ddir \
-c chr7 -s 200 -e 600 -r 40000 -ext png &
bash ~/codes/hic_plotter.sh -w $wdir -d $ddir \
-c chr2 -s 100 -e 600 -r 40000 -ext png &

####################20180413 tad di is:
alias plot="bash ~/codes/hic_plotter_tad.sh -w ~/workspace/8.NT-HiC/4.HiCPlotter/14.20150515_tad \
-d ~/workspace/8.NT-HiC/5.maps/4.except/6.norm_to_mean \
-tad ~/workspace/8.NT-HiC/g.DI_ALL/2.all_except_norm_mean_20180514/3.plot_tads \
-di ~/workspace/8.NT-HiC/g.DI_ALL/2.all_except_norm_mean_20180514/4.directional_index \
-is ~/workspace/8.NT-HiC/f.IS_ALL/3.all_except_norm_mean_20180514/cat_insulation"
plot -c chr7 -s 100 -e 400 -r 40000 &
plot -c chr7 -s 100 -e 700 -r 100000 &
plot -c chr19 -r 500000 &
plot -c chr2 -s 100 -e 700 -r 40000 &
plot -c chr2 -s 100 -e 700 -r 100000 &


####################不要TAD plot is
alias plot="bash ~/codes/hic_plotter_is.sh -w ~/workspace/8.NT-HiC/4.HiCPlotter/15.20150522_tad \
-d ~/workspace/8.NT-HiC/5.maps/4.except/6.norm_to_mean \
-di ~/workspace/8.NT-HiC/g.DI_ALL/6.except_norm_mean_20180521/4.directional_index \
-is ~/workspace/8.NT-HiC/f.IS_ALL/3.all_except_norm_mean_20180514/cat_insulation"
plot -c chr7 -s 100 -e 500 -r 40000 &
plot -c chr7 -s 100 -e 500 -r 100000 &
plot -c chr7 -r 500000 &


plot -c chr19 -s 100 -e 300 -r 40000 &
plot -c chr19 -s 100 -e 300 -r 100000 &
plot -c chr19 -r 500000 &


########################################################
######2018年8月2日
########################################################

alias plot="bash ~/codes/hic_plotter_is.sh -w ~/workspace/8.NT-HiC/4.HiCPlotter/a.20180802_tad \
-d ~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M \
-di ~/workspace/8.NT-HiC/g.DI_ALL/10.except_norm_1e7_20180727/4.directional_index \
-is ~/workspace/8.NT-HiC/f.IS_ALL/7.except_norm_1e8_20180801/cat_insulation"
plot -c chr19 -r 100000 -ext pdf &

#每个sample单独画，作图
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/c.20180813_each_sample
di=~/workspace/8.NT-HiC/g.DI_ALL/11.except_norm_1e8_20180813/4.directional_index
is=~/workspace/8.NT-HiC/f.IS_ALL/12.except_norm_1e8_20180813/cat_insulation
mkdir -p $wdir
res=40000
for i in ${keys[@]}; do 
python ~/codes/HiCPlotter.py -f $ddir/${i}_${res}_iced.matrix \
-hist $di/$i.bedGraph,$is/$i.bedGraph -hl Directional_Index,Insulation_Score \
-bed $ConfigHP/${res}_mm10.bed -n $i -tri 1 \
-chr chr19 -o $wdir/$i -ext pdf -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 -s 875 -e 1125 &
done

#all DI IS:
for i in ${keys[@]}; do  echo -n "\$ddir/${i}_${res}_iced.matrix "; done
for i in ${keys[@]}; do  echo -n "\$di/$i.bedGraph,\$is/$i.bedGraph "; done
for i in ${keys[@]}; do  echo -n "DI,IS "; done
for i in ${keys[@]}; do  echo -n "$i "; done
res=40000
python ~/codes/HiCPlotter.py -f `for i in ${keys[@]}; do  echo -n "\$ddir/${i}_${res}_iced.matrix "; done` \
-hist `for i in ${keys[@]}; do  echo -n "\$di/$i.bedGraph,\$is/$i.bedGraph "; done` \
-hl `for i in ${keys[@]}; do  echo -n "DI,IS "; done` \
-n `for i in ${keys[@]}; do  echo -n "$i "; done` \
-bed $ConfigHP/${res}_mm10.bed -tri 1 \
-chr chr19 -o $wdir/DI-IS -ext pdf -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 -s 875 -e 1125 &

res=100000
python ~/codes/HiCPlotter.py -f `for i in ${keys[@]}; do  echo -n "\$ddir/${i}_${res}_iced.matrix "; done` \
-hist `for i in ${keys[@]}; do  echo -n "\$di/$i.bedGraph,\$is/$i.bedGraph "; done` \
-hl `for i in ${keys[@]}; do  echo -n "DI,IS "; done` \
-n `for i in ${keys[@]}; do  echo -n "$i "; done` \
-bed $ConfigHP/${res}_mm10.bed -tri 1 \
-chr chr19 -o $wdir/DI-IS -ext pdf -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 &


##########whole genome
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/e.20180908_wholeGenome
mkdir -p $wdir
res=1000000
for i in ${keys[@]}; do 
python ~/codes/HiCPlotter.py -f $ddir/${i}_${res}_iced.matrix \
--bedFile $ConfigHP/${res}_mm10.bed -n $i --tripleColumn 1 --wholeGenome 1 \
-chr chr19 -o $wdir/test1_$i -ext pdf -dpi 360 -r ${res} --spine 1 --heatmapColor 1 --matrixMax 15 &
done
#test1 --matrixMax 100

##############chr14 
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/f.20180913_chr4
di=~/workspace/8.NT-HiC/g.DI_ALL/2.except_norm_100M_20180825/4.directional_index
is=~/workspace/8.NT-HiC/f.IS_ALL/12.except_norm_1e8_20180813/cat_insulation
mkdir -p $wdir
res=100000
for i in ${keys[@]}; do 
python ~/codes/HiCPlotter.py -f $ddir/${i}_${res}_iced.matrix \
-hist $di/$i.bedGraph,$is/$i.bedGraph -hl Directional_Index,Insulation_Score \
-bed $ConfigHP/${res}_mm10.bed -n $i -tri 1 \
-chr chr4 -o $wdir/$i -ext pdf -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 &
done
for i in ${keys[@]}; do 
python ~/codes/HiCPlotter.py -f $ddir/${i}_${res}_iced.matrix \
-hist $di/$i.bedGraph,$is/$i.bedGraph -hl Directional_Index,Insulation_Score \
-bed $ConfigHP/${res}_mm10.bed -n $i -tri 1 \
-chr chr4 -o $wdir/$i -ext pdf -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 -s  -e &
done

res=40000
for i in ${keys[@]}; do 
python ~/codes/HiCPlotter.py -f $ddir/${i}_${res}_iced.matrix \
-hist $di/$i.bedGraph,$is/$i.bedGraph -hl Directional_Index,Insulation_Score \
-bed $ConfigHP/${res}_mm10.bed -n $i -tri 1 \
-chr chr4 -o $wdir/$i -ext pdf -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 -s 2250 -e 2500 &
done


###chr6 122,280,000-122,380,000
res=10000
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/
python ~/codes/HiCPlotter.py -f `for i in ${keys[@]}; do  echo -n "\$ddir/${i}_${res}_iced.matrix "; done` \
-bed $ConfigHP/${res}_mm10.bed -n `for i in ${keys[@]}; do  echo -n "$i "; done` -tri 1 \
-chr chr6 -o $wdir/chr6-Phc1 -ext pdf -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 -s 12228 -e 12238 &


##8-cell rep1，rep2，rep3；chr19 1-61Mb res=100kb；chr19 35-45Mb res=40kb。

ddir=~/workspace/8.NT-HiC/5.maps/2.replicate/4.norm_to_100M/
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/g.8cell
mkdir -p $wdir

for key in 8cell 2cell; do 
res=40000
python ~/codes/HiCPlotter.py -f $ddir/${key}_rep*_${res}_iced.matrix \
-bed $ConfigHP/${res}_mm10.bed -n rep1 rep2 rep3 -tri 1 \
-chr chr19 -o $wdir/${key} -ext pdf -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 -s 875 -e 1125 &
res=100000
python ~/codes/HiCPlotter.py -f $ddir/${key}_rep[123]_${res}_iced.matrix \
-bed $ConfigHP/${res}_mm10.bed -n rep1 rep2 rep3 -tri 1 \
-chr chr19 -o $wdir/${key} -ext pdf -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 &
done



###e2cell rep1，rep2，rep3；chr19 1-61Mb res=100kb；chr19 35-45Mb res=40kb。
key=icm
res=40000
python ~/codes/HiCPlotter.py -f $ddir/${key}_rep*_${res}_iced.matrix \
-bed $ConfigHP/${res}_mm10.bed -n rep1 rep2 -tri 1 \
-chr chr19 -o $wdir/${key} -ext pdf -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 -s 875 -e 1125 &
res=100000
python ~/codes/HiCPlotter.py -f $ddir/${key}_rep*_${res}_iced.matrix \
-bed $ConfigHP/${res}_mm10.bed -n rep1 rep2 -tri 1 \
-chr chr19 -o $wdir/${key} -ext pdf -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 &


#for CMM change
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
ddir2=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/6.norm_to_100M
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/h.cmm
res=40000
mkdir -p $wdir

python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir2/ICM_${res}_iced.matrix $ddir/icm_${res}_iced.matrix \
-bed $ConfigHP/${res}_mm10.bed -n cc ICM icm -tri 1 \
-chr chr15 -o $wdir/cmm -ext pdf -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 -s 275 -e 400 &




chr15:26,672,807-29,403,952

ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
ddir2=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/6.norm_to_100M
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/h.cmm
res=40000
#for i in chr1:3800-3925 chr2:250-375 chr4:1225-1350 chr4:2875-3000 chr5:3325-3450; do chr=${i%%:*};j=${i##*:};s=${j%%-*};e=${j##*-}
for i in chr8:1950-2075 chr9:0-250 chr10:2050-2175 chr15:500-750; do chr=${i%%:*};j=${i##*:};s=${j%%-*};e=${j##*-}
python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir2/ICM_${res}_iced.matrix $ddir/icm_${res}_iced.matrix \
-bed $ConfigHP/${res}_mm10.bed -n cc ICM icm -tri 1 \
-chr $chr -o $wdir/cmm -ext pdf -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 -s $s -e $e &
done


#############test TAD in 
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/i.for_GO
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
ddir2=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/6.norm_to_100M
tdir=~/workspace/8.NT-HiC/g.DI_ALL/2.except_norm_100M_20180825/3.plot_tads
didir=~/workspace/8.NT-HiC/g.DI_ALL/2.except_norm_100M_20180825/4.directional_index
isdir=~/workspace/8.NT-HiC/f.IS_ALL/1.except_res40k_is1M_ids200k_nt025/2.cat_insulation

res=40000
python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir/6h_${res}_iced.matrix $ddir/12h_${res}_iced.matrix \
$ddir/e2cell_${res}_iced.matrix $ddir/l2cell_${res}_iced.matrix $ddir/8cell_${res}_iced.matrix $ddir/icm_${res}_iced.matrix \
$ddir2/PN3_zygote_${res}_iced.matrix $ddir2/PN5_zygote_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix $ddir2/late_2cell_${res}_iced.matrix \
$ddir2/8cell_${res}_iced.matrix $ddir2/ICM_${res}_iced.matrix \
-bed $ConfigHP/${res}_mm10.bed -n cc 6h 12h e2cell l2cell 8cell icm PN3 PN5 e2 l2 8c ICM -tri 1 \
-chr chr7 -s 3500 -e 3625 -o $wdir/nt-nf -ext pdf -dpi 360 -r 40000 -spi 1 -hmc 1 -ptr 1 &

res=40000
python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir/6h_${res}_iced.matrix $ddir/12h_${res}_iced.matrix \
$ddir/e2cell_${res}_iced.matrix $ddir/l2cell_${res}_iced.matrix $ddir/8cell_${res}_iced.matrix $ddir/icm_${res}_iced.matrix \
$ddir2/PN3_zygote_${res}_iced.matrix $ddir2/PN5_zygote_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix $ddir2/late_2cell_${res}_iced.matrix \
$ddir2/8cell_${res}_iced.matrix $ddir2/ICM_${res}_iced.matrix \
-bed $ConfigHP/${res}_mm10.bed -n cc 6h 12h e2cell l2cell 8cell icm PN3 PN5 e2 l2 8c ICM -tri 1 \
-chr chr7 -s 750 -e 875 -o $wdir/nt-nf -ext pdf -dpi 360 -r 40000 -spi 1 -hmc 1 -ptr 1 -trh 50 &
#chrX


res=500000
python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir/6h_${res}_iced.matrix $ddir/12h_${res}_iced.matrix \
$ddir/e2cell_${res}_iced.matrix $ddir/l2cell_${res}_iced.matrix $ddir/8cell_${res}_iced.matrix $ddir/icm_${res}_iced.matrix \
$ddir2/PN3_zygote_${res}_iced.matrix $ddir2/PN5_zygote_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix $ddir2/late_2cell_${res}_iced.matrix \
$ddir2/8cell_${res}_iced.matrix $ddir2/ICM_${res}_iced.matrix \
-bed $ConfigHP/${res}_mm10.bed -n cc 6h 12h e2cell l2cell 8cell icm PN3 PN5 e2 l2 8c ICM -tri 1 \
-chr chrX -o $wdir/nt-nf -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &

res=40000
python ~/codes/HiCPlotter.py -f $ddir2/ICM_${res}_iced.matrix \
-bed $ConfigHP/${res}_mm10.bed -n ICM -tri 1 \
-pcd 1 \
-pcdf ~/workspace/8.NT-HiC/g.DI_ALL/a.Consolidation_analysis/6.NF_20180803/1.NF/ICM.tad \
-chr chr7 -s 500 -e 625 -o $wdir/ICM -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &

isdir=~/workspace/8.NT-HiC/f.IS_ALL/1.except_res40k_is1M_ids200k_nt025/2.cat_insulation
isdir2=~/workspace/8.NT-HiC/3.public_data/f.XW_analysis/9.IS_res40k_is1M_ids200k_nt025/2.cat_insulation
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/m.loops
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
ddir2=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/6.norm_to_100M
ddir3=~/workspace/8.NT-HiC/3.public_data/a.HiC-Pro_lj/6.norm_to_100M
res=100000
python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir/e2cell_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix \
$ddir/l2cell_${res}_iced.matrix $ddir2/late_2cell_${res}_iced.matrix -bed $ConfigHP/${res}_mm10.bed -n cc NT-e2 NF-e2 NT-l2 NF-l2 \
-tri 1 -chr chr7 -s 100 -e 140 -o $wdir/zscan4d-ep -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &
res=10000
python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir/e2cell_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix $ddir/l2cell_${res}_iced.matrix $ddir2/late_2cell_${res}_iced.matrix -bed $ConfigHP/${res}_mm10.bed -n cc
NT-e2 NF-e2 NT-l2 NF-l2 -tri 1 -chr chr7 -s 1000 -e 1500 -o $wdir/zscan4d-ep -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &
res=10000
python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir/e2cell_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix $ddir/l2cell_${res}_iced.matrix $ddir2/late_2cell_${res}_iced.matrix -bed $ConfigHP/${res}_mm10.bed -n cc
NT-e2 NF-e2 NT-l2 NF-l2 -tri 1 -chr chr7 -s 15500 -e 16000 -o $wdir/Yy2-ep -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &
res=25000
python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir/e2cell_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix $ddir/l2cell_${res}_iced.matrix $ddir2/late_2cell_${res}_iced.matrix -bed $ConfigHP/${res}_mm10.bed -n cc
NT-e2 NF-e2 NT-l2 NF-l2 -tri 1 -chr chr7 -s 400 -e 600 -o $wdir/zscan4d-ep -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &
res=40000
python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir/e2cell_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix $ddir/l2cell_${res}_iced.matrix $ddir2/late_2cell_${res}_iced.matrix -bed $ConfigHP/${res}_mm10.bed -n cc
NT-e2 NF-e2 NT-l2 NF-l2 -tri 1 -chr chr7 -s 250 -e 375 -o $wdir/zscan4d-ep -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &

res=100000
python ~/codes/HiCPlotter.py -f $ddir/e2cell_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix -n NT-e2 NF-e2 -bed $ConfigHP/${res}_mm10.bed \
-tri 1 -chr chr16 -s 160 -e 240 -o $wdir/St6gal1-ep -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &

#加LJ的2cell
res=100000
python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir/e2cell_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix $ddir3/2cell_${res}_iced.matrix -bed $ConfigHP/${res}_mm10.bed -n cc NT-e2 XW-e2 LJ-2cell -tri 1 -chr chr7 -s 100 -e 140 -o $wdir/zscan4d-ep-3 -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &


python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir/e2cell_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix $ddir3/2cell_${res}_iced.matrix -bed $ConfigHP/${res}_mm10.bed -n cc NT-e2 XW-e2 LJ-2cell -tri 1 -chr chr7 -s 100 -e 140 -o $wdir/zscan4d-ep-mm$i -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 -mm $i &

#分rep
ddir=~/workspace/8.NT-HiC/5.maps/2.replicate/2.iced
res=100000
python ~/codes/HiCPlotter.py -f $ddir/e2cell_*_${res}_iced.matrix -bed $ConfigHP/${res}_mm10.bed -n rep1 rep2 rep3 -tri 1 -chr chr7 -s 100 -e 140 -o $wdir/zscan4d-ep-NT-e2-reps-unnorm -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &
python ~/codes/HiCPlotter.py -f $ddir/cc_*_${res}_iced.matrix -bed $ConfigHP/${res}_mm10.bed -n rep1 rep2 rep3 rep4 -tri 1 -chr chr7 -s 100 -e 140 -o $wdir/zscan4d-ep-NT-CC-reps-unnorm -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &

ddir1=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/7.replicate/4.iced
res=100000
python ~/codes/HiCPlotter.py -f $ddir1/early_2cell*_${res}_iced.matrix -bed $ConfigHP/${res}_mm10.bed -n rep1 rep2 rep3 rep4 rep5 -tri 1 -chr chr7 -s 100 -e 140 -o $wdir/zscan4d-ep-xw-e2-reps-unnorm -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &

######################ICM loops
isdir=~/workspace/8.NT-HiC/f.IS_ALL/1.except_res40k_is1M_ids200k_nt025/2.cat_insulation
isdir2=~/workspace/8.NT-HiC/3.public_data/f.XW_analysis/9.IS_res40k_is1M_ids200k_nt025/2.cat_insulation
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/m.loops
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
ddir2=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/6.norm_to_100M
res=100000
python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir/e2cell_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix $ddir/icm_${res}_iced.matrix $ddir2/ICM_${res}_iced.matrix -bed $ConfigHP/${res}_mm10.bed -n cc NT-e2 NF-e2 NT-ICM NF-ICM -tri 1 -chr chr16 -s 135 -e 945 -o $wdir/Chaf1b-ep -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &

python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir/e2cell_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix $ddir/icm_${res}_iced.matrix $ddir2/ICM_${res}_iced.matrix -bed $ConfigHP/${res}_mm10.bed -n cc NT-e2 NF-e2 NT-ICM NF-ICM -tri 1 -chr chr10 -s 600 -e 700 -o $wdir/Tet1-ep -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 -trh 1000 &

python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir/e2cell_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix $ddir/icm_${res}_iced.matrix $ddir2/ICM_${res}_iced.matrix -bed $ConfigHP/${res}_mm10.bed -n cc NT-e2 NF-e2 NT-ICM NF-ICM -tri 1 -chr chr19 -s 68 -e 103 -o $wdir/ICM-loop -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 -trh 10000 -hist $isdir/cc.bedGraph $isdir/e2cell.bedGraph $isdir2/early_2cell.bedGraph $isdir/icm.bedGraph $isdir2/ICM.bedGraph -hl IS IS IS IS IS &

python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir/e2cell_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix $ddir/icm_${res}_iced.matrix $ddir2/ICM_${res}_iced.matrix -bed $ConfigHP/${res}_mm10.bed -n cc NT-e2 NF-e2 NT-ICM NF-ICM -tri 1 -chr chr4 -s 1370 -e 1380 -o $wdir/ICM-loop -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 -trh 10000 -hist $isdir/cc.bedGraph $isdir/e2cell.bedGraph $isdir2/early_2cell.bedGraph $isdir/icm.bedGraph $isdir2/ICM.bedGraph -hl IS IS IS IS IS &

python ~/codes/HiCPlotter.py -f $ddir/e2cell_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix -bed $ConfigHP/${res}_mm10.bed -a $wdir/NF-loops.bed $wdir/NF-loops.bed -al NF-e2 NF-e2 -n NT-e2 NF-e2  -tri 1 -chrchr7 -s 100 -e 140 -o $wdir/zscan4d-ep2 -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &

awk '($4=="chr19-7568"||$4=="chr19-6041")&&($6>7000000&&$7<9000000){print $1,($2+$3)/2,($6+$7)/2 > $9".loops"}' ~/workspace/8.NT-HiC/k.fithic/2.xw_100k_20181202/ICM/d.SE_pairs_normalized/all.SEP
python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir/e2cell_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix $ddir/icm_${res}_iced.matrix $ddir2/ICM_${res}_iced.matrix -bed $ConfigHP/${res}_mm10.bed -n cc NT-e2 NF-e2 NT-ICM NF-ICM -tri 1 -chr chr19 -s 70 -e 90 -o $wdir/ICM-loop -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 -trh 10000 -hist $isdir/cc.bedGraph $isdir/e2cell.bedGraph $isdir2/early_2cell.bedGraph $isdir/icm.bedGraph $isdir2/ICM.bedGraph -hl IS IS IS IS IS -a $wdir/both.loops $wdir/NF.loops $wdir/NF.loops $wdir/NF.loops $wdir/NF.loops -al NF NF NF NF NF &
###2cell loops

res=100000
python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir/e2cell_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix $ddir/l2cell_${res}_iced.matrix $ddir2/late_2cell_${res}_iced.matrix -bed $ConfigHP/${res}_mm10.bed -n cc NT-e2 NF-e2 NT-l2 NF-l2 -tri 1 -chr chr13 -s 540 -e 590 -o $wdir/2cell-loop -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 -trh 1000 -hist $isdir/cc.bedGraph $isdir/e2cell.bedGraph $isdir2/early_2cell.bedGraph $isdir/l2cell.bedGraph $isdir2/late_2cell.bedGraph -hl IS IS IS IS IS &

python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir/e2cell_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix $ddir/l2cell_${res}_iced.matrix $ddir2/late_2cell_${res}_iced.matrix -bed $ConfigHP/${res}_mm10.bed -n cc NT-e2 NF-e2 NT-l2 NF-l2 -tri 1 -chr chr16 -s 100 -e 140 -o $wdir/2cell-loop -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 -trh 1000 -hist $isdir/cc.bedGraph $isdir/e2cell.bedGraph $isdir2/early_2cell.bedGraph $isdir/l2cell.bedGraph $isdir2/late_2cell.bedGraph -hl IS IS IS IS IS &

python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir/e2cell_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix $ddir/l2cell_${res}_iced.matrix $ddir2/late_2cell_${res}_iced.matrix -bed $ConfigHP/${res}_mm10.bed -n cc NT-e2 NF-e2 NT-l2 NF-l2 -tri 1 -chr chr2 -s 1100 -e 1150 -o $wdir/2cell-loop -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 -trh 1000 &

awk '($4=="chr7-13538"||$4=="chr7-44898")&&($6>7000000&&$7<14000000){printf("%s\t%d\t%d\n",$1,($2+$3)/2,($6+$7)/2) > $9".loops"}' ~/workspace/8.NT-HiC/k.fithic/2.xw_100k_20181202/early_2cell/d.SE_pairs_normalized/all.SEP
python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir/e2cell_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix $ddir/l2cell_${res}_iced.matrix $ddir2/late_2cell_${res}_iced.matrix -bed $ConfigHP/${res}_mm10.bed -n cc NT-e2 NF-e2 NT-l2 NF-l2 -tri 1 -chr chr7 -s 100 -e 140 -o $wdir/zscan4d-ep -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 -trh 1000 -hist $isdir/cc.bedGraph $isdir/e2cell.bedGraph $isdir2/early_2cell.bedGraph $isdir/l2cell.bedGraph $isdir2/late_2cell.bedGraph -hl IS IS IS IS IS -hm 1 1 1 1 1 -a $wdir/NF.loops $wdir/NF.loops $wdir/NF.loops $wdir/NF.loops $wdir/NF.loops -al NF NF NF NF NF &



####################delta RTI
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/i.for_GO
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
ddir2=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/6.norm_to_100M

res=40000
pos="-chr chr16 -s 200 -e 300"
python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir/6h_${res}_iced.matrix $ddir2/PN3_zygote_${res}_iced.matrix \
$ddir/12h_${res}_iced.matrix $ddir2/PN5_zygote_${res}_iced.matrix \
$ddir/e2cell_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix \
$ddir/l2cell_${res}_iced.matrix $ddir2/late_2cell_${res}_iced.matrix $ddir/icm_${res}_iced.matrix $ddir2/ICM_${res}_iced.matrix \
-bed $ConfigHP/${res}_mm10.bed -n cc NT-6h NF-PN3 NT-12h NF-PN5 NT-e2 NF-e2 NT-l2 NF-l2 NT-ICM NF-ICM -tri 1 $pos -o $wdir/for-delta-RTI \
-ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 -trh 100 &



##############################2019年9月6日
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/n.kdm4d
ddir=~/workspace/8.NT-HiC/5.maps/2.replicate/4.norm_to_100M
res=100000
python ~/codes/HiCPlotter.py -f $ddir/*_${res}_iced.matrix -bed $ConfigHP/${res}_mm10.bed -n 6h_rep6 e2cell_rep4 kdm4d-e2cell_rep1 kdm4d-e2cell_rep2 sertoli-e2cell_rep1 -tri 1 -chr chr16 -s 80 -e 120 -o $wdir/reps -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &

wdir=~/workspace/8.NT-HiC/4.HiCPlotter/n.kdm4d
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
ddir2=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/6.norm_to_100M
res=100000
python ~/codes/HiCPlotter.py -f $ddir/cc_${res}_iced.matrix $ddir/e2cell_${res}_iced.matrix $ddir/kdm4d-e2cell_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix -bed $ConfigHP/${res}_mm10.bed -n cc NT-e2 kdm4d-e2 NF-e2 \
-tri 1 -chr chr7 -s 100 -e 140 -o $wdir/zscan4d-ep -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &

python ~/codes/HiCPlotter.py -f $ddir/e2cell_${res}_iced.matrix $ddir/kdm4d-e2cell_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix $ddir/partheno-e2cell_${res}_iced.matrix -bed $ConfigHP/${res}_mm10.bed -n NT kdm4d NF partheno \
-tri 1 -chr chr7 -s 100 -e 140 -o $wdir/zscan4d-ep-new -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &

python ~/codes/HiCPlotter.py -f $ddir/e2cell_${res}_iced.matrix $ddir/kdm4d-e2cell_${res}_iced.matrix $ddir/tsa-e2cell_${res}_iced.matrix \
$ddir2/early_2cell_${res}_iced.matrix $ddir/partheno-e2cell_${res}_iced.matrix $ddir/partheno-1hpa_${res}_iced.matrix \
-bed $ConfigHP/${res}_mm10.bed -n NT kdm4d tsa NF partheno-e2cell partheno-1hpa \
-tri 1 -chr chr7 -s 100 -e 140 -o $wdir/zscan4d-ep-new -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &

python ~/codes/HiCPlotter.py -f $ddir/e2cell_${res}_iced.matrix $ddir/kdm4d-e2cell_${res}_iced.matrix $ddir/tsa-e2cell_${res}_iced.matrix \
$ddir2/early_2cell_${res}_iced.matrix $ddir/partheno-e2cell_${res}_iced.matrix \
-bed $ConfigHP/${res}_mm10.bed -n NT kdm4d tsa NF partheno-e2cell \
-tri 1 -chr chr16 -s 80 -e 120 -o $wdir/eliminate-TAD -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &

wdir=~/workspace/8.NT-HiC/4.HiCPlotter/n.kdm4d
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
ddir2=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/6.norm_to_100M
res=100000
python ~/codes/HiCPlotter.py -f $ddir/e2cell_${res}_iced.matrix $ddir/sertoli-e2cell_${res}_iced.matrix $ddir/kdm4d-e2cell_${res}_iced.matrix $ddir/tsa-e2cell_${res}_iced.matrix \
$ddir2/early_2cell_${res}_iced.matrix $ddir/partheno-e2cell_${res}_iced.matrix \
-bed $ConfigHP/${res}_mm10.bed -n NT sertoli kdm4d tsa NF partheno  \
-tri 1 -chr chr5 -o $wdir/all-e2cell -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &
python ~/codes/HiCPlotter.py -f $ddir/e2cell_${res}_iced.matrix $ddir/kdm4d-e2cell_${res}_iced.matrix $ddir/tsa-e2cell_${res}_iced.matrix \
$ddir2/early_2cell_${res}_iced.matrix \
-bed $ConfigHP/${res}_mm10.bed -n NT kdm4d tsa NF \
-tri 1 -chr chr15 -o $wdir/all-e2cell -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &
python ~/codes/HiCPlotter.py -f $ddir/e2cell_${res}_iced.matrix $ddir/kdm4d-e2cell_${res}_iced.matrix $ddir/tsa-e2cell_${res}_iced.matrix \
$ddir2/early_2cell_${res}_iced.matrix \
-bed $ConfigHP/${res}_mm10.bed -n NT kdm4d tsa NF \
-tri 1 -chr chr15 -s 800 -e 1040 -o $wdir/all-e2cell -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &

res=500000
python ~/codes/HiCPlotter.py -f $ddir/e2cell_${res}_iced.matrix $ddir/sertoli-e2cell_${res}_iced.matrix $ddir/kdm4d-e2cell_${res}_iced.matrix $ddir/tsa-e2cell_${res}_iced.matrix \
$ddir2/early_2cell_${res}_iced.matrix $ddir/partheno-e2cell_${res}_iced.matrix \
-bed $ConfigHP/${res}_mm10.bed -n NT sertoli kdm4d tsa NF partheno  \
-tri 1 -chr chr5 -o $wdir/all-e2cell -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &



##################partheno chr19
### chr19 1-61Mb res=100kb；chr19 35-45Mb res=40kb。
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/n.kdm4d
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
key=partheno
res=40000
python ~/codes/HiCPlotter.py -f $ddir/partheno-1hpa_${res}_iced.matrix $ddir/partheno-e2cell_${res}_iced.matrix \
-bed $ConfigHP/${res}_mm10.bed -n 1hpa e2cell -tri 1 \
-chr chr19 -o $wdir/${key} -ext pdf -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 -s 875 -e 1125 &
res=100000
python ~/codes/HiCPlotter.py -f $ddir/partheno-1hpa_${res}_iced.matrix $ddir/partheno-e2cell_${res}_iced.matrix \
-bed $ConfigHP/${res}_mm10.bed -n 1hpa e2cell -tri 1 \
-chr chr19 -o $wdir/${key} -ext pdf -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 &


wdir=~/workspace/8.NT-HiC/4.HiCPlotter/n.kdm4d
ddir=~/workspace/8.NT-HiC/5.maps/2.replicate/2.iced
res=40000
python ~/codes/HiCPlotter.py -f $ddir/6h_rep3_${res}_iced.matrix $ddir/6h_rep5_${res}_iced.matrix $ddir/6h_rep6_${res}_iced.matrix \
-bed $ConfigHP/${res}_mm10.bed -n 6h_rep3 6h_rep5 6h_rep6 -tri 1 \
-chr chr19 -o $wdir/6h -ext pdf -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 -s 875 -e 1125 &

#ES zscan4d
ddir=~/workspace/8.NT-HiC/e.ES/0.links/4.norm_to_100M
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/n.kdm4d
res=100000
python ~/codes/HiCPlotter.py -f $ddir/xiewei_${res}_iced.matrix $ddir/es500_rep1_${res}_iced.matrix $ddir/es500_rep2_${res}_iced.matrix \
-bed $ConfigHP/${res}_mm10.bed -n xiewei es500-1 es500-2 \
-tri 1 -chr chr7 -s 100 -e 140 -o $wdir/es-zscan4d-ep -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &

#dux in all 2cell 
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/n.kdm4d
ddir=~/workspace/8.NT-HiC/5.maps/4.except/4.norm_to_100M
ddir2=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/6.norm_to_100M
res=100000
python ~/codes/HiCPlotter.py -f $ddir/e2cell_${res}_iced.matrix $ddir/sertoli-e2cell_${res}_iced.matrix \
$ddir/kdm4d-e2cell_${res}_iced.matrix $ddir/tsa-e2cell_${res}_iced.matrix $ddir2/early_2cell_${res}_iced.matrix $ddir/partheno-e2cell_${res}_iced.matrix \
-bed $ConfigHP/${res}_mm10.bed -n NT sertoli kdm4d tsa NF partheno  \
-tri 1 -chr chr10 -s 550 -e 610 -o $wdir/dux-e2cell -ext pdf -dpi 360 -r $res -spi 1 -hmc 1 -ptr 1 &




########################
ddir=~/workspace/8.NT-HiC/3.public_data/c.HiC-Pro_xw/6.norm_to_100M
wdir=~/workspace/8.NT-HiC/4.HiCPlotter/8.xw_development_100k_20180127/1.Zscan4d
di=~/workspace/8.NT-HiC/3.public_data/f.XW_analysis/5.DI/4.directional_index
is=~/workspace/8.NT-HiC/3.public_data/f.XW_analysis/4.IS/2.cat_insulation

keys=(early_2cell late_2cell ICM)
res=40000
python ~/codes/HiCPlotter.py -f `for i in ${keys[@]}; do  echo -n "\$ddir/${i}_${res}_iced.matrix "; done` \
-hist `for i in ${keys[@]}; do  echo -n "\$di/$i.bedGraph,\$is/$i.bedGraph "; done` \
-hl `for i in ${keys[@]}; do  echo -n "DI,IS "; done` \
-n `for i in ${keys[@]}; do  echo -n "$i "; done` \
-bed $ConfigHP/${res}_mm10.bed -tri 1 \
-chr chr7 -s 250 -e 350 -o $wdir/zscan4d -ext pdf -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 &
res=100000
python ~/codes/HiCPlotter.py -f `for i in ${keys[@]}; do  echo -n "\$ddir/${i}_${res}_iced.matrix "; done` \
-hist `for i in ${keys[@]}; do  echo -n "\$di/$i.bedGraph,\$is/$i.bedGraph "; done` \
-hl `for i in ${keys[@]}; do  echo -n "DI,IS "; done` \
-n `for i in ${keys[@]}; do  echo -n "$i "; done` \
-bed $ConfigHP/${res}_mm10.bed -tri 1 \
-chr chr7 -s 100 -e 140 -o $wdir/zscan4d -ext pdf -dpi 360 -r ${res} -spi 1 -hmc 1 -ptr 1 &





