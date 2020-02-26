QC=/usr/local/software/HiC-software/3DChromatin_ReplicateQC


selects=("cc_rep1" "cc_rep2" "cc_rep3" "cc_rep4" "05h_rep2" "05h_rep3" "1h_rep1" "1h_rep2" "2h_rep1" "2h_rep2" "6h_rep5" "12h_rep1" "12h_rep2" "12h_rep3" "e2cell_rep3" "l2cell_rep1" "l2cell_rep2" "4cell_rep2" "4cell_rep3" "8cell_rep1" "8cell_rep2" "8cell_rep3" "morula_rep1" "morula_rep3" "morula_rep4" "icm_rep1" "icm_rep2" "te_rep1" "te_rep2" "te_rep3") 
selects=("e2cell_rep4" "6h_rep6" "kdm4d-e2cell_rep1" "kdm4d-e2cell_rep2" "partheno-1hpa_rep1" "partheno-1hpa_rep2" "partheno-e2cell_rep1" "partheno-e2cell_rep2" "sertoli-e2cell_rep1" "tsa-e2cell_rep1") #20190920
selects=("sertoli-e2cell_rep2" "tsa-e2cell_rep2")
wdir=~/workspace/8.NT-HiC/o.HiCRep/4.3DRepQC/1.NT/0.fithic
ddir=~/workspace/8.NT-HiC/5.maps/2.replicate/1.raw
bdir=~/workspace/8.NT-HiC/5.maps/2.replicate/3.biases
res=100000
for key in ${selects[@]}; do
(mkdir -p $wdir/$key/data $wdir/$key/output $wdir/$key/logs
nohup python $UtilsHP/hicpro2fithic.py -i $ddir/${key}_${res}.matrix -s $bdir/${key}_${res}_iced.matrix.biases \
-b $ConfigHP/${res}_mm10.bed -o $wdir/$key/data -r ${res} > $wdir/$key/logs/hicpro2fithic.log )&
done


ddir=~/workspace/8.NT-HiC/o.HiCRep/4.3DRepQC/1.NT/0.fithic/*/data/fithic.interactionCounts.gz
rm metadata.samples
for i in ~/workspace/8.NT-HiC/o.HiCRep/4.3DRepQC/1.NT/0.fithic/*/data/fithic.interactionCounts.gz; do key=${i##*0.fithic/};key=${key%%/data*};
echo -e $key"\t"$i >> metadata.samples
done

selects=("cc_rep1" "cc_rep2" "cc_rep3" "cc_rep4" "05h_rep2" "05h_rep3" "1h_rep1" "1h_rep2" "2h_rep1" "2h_rep2" "6h_rep5" "6h_rep6" "12h_rep1" "12h_rep2" "12h_rep3" "e2cell_rep3" "e2cell_rep4" "l2cell_rep1" "l2cell_rep2" "4cell_rep2" "4cell_rep3" "8cell_rep1" "8cell_rep2" "8cell_rep3" "morula_rep1" "morula_rep3" "morula_rep4" "icm_rep1" "icm_rep2" "te_rep1" "te_rep2" "te_rep3"  "kdm4d-e2cell_rep1" "kdm4d-e2cell_rep2" "partheno-1hpa_rep1" "partheno-1hpa_rep2" "partheno-e2cell_rep1" "partheno-e2cell_rep2" "sertoli-e2cell_rep1" "sertoli-e2cell_rep2" "tsa-e2cell_rep1" "tsa-e2cell_rep2")

for i in ${selects[@]}; do echo $i;done | awk -v FS="_" 'BEGIN{k="";temp[0]} k!=$1{k=$1;delete temp;temp[$0];next} k==$1{for(i in temp){print i,$0};temp[$0]}' > metadata.pairs

awk '{print $1,$2,$3,int(($2+$3)/2)}' $ConfigHP/100000_mm10.bed > $wdir/Bins_100000.bed
gzip $wdir/Bins_100000.bed
wdir=~/workspace/8.NT-HiC/o.HiCRep/4.3DRepQC/1.NT
3DChromatin_ReplicateQC run_all --metadata_samples $wdir/metadata.samples --metadata_pairs $wdir/metadata.pairs --bins $wdir/Bins_100000.bed.gz --outdir $wdir/1.output --parameters_file $wdir/example_parameters.txt


cat <<\EOF > $wdir/parameters2.txt
GenomeDISCO|subsampling	lowest
GenomeDISCO|tmin	1
GenomeDISCO|tmax	5
GenomeDISCO|norm	uniform
GenomeDISCO|scoresByStep	yes
GenomeDISCO|removeDiag	yes
GenomeDISCO|transition	yes
HiCRep|h	2
HiCRep|maxdist	5000000
HiC-Spector|n	10
QuASAR|rebinning	5000000
SGE|text	"-l h_vmem=10G"
slurm|text	"--mem 3G"
EOF
wdir=~/workspace/8.NT-HiC/o.HiCRep/4.3DRepQC/1.NT
for i in {1..4}; do 
nohup 3DChromatin_ReplicateQC run_all --metadata_samples $wdir/metadata.samples --metadata_pairs $wdir/metadata.pairs --bins $wdir/Bins_100000.bed.gz --outdir $wdir/$i.output --parameters_file $wdir/parameters/$i.parameters.txt > $wdir/logs/$i.log &
done

for i in *output/results/reproducibility/pdfs/*.pdf; do k=${i%%.*}
cp $i pdfs/$k.$(basename $i)
done


###each 2 rep do for heatmap:
for i in ${selects[@]}; do echo $i;done | awk '{a[NR]=$1;b=NR} END{for(i=1;i<NR;i++){for(j=i+1;j<=NR;j++){print a[i],a[j]}}}' > metadata.pairs2
wdir=~/workspace/8.NT-HiC/o.HiCRep/4.3DRepQC/1.NT
mkd $wdir/5.each2/logs
for i in {1..4}; do 
nohup 3DChromatin_ReplicateQC run_all --metadata_samples $wdir/metadata.samples --metadata_pairs $wdir/metadata.pairs2 --bins $wdir/Bins_100000.bed.gz --outdir $wdir/5.each2/$i.output --parameters_file $wdir/parameters/$i.parameters.txt > $wdir/5.each2/logs/$i.log &
done


##############################################################################
**GenomeDISCO parameters**
- `GenomeDISCO|subsampling` This allows subsampling the datasets to a specific desired sequencing depth. Possible values are: `lowest` (subsample to the depth of the sample with the lower sequencing depth from the pair being compared),
`<samplename>` where <samplename> is the name of the sample that is used to determine the sequencing depth to subsample from.

- `GenomeDISCO|tmin` The minimum number of steps of random walk to perform. Integer, > 0.

- `GenomeDISCO|tmax` The max number of steps of random walk to perform. Integer, > tmin.

- `GenomeDISCO|norm` The normalization to use on the data when running GenomeDISCO. Possible values include: `uniform` (no normalization), `sqrtvc`.

- `GenomeDISCO|scoresByStep` Whether to report the score at each t. By default (GenomeDISCO|scoresByStep no), only the final reproducibility score is returned.

- `GenomeDISCO|removeDiag` Whether to set the diagonal to entries in the contact map to 0. By default (GenomeDISCO|removeDiag yes), the diagonal entries are set to 0.

- `GenomeDISCO|transition` Whether to convert the normalized contact map to an appropriate transition matrix before running the random walks. By default (GenomeDISCO|transition yes) the normalized contact map is converted to a proper tr
ansition matrix, such that all rows sum to 1 exactly.

**HiCRep parameters**
- `HiCRep|h` The h parameter in HiCRep that determines the extent of 2D smoothing. See the HiCRep paper (http://genome.cshlp.org/content/early/2017/08/30/gr.220640.117) for details. Integer, >=0.

- `HiCRep|maxdist` The maximum distance to consider when computing the HiCRep score. Integer, should be a amultiple of the resolution of the data.

**HiC-Spector parameters**
- `HiC-Spector|n` The number of eigenvectors to use for HiC-Spector. Integer, > 0.

**QuASAR parameters**
- `QuASAR|rebinning` The rebinning distance. See the QuASAR paper (https://www.biorxiv.org/content/early/2017/10/17/204438) for details. Integer.

**Job submission parameters**
- `SGE|text` Text to append to the job submission for SGE. The default is "-l h_vmem=3G".

- `slurm|text` Text to append to the job submission for slurm. The default is "--mem 3G".

**Note about normalization**: At the moment, the different methods operate on different types of normalizations. For GenomeDISCO, the user can specify the desired normalization. For HiCRep and HiC-Spector the scores are computed on the
provided data, without normalization.
Thus, if you have normalized data, then you can provide that as an input, and set `GenomeDISCO|norm` to uniform. If you have raw data, then your HiCRep and HiC-Spector scores will be run on the raw data, and GenomeDISCO will be run on t
he normalization you specify with `GenomeDISCO|norm`.

###################################################################################################################################################40k
selects=("cc_rep1" "cc_rep2" "cc_rep3" "cc_rep4" "05h_rep2" "05h_rep3" "1h_rep1" "1h_rep2" "2h_rep1" "2h_rep2" "6h_rep5" "12h_rep1" "12h_rep2" "12h_rep3" "e2cell_rep3" "l2cell_rep1" "l2cell_rep2" "4cell_rep2" "4cell_rep3" "8cell_rep1" "8cell_rep2" "8cell_rep3" "morula_rep1" "morula_rep3" "morula_rep4" "icm_rep1" "icm_rep2" "te_rep1" "te_rep2" "te_rep3") 
selects=("e2cell_rep4" "6h_rep6" "kdm4d-e2cell_rep1" "kdm4d-e2cell_rep2" "partheno-1hpa_rep1" "partheno-1hpa_rep2" "partheno-e2cell_rep1" "partheno-e2cell_rep2" "sertoli-e2cell_rep1" "tsa-e2cell_rep1") #20190920
wdir=~/workspace/8.NT-HiC/o.HiCRep/4.3DRepQC/2.NT-40k/0.fithic
ddir=~/workspace/8.NT-HiC/5.maps/2.replicate/1.raw
bdir=~/workspace/8.NT-HiC/5.maps/2.replicate/3.biases
res=40000
for key in ${selects[@]}; do
(mkdir -p $wdir/$key/data $wdir/$key/output $wdir/$key/logs
nohup python $UtilsHP/hicpro2fithic.py -i $ddir/${key}_${res}.matrix -s $bdir/${key}_${res}_iced.matrix.biases \
-b $ConfigHP/${res}_mm10.bed -o $wdir/$key/data -r ${res} > $wdir/$key/logs/hicpro2fithic.log )&
done
rm 0.fithic/6h_rep3/ 0.fithic/e2cell_rep[12]

ddir=~/workspace/8.NT-HiC/o.HiCRep/4.3DRepQC/2.NT-40k/0.fithic/*/data/fithic.interactionCounts.gz
for i in ~/workspace/8.NT-HiC/o.HiCRep/4.3DRepQC/2.NT-40k/0.fithic/*/data/fithic.interactionCounts.gz; do key=${i##*0.fithic/};key=${key%%/data*};
echo -e $key"\t"$i >> metadata.samples
done

selects=("cc_rep1" "cc_rep2" "cc_rep3" "cc_rep4" "05h_rep2" "05h_rep3" "1h_rep1" "1h_rep2" "2h_rep1" "2h_rep2" "6h_rep5" "6h_rep6" "12h_rep1" "12h_rep2" "12h_rep3" "e2cell_rep3" "e2cell_rep4" "l2cell_rep1" "l2cell_rep2" "4cell_rep2" "4cell_rep3" "8cell_rep1" "8cell_rep2" "8cell_rep3" "morula_rep1" "morula_rep3" "morula_rep4" "icm_rep1" "icm_rep2" "te_rep1" "te_rep2" "te_rep3"  "kdm4d-e2cell_rep1" "kdm4d-e2cell_rep2" "partheno-1hpa_rep1" "partheno-1hpa_rep2" "partheno-e2cell_rep1" "partheno-e2cell_rep2" "sertoli-e2cell_rep1" "tsa-e2cell_rep1")

for i in ${selects[@]}; do echo $i;done | awk -v FS="_" 'BEGIN{k="";temp[0]} k!=$1{k=$1;delete temp;temp[$0];next} k==$1{for(i in temp){print i,$0};temp[$0]}' > metadata.pairs

#awk '{print $1,$2,$3,int(($2+$3)/2)}' $ConfigHP/40000_mm10.bed > $wdir/Bins_40000.bed
#gzip $wdir/Bins_40000.bed
wdir=~/workspace/8.NT-HiC/o.HiCRep/4.3DRepQC/2.NT-40k
for i in {1..6}; do 
nohup 3DChromatin_ReplicateQC run_all --metadata_samples $wdir/metadata.samples --metadata_pairs $wdir/metadata.pairs --bins $wdir/Bins_40000.bed.gz --outdir $wdir/$i.output --parameters_file $wdir/parameters/$i.parameters.txt > $wdir/logs/$i.log &
done


###########################################ES
wdir=~/workspace/8.NT-HiC/o.HiCRep/4.3DRepQC/4.ES/0.fithic
ddir=~/workspace/8.NT-HiC/e.ES/0.links/1.raw
bdir=~/workspace/8.NT-HiC/e.ES/0.links/6.biases
res=100000
mkd $wdir
for i in $ddir/xw_rep*${res}.matrix $ddir/es[125]00_rep*${res}.matrix; do key=$(basename $i _${res}.matrix);
(mkdir -p $wdir/$key/data $wdir/$key/output $wdir/$key/logs
nohup python $UtilsHP/hicpro2fithic.py -i $ddir/${key}_${res}.matrix -s $bdir/${key}_${res}_iced.matrix.biases \
-b $ConfigHP/${res}_mm10.bed -o $wdir/$key/data -r ${res} > $wdir/$key/logs/hicpro2fithic.log )&
done

wdir=~/workspace/8.NT-HiC/o.HiCRep/4.3DRepQC/4.ES
for i in $wdir/0.fithic/*/data/fithic.interactionCounts.gz; do key=${i##*0.fithic/};key=${key%%/data*};
echo -e $key"\t"$i >> $wdir/metadata.samples
done

for i in $wdir/0.fithic/*/data/fithic.interactionCounts.gz; do key=${i##*0.fithic/};key=${key%%/data*};echo $key;done | \
awk '{a[NR]=$1;b=NR} END{for(i=1;i<b;i++){for(j=i+1;j<=b;j++){print a[i],a[j]}}}' > $wdir/metadata.pairs

cp -r $wdir/../1.NT/parameters/ $wdir/
cp ../1.NT/Bins_100000.bed.gz ./
mkd $wdir/logs $wdir/pdfs $wdir/genome.pdfs
for i in {1..6}; do 
mkd $wdir/$i.output
nohup 3DChromatin_ReplicateQC run_all --metadata_samples $wdir/metadata.samples --metadata_pairs $wdir/metadata.pairs --bins $wdir/Bins_${res}.bed.gz \
--outdir $wdir/$i.output --parameters_file $wdir/parameters/$i.parameters.txt > $wdir/logs/$i.log &
done

awk 'FNR>1{print }' $wdir/*.output/scores/reproducibility.genomewide.txt > 
