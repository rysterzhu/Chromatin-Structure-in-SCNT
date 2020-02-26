wdir=~/workspace/8.NT-HiC/k.fithic/2.xw_100k_20181202/early_2cell/4.distal_100k_20181205
wdir=~/workspace/8.NT-HiC/k.fithic/2.xw_100k_20181202/late_2cell/3.distal
wdir=~/workspace/8.NT-HiC/k.fithic/2.xw_100k_20181202/ICM/3.distal
wdir=~/workspace/8.NT-HiC/k.fithic/2.xw_100k_20181202/ICM/4.NT-NF

mkdir -p $wdir/1.allPairs_analysis/1.distal_4M-40M $wdir/1.allPairs_analysis/2.distal_400k-1M $wdir/1.allPairs_analysis/3.distal_200k-2M
ln -s $wdir/temp* $wdir/1.allPairs_analysis/
ln -s $wdir/*all.Pairs $wdir/1.allPairs_analysis/

#Enhancer-Promoter
cd $wdir/1.allPairs_analysis/
awk 'FILENAME=="temp"{a[$0]} FILENAME=="temp3"{b[$0]} FILENAME=="NF.all.Pairs"{left=$1"\t"$2"\t"$3"\t"$7;right=$4"\t"$5"\t"$6"\t"$7;
if((left in a)&&(right in b)){print left,right > "NF.ep.pairs"};
if((left in b)&&(right in a)){print right,left > "NF.ep.pairs"}
if((left in a)&&(right in a)){print left,right > "NF.pp.pairs";print right,left > "NF.pp.pairs"};
if((left in a)&&!(right in b)&&!(right in a)){print left,right > "NF.np.pairs"};
if((right in a)&&!(left in b)&&!(left in a)){print right,left > "NF.np.pairs"}
if(!(left in a)&&!(right in b)&&!(left in b)&&!(right in a)){print left,right > "NF.null.pairs"}
}' temp temp3 NF.all.Pairs &
awk 'FILENAME=="tempNT"{a[$0]} FILENAME=="temp3NT"{b[$0]} FILENAME=="NT.all.Pairs"{left=$1"\t"$2"\t"$3"\t"$7;right=$4"\t"$5"\t"$6"\t"$7;
if((left in a)&&(right in b)){print left,right > "NT.ep.pairs"};
if((left in b)&&(right in a)){print right,left > "NT.ep.pairs"}
if((left in a)&&(right in a)){print left,right > "NT.pp.pairs";print right,left > "NT.pp.pairs"};
if((left in a)&&!(right in b)&&!(right in a)){print left,right > "NT.np.pairs"};
if((right in a)&&!(left in b)&&!(left in a)){print right,left > "NT.np.pairs"}
if(!(left in a)&&!(right in b)&&!(left in b)&&!(right in a)){print left,right > "NT.null.pairs"}
}' tempNT temp3NT NT.all.Pairs &
wait

for i in *airs; do 
awk '{a=($3+$2)/2-($6+$7)/2;if(a<0){a=-a}} a>4e6&&a<4e7{print}' $i > 1.distal_4M-40M/$i &
done
for i in *airs; do 
awk '{a=($3+$2)/2-($6+$7)/2;if(a<0){a=-a}} a>4e5&&a<1e6{print}' $i > 2.distal_400k-1M/$i &
done
for i in *airs; do 
awk '{a=($3+$2)/2-($6+$7)/2;if(a<0){a=-a}} a>2e5&&a<2e6{print}' $i > 3.distal_200k-2M/$i &
done

#确定pair在A还是B compartmen上
awk 'NR>1{if($4>0){print $1,$2,$3,"A"};if($4<0){print $1,$2,$3,"B"}}' ~/workspace/8.NT-HiC/3.public_data/f.XW_analysis/7.newHomer/5.PC1_100-400k/late_2cell.PC1.bedGraph > NF.compartment  
cut -f 1-4 NF.ep.pairs | intersectBed -a - -b NF.compartment -r -e -f 0.5 -wao > temp2
cut -f 5-8 NF.ep.pairs | intersectBed -a - -b NF.compartment -r -e -f 0.5 -wao | \
awk 'NR==FNR{a[$4]=$8} NR>FNR{print $4,a[$4],$8}' temp2 - | awk '$2=="."||$3=="."{a+=1;next} $2!=$3{b+=1;next} $2=="A"{c+=1;next} $2=="B"{d+=1} END{print a,b,c,d}'  
cut -f 1-4 NF.null.pairs | intersectBed -a - -b NF.compartment -r -e -f 0.5 -wao > temp1
cut -f 5-8 NF.null.pairs | intersectBed -a - -b NF.compartment -r -e -f 0.5 -wao | \
awk 'NR==FNR{a[$4]=$8} NR>FNR{print $4,a[$4],$8}' temp1 - | awk '$2=="."||$3=="."{a+=1;next} $2!=$3{b+=1;next} $2=="A"{c+=1;next} $2=="B"{d+=1} END{print a,b,c,d}'
cut -f 1-4 NF.np.pairs | intersectBed -a - -b NF.compartment -r -e -f 0.5 -wao > temp2
cut -f 5-8 NF.np.pairs | intersectBed -a - -b NF.compartment -r -e -f 0.5 -wao | \
awk 'NR==FNR{a[$4]=$8} NR>FNR{print $4,a[$4],$8}' temp2 - | awk '$2=="."||$3=="."{a+=1;next} $2!=$3{b+=1;next} $2=="A"{c+=1;next} $2=="B"{d+=1} END{print a,b,c,d}'  
cut -f 1-4 NF.pp.pairs | intersectBed -a - -b NF.compartment -r -e -f 0.5 -wao > temp2
cut -f 5-8 NF.pp.pairs | intersectBed -a - -b NF.compartment -r -e -f 0.5 -wao | \
awk 'NR==FNR{a[$4]=$8} NR>FNR{print $4,a[$4],$8}' temp2 - | awk '$2=="."||$3=="."{a+=1;next} $2!=$3{b+=1;next} $2=="A"{c+=1;next} $2=="B"{d+=1} END{print a,b,c,d}' 
#NT
awk 'NR>1{if($4>0){print $1,$2,$3,"A"};if($4<0){print $1,$2,$3,"B"}}' ~/workspace/8.NT-HiC/h.homer_ALL/7.except_newHomer_20180705/5.PC1_100-400k/l2cell.PC1.bedGraph > NT.compartment  
cut -f 1-4 NT.ep.pairs | intersectBed -a - -b NT.compartment -r -e -f 0.5 -wao > temp2
cut -f 5-8 NT.ep.pairs | intersectBed -a - -b NT.compartment -r -e -f 0.5 -wao | \
awk 'NR==FNR{a[$4]=$8} NR>FNR{print $4,a[$4],$8}' temp2 - | awk '$2=="."||$3=="."{a+=1;next} $2!=$3{b+=1;next} $2=="A"{c+=1;next} $2=="B"{d+=1} END{print a,b,c,d}' > test
cut -f 1-4 NT.null.pairs | intersectBed -a - -b NT.compartment -r -e -f 0.5 -wao > temp1
cut -f 5-8 NT.null.pairs | intersectBed -a - -b NT.compartment -r -e -f 0.5 -wao | \
awk 'NR==FNR{a[$4]=$8} NR>FNR{print $4,a[$4],$8}' temp1 - | awk '$2=="."||$3=="."{a+=1;next} $2!=$3{b+=1;next} $2=="A"{c+=1;next} $2=="B"{d+=1} END{print a,b,c,d}' >> test
cut -f 1-4 NT.np.pairs | intersectBed -a - -b NT.compartment -r -e -f 0.5 -wao > temp2
cut -f 5-8 NT.np.pairs | intersectBed -a - -b NT.compartment -r -e -f 0.5 -wao | \
awk 'NR==FNR{a[$4]=$8} NR>FNR{print $4,a[$4],$8}' temp2 - | awk '$2=="."||$3=="."{a+=1;next} $2!=$3{b+=1;next} $2=="A"{c+=1;next} $2=="B"{d+=1} END{print a,b,c,d}' >> test
cut -f 1-4 NT.pp.pairs | intersectBed -a - -b NT.compartment -r -e -f 0.5 -wao > temp2
cut -f 5-8 NT.pp.pairs | intersectBed -a - -b NT.compartment -r -e -f 0.5 -wao | \
awk 'NR==FNR{a[$4]=$8} NR>FNR{print $4,a[$4],$8}' temp2 - | awk '$2=="."||$3=="."{a+=1;next} $2!=$3{b+=1;next} $2=="A"{c+=1;next} $2=="B"{d+=1} END{print a,b,c,d}' >> test


#EP与NP对应的基因的表达量
intersectBed -a NF.ep.pairs -b ~/ann/mm10_promoter_1000-1000.bed -r -e -f 0.5 -wo > NF.ep.genes
intersectBed -a NF.np.pairs -b ~/ann/mm10_promoter_1000-1000.bed -r -e -f 0.5 -wo > NF.np.genes
intersectBed -a NF.pp.pairs -b ~/ann/mm10_promoter_1000-1000.bed -r -e -f 0.5 -wo > NF.pp.genes

intersectBed -a NT.ep.pairs -b ~/ann/mm10_promoter_1000-1000.bed -r -e -f 0.5 -wo > NT.ep.genes
intersectBed -a NT.np.pairs -b ~/ann/mm10_promoter_1000-1000.bed -r -e -f 0.5 -wo > NT.np.genes
intersectBed -a NT.pp.pairs -b ~/ann/mm10_promoter_1000-1000.bed -r -e -f 0.5 -wo > NT.pp.genes
#
awk 'FILENAME=="NF.ep.genes"{a[$13]}
FILENAME=="NF.pp.genes"{b[$13]}
FILENAME=="NF.np.genes"{c[$13]}
END{for(i in a){print i,"EP"};for(i in b){if(!(i in a)){print i,"PP"}};for(i in c){if(!(i in a)&&!(i in b)){print i,"NP"}}}' NF.ep.genes NF.pp.genes NF.np.genes > NF.genes
awk 'FILENAME=="NT.ep.genes"{a[$13]}
FILENAME=="NT.pp.genes"{b[$13]}
FILENAME=="NT.np.genes"{c[$13]}
END{for(i in a){print i,"EP"};for(i in b){if(!(i in a)){print i,"PP"}};for(i in c){if(!(i in a)&&!(i in b)){print i,"NP"}}}' NT.ep.genes NT.pp.genes NT.np.genes > NT.genes

awk 'NR==FNR{a[$13]} NR>FNR{b[$13]} END{for(i in b){if(!(i in a)){print i,"PP"}};for(i in a){if(!(i in b)){print i,"EP"}else{print i,"both"}}}' NF.ep.genes NF.pp.genes > NF2.genes
awk 'NR==FNR{a[$13]} NR>FNR{b[$13]} END{for(i in b){if(!(i in a)){print i,"PP"}};for(i in a){if(!(i in b)){print i,"EP"}else{print i,"both"}}}' NT.ep.genes NT.pp.genes > NT2.genes


for k in NF.ep NT.ep NF.pp NT.pp;do
awk '{print $1,$2,$3;print $5,$6,$7}' $k.pairs | nohup findMotifsGenome.pl - mm10 $k.motif -size 500 -mask -p 8 &
done
awk 'NR==FNR{a[$13]} NR>FNR{b[$13]} END{for(i in b){if(i in a){print i,"both"}else{print i,"NT"}}for(i in a){if(!(i in b)){print i,"NF"}}}' NF.ep.genes NT.ep.genes > all.ep.genes