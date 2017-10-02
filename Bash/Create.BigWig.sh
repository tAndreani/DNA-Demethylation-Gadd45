#From the methylation extracted files (all of them with 43735674 CG extracted ) create a coverage file
zcat ES.WT.myCpG.gz | awk '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print $1"\t"$2-1"\t"$2"\t"$4/($4+$5),
$4+$5;}' > ES.WT.myCpG.bed &
zcat G45.TKO1.myCpG.gz | awk '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print $1"\t"$2-1"\t"$2"\t"$3"\t"$4/($4+$5),
$4+$5;}' > G45.TKO1.myCpG.bed &
zcat G45.TKO2.myCpG.gz | awk '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print $1"\t"$2-1"\t"$2"\t"$3"\t"$4/($4+$5),
$4+$5;}' > G45.TKO2.myCpG.bed &
zcat G45.TKO3.myCpG.gz | awk '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print $1"\t"$2-1"\t"$2"\t"$3"\t"$4/($4+$5),
$4+$5;}' > G45.TKO3.myCpG.bed &

#Filter CG with Coverage >= 10
cat ES.WT.myCpG.bed | awk '$5 >= 10' > ES.WT.cov.major.egual.10.bed &
cat G45.TKO1.myCpG.bed | awk '$5 >= 10' > G45.TKO1.cov.major.egual.10.bed &
cat G45.TKO2.myCpG.bed | awk '$5 >= 10' > G45.TKO2.cov.major.egual.10.bed &
cat G45.TKO3.myCpG.bed | awk '$5 >= 10' > G45.TKO3.cov.major.egual.10.bed &


#Download bedGraphToBigWig

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig


#Convert the files in wig format #chr, start, end, meth	values
for i in *egual.10.bed; 
do
        cut -f 1,2,3,4 $i > $i.bedgraph; 
done &


#Convert to BigWig
./bedGraphToBigWig ES.WT.cov.major.egual.10.bed.bedgraph mm10.chr.sizes ES.WT.myCpG.cov.major.egual.10.bigWig &
./bedGraphToBigWig G45.TKO1.cov.major.egual.10.bed.bedgraph mm10.chr.sizes G45.TKO1.cov.major.egual.10.bigWig &
./bedGraphToBigWig G45.TKO2.cov.major.egual.10.bed.bedgraph mm10.chr.sizes G45.TKO2.cov.major.egual.10.bigWig &
./bedGraphToBigWig G45.TKO3.cov.major.egual.10.bed.bedgraph mm10.chr.sizes G45.TKO3.cov.major.egual.10.bigWig &



###Now for coverage >= 5

#Filter CG with Coverage >= 10
cat ES.WT.myCpG.bed | awk '$5 >= 5' > ES.WT.cov.major.egual.5.bed &
cat G45.TKO1.myCpG.bed | awk '$5 >= 5' > G45.TKO1.cov.major.egual.5.bed &
cat G45.TKO2.myCpG.bed | awk '$5 >= 5' > G45.TKO2.cov.major.egual.5.bed &
cat G45.TKO3.myCpG.bed | awk '$5 >= 5' > G45.TKO3.cov.major.egual.5.bed &


#Download bedGraphToBigWig

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig


#Convert the files in wig format #chr, start, end, meth	values
for i in *egual.5.bed; 
do
        cut -f 1,2,3,4 $i > $i.bedgraph; 
done &


#Convert to BigWig
./bedGraphToBigWig ES.WT.cov.major.egual.5.bed.bedgraph mm10.chr.sizes ES.WT.myCpG.cov.major.egual.5.bigWig &
./bedGraphToBigWig G45.TKO1.cov.major.egual.5.bed.bedgraph mm10.chr.sizes G45.TKO1.cov.major.egual.5.bigWig &
./bedGraphToBigWig G45.TKO2.cov.major.egual.5.bed.bedgraph mm10.chr.sizes G45.TKO2.cov.major.egual.5.bigWig &
./bedGraphToBigWig G45.TKO3.cov.major.egual.5.bed.bedgraph mm10.chr.sizes G45.TKO3.cov.major.egual.5.bigWig &
