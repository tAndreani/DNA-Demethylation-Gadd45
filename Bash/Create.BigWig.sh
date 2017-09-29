#Compute the coverage
zcat Gadd45.tko1.bismark.cov.gz | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5+$6}' > Gadd45.tko1.bismark.cov.bed &
zcat Gadd45.tko2.fast.bismark.cov.gz | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5+$6}' > Gadd45.tko2.fast.bismark.cov.bed &
zcat Gadd45.tko3.fast.bismark.cov.gz | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5+$6}' > Gadd45.tko3.fast.bismark.cov.bed &
zcat mESC.bismark.cov.gz | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5+$6}' > mESC.bismark.cov.bed &


#Filter CG with Coverage >= 10
cat Gadd45.tko1.bismark.cov.bed | awk '$5 >= 10' > Gadd45.tko1.bismark.cov.major.egual.10.bed &
cat Gadd45.tko2.fast.bismark.cov.bed | awk '$5 >= 10' > Gadd45.tko2.fast.bismark.cov.major.egual.10.bed &
cat Gadd45.tko3.fast.bismark.cov.bed | awk '$5 >= 10' > Gadd45.tko3.fast.bismark.cov.major.egual.10.bed &
cat mESC.bismark.cov.bed | awk '$5 >= 10' > mESC.bismark.cov.major.egual.10.bed &


#Download wigToGiWig

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig


#Convert the files in wig format #chr, start, end, meth	values
for i in *egual.10.bed; 
do
        cut -f 1,2,3,4 $i > $i.wig; 
done &


#Convert to BigWig
./wigToBigWig Gadd45.tko1.bismark.cov.major.egual.10.bed.wig mm10.chr.sizes Gadd45.tko1.bismark.cov.major.egual.10.bed.bigWig &
./wigToBigWig Gadd45.tko2.fast.bismark.cov.major.egual.10.bed.wig mm10.chr.sizes Gadd45.tko2.fast.bismark.cov.major.egual.10.bed.bigWig &
./wigToBigWig Gadd45.tko3.fast.bismark.cov.major.egual.10.bed.wig mm10.chr.sizes Gadd45.tko3.fast.bismark.cov.major.egual.10.bed.bigWig &
./wigToBigWig mESC.bismark.cov.major.egual.10.bed.wig mm10.chr.sizes mESC.bismark.cov.major.egual.10.bed.bigWig