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


#Download bedGraphToBigWig

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig


#Convert the files in wig format #chr, start, end, meth	values
for i in *egual.10.bed; 
do
        cut -f 1,2,3,4 $i > $i.bedgraph; 
done &

## Sort he files
LC_ALL=C sort -T Temporary -k1,1 -k2,2n -o Gadd45.tko1.bismark.cov.major.egual.10.sort.bedgraph Gadd45.tko1.bismark.cov.major.egual.10.bedgraph
LC_ALL=C sort -T Temporary -k1,1 -k2,2n -o Gadd45.tko2.fast.bismark.cov.major.egual.10.sort.bedgraph Gadd45.tko2.fast.bismark.cov.major.egual.10.bedgraph
LC_ALL=C sort -T Temporary -k1,1 -k2,2n -o Gadd45.tko3.fast.bismark.cov.major.egual.10.sort.bedgraph Gadd45.tko3.fast.bismark.cov.major.egual.10.bedgraph
LC_ALL=C sort -T Temporary -k1,1 -k2,2n -o mESC.bismark.cov.major.egual.10.sort.bedgraph mESC.bismark.cov.major.egual.10.bedgraph

#Correct coordinates
for i in *.cov.major.egual.10.sort.bedgraph;
do
       awk -F "\t" '{print $1"\t"$2-1"\t"$2"\t"$4}' $i > $i.coordinates;
done &


#Convert to BigWig
./bedGraphToBigWig Gadd45.tko1.bismark.cov.major.egual.10.bed.bedgraph.coordinates mm10.chr.sizes Gadd45.tko1.bismark.cov.major.egual.10.bed.bigWig &
./bedGraphToBigWig Gadd45.tko2.fast.bismark.cov.major.egual.10.bed.bedgraph.coordinates mm10.chr.sizes Gadd45.tko2.fast.bismark.cov.major.egual.10.bed.bigWig &
./bedGraphToBigWig Gadd45.tko3.fast.bismark.cov.major.egual.10.bed.bedgraph.coordinates mm10.chr.sizes Gadd45.tko3.fast.bismark.cov.major.egual.10.bed.bigWig &
./bedGraphToBigWig mESC.bismark.cov.major.egual.10.bed.bedgraph.coordinates mm10.chr.sizes mESC.bismark.cov.major.egual.10.bed.bigWig