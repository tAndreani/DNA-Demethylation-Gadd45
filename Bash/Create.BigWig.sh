#############
#G45 Project#
#############

#From the methylation extracted files (all of them with 43735674 CG extracted ) create a coverage file
zcat ES.WT.myCpG.gz | awk '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print $1"\t"$2-1"\t"$2"\t"$4/($4+$5),
$4+$5;}' > ES.WT.myCpG.bed &
zcat G45.TKO1.myCpG.gz | awk '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print $1"\t"$2-1"\t"$2"\t"$4/($4+$5),
$4+$5;}' > G45.TKO1.myCpG.bed &
zcat G45.TKO2.myCpG.gz | awk '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print $1"\t"$2-1"\t"$2"\t"$4/($4+$5),
$4+$5;}' > G45.TKO2.myCpG.bed &
zcat G45.TKO3.myCpG.gz | awk '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print $1"\t"$2-1"\t"$2"\t"$4/($4+$5),
$4+$5;}' > G45.TKO3.myCpG.bed &

#Filter CG with Coverage >= 10
cat ES.WT.myCpG.bed | awk '$5 >= 10' > ES.WT.cov.major.egual.10.bed &
cat G45.TKO1.myCpG.bed | awk '$5 >= 10' > G45.TKO1.cov.major.egual.10.bed &
cat G45.TKO2.myCpG.bed | awk '$5 >= 10' > G45.TKO2.cov.major.egual.10.bed &
cat G45.TKO3.myCpG.bed | awk '$5 >= 10' > G45.TKO3.cov.major.egual.10.bed &


#Download bedGraphToBigWig

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig


#Convert the files in bedgraph format #chr, start, end, meth	values
for i in *egual.10.bed; 
do
        cut -f 1,2,3,4 $i > $i.bedgraph; 
done &


#Convert to BigWig
./bedGraphToBigWig ES.WT.cov.major.egual.10.bed.bedgraph mm10.chr.sizes ES.WT.cov.major.egual.10.bigWig &
./bedGraphToBigWig G45.TKO1.cov.major.egual.10.bed.bedgraph mm10.chr.sizes G45.TKO1.cov.major.egual.10.bigWig &
./bedGraphToBigWig G45.TKO2.cov.major.egual.10.bed.bedgraph mm10.chr.sizes G45.TKO2.cov.major.egual.10.bigWig &
./bedGraphToBigWig G45.TKO3.cov.major.egual.10.bed.bedgraph mm10.chr.sizes G45.TKO3.cov.major.egual.10.bigWig &



###Now for coverage >= 5

#Filter CG with Coverage >= 5
cat ES.WT.myCpG.bed | awk '$5 >= 5' > ES.WT.cov.major.egual.5.bed &
cat G45.TKO1.myCpG.bed | awk '$5 >= 5' > G45.TKO1.cov.major.egual.5.bed &
cat G45.TKO2.myCpG.bed | awk '$5 >= 5' > G45.TKO2.cov.major.egual.5.bed &
cat G45.TKO3.myCpG.bed | awk '$5 >= 5' > G45.TKO3.cov.major.egual.5.bed &


#Download bedGraphToBigWig

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig


#Convert the files in bedgraph format #chr, start, end, meth	values
for i in *egual.5.bed; 
do
        cut -f 1,2,3,4 $i > $i.bedgraph; 
done &


#Convert to BigWig
./bedGraphToBigWig ES.WT.cov.major.egual.5.bed.bedgraph mm10.chr.sizes ES.WT.cov.major.egual.5.bigWig &
./bedGraphToBigWig G45.TKO1.cov.major.egual.5.bed.bedgraph mm10.chr.sizes G45.TKO1.cov.major.egual.5.bigWig &
./bedGraphToBigWig G45.TKO2.cov.major.egual.5.bed.bedgraph mm10.chr.sizes G45.TKO2.cov.major.egual.5.bigWig &
./bedGraphToBigWig G45.TKO3.cov.major.egual.5.bed.bedgraph mm10.chr.sizes G45.TKO3.cov.major.egual.5.bigWig &


###########
##Neil DKO#
###########


#From the methylation extracted files (all of them with 43735674 CG extracted ) create a coverage file
cat ES.WT.myCpG.txt | awk '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print $1"\t"$2-1"\t"$2"\t"$4/($4+$5),
$4+$5;}' > ES.WT.myCpG.bed &
cat Neil.DKO.6.myCpG.txt | awk '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print $1"\t"$2-1"\t"$2"\t"$4/($4+$5),
$4+$5;}' > Neil.DKO.6.myCpG.bed &
cat Neil.DKO.7.myCpG.txt | awk '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print $1"\t"$2-1"\t"$2"\t"$4/($4+$5),
$4+$5;}' > Neil.DKO.7.myCpG.bed &
cat Neil.DKO.8.myCpG.txt | awk '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print $1"\t"$2-1"\t"$2"\t"$4/($4+$5),
$4+$5;}' > Neil.DKO.8.myCpG.bed &

#Filter CG with Coverage >= 10
cat ES.WT.myCpG.bed | awk '$5 >= 10' > ES.WT.cov.major.egual.10.bed &
cat Neil.DKO.6.myCpG.bed | awk '$5 >= 10' > Neil.DKO.6.cov.major.egual.10.bed &
cat Neil.DKO.7.myCpG.bed | awk '$5 >= 10' > Neil.DKO.7.cov.major.egual.10.bed &
cat Neil.DKO.8.myCpG.bed | awk '$5 >= 10' > Neil.DKO.8.cov.major.egual.10.bed &


#Download bedGraphToBigWig

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig


#Convert the files in bedgraph format #chr, start, end, meth	values
for i in *egual.10.bed; 
do
        cut -f 1,2,3,4 $i > $i.bedgraph; 
done &


#Convert to BigWig
./bedGraphToBigWig ES.WT.cov.major.egual.10.bed.bedgraph mm10.chr.sizes ES.WT.cov.major.egual.10.bigWig &
./bedGraphToBigWig Neil.DKO.6.cov.major.egual.10.bed.bedgraph mm10.chr.sizes Neil.DKO.6.cov.major.egual.10.bigWig &
./bedGraphToBigWig Neil.DKO.7.cov.major.egual.10.bed.bedgraph mm10.chr.sizes Neil.DKO.7.cov.major.egual.10.bigWig &
./bedGraphToBigWig Neil.DKO.8.cov.major.egual.10.bed.bedgraph mm10.chr.sizes Neil.DKO.8.cov.major.egual.10.bigWig &



###Now for coverage >= 5

#Filter CG with Coverage >= 5
cat ES.WT.myCpG.bed | awk '$5 >= 5' > ES.WT.cov.major.egual.5.bed &
cat Neil.DKO.6.myCpG.bed | awk '$5 >= 5' > Neil.DKO.6.cov.major.egual.5.bed &
cat Neil.DKO.7.myCpG.bed | awk '$5 >= 5' > Neil.DKO.7.cov.major.egual.5.bed &
cat Neil.DKO.8.myCpG.bed | awk '$5 >= 5' > Neil.DKO.8.cov.major.egual.5.bed &


#Download bedGraphToBigWig

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig


#Convert the files in bedgraph format #chr, start, end, meth	values
for i in *egual.5.bed; 
do
        cut -f 1,2,3,4 $i > $i.bedgraph; 
done &


#Convert to BigWig
./bedGraphToBigWig ES.WT.cov.major.egual.5.bed.bedgraph mm10.chr.sizes ES.WT.cov.major.egual.5.bigWig &
./bedGraphToBigWig Neil.DKO.6.cov.major.egual.5.bed.bedgraph mm10.chr.sizes Neil.DKO.6.cov.major.egual.5.bigWig &
./bedGraphToBigWig Neil.DKO.7.cov.major.egual.5.bed.bedgraph mm10.chr.sizes Neil.DKO.7.cov.major.egual.5.bigWig &
./bedGraphToBigWig Neil.DKO.8.cov.major.egual.5.bed.bedgraph mm10.chr.sizes Neil.DKO.8.cov.major.egual.5.bigWig &
