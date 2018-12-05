###################################################
###Extract methylation values with coverage >= 10##
###then take the common CG shared among all the###
###replicates######################################



#From the methylation extracted files (all of them with 43735674 CG extracted ) create a coverage file with methylation values

zcat ES.WT1.myCpG.gz | awk '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print $1"\t"$2-1"\t"$2"\t"$4/($4+$5),
$4+$5;}' > ES.WT1.myCpG.bed &
zcat ES.WT2.myCpG.gz | awk '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print $1"\t"$2-1"\t"$2"\t"$4/($4+$5),
$4+$5;}' > ES.WT2.myCpG.bed &
zcat G45.TKO2.myCpG.gz | awk '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print $1"\t"$2-1"\t"$2"\t"$4/($4+$5),
$4+$5;}' > G45.TKO2.myCpG.bed &
zcat G45.TKO3.myCpG.gz | awk '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print $1"\t"$2-1"\t"$2"\t"$4/($4+$5),
$4+$5;}' > G45.TKO3.myCpG.bed &

#Filter CG with Coverage >= 10
cat ES.WT1.myCpG.bed | awk '$5 >= 10' > ES.WT1.cov.major.egual.10.bed &
cat ES.WT2.myCpG.bed | awk '$5 >= 10' > ES.WT2.cov.major.egual.10.bed &
cat G45.TKO2.myCpG.bed | awk '$5 >= 10' > G45.TKO2.cov.major.egual.10.bed &
cat G45.TKO3.myCpG.bed | awk '$5 >= 10' > G45.TKO3.cov.major.egual.10.bed &


#Extract the common CG with coverage 10 that are present in all the 5 samples

#First create the Id values
for i in *.cov.major.egual.10;
do awk -F "\t" '{print $1"_"$2"_"$3}' $i > $i.Id;
done &




#Take common Id
comm -12 G45.TKO2.cov.major.egual.10.Id  G45.TKO3.cov.major.egual.10.Id | comm -12 - ES.WT1.cov.major.egual.10.Id |
comm -12 -  ES.WT1.cov.major.egual.10..Id  > Common.CG.all.sampels.txt

#Format the Coordinates and sort again to be sure
cat Common.CG.all.sampels.sort.txt | awk -F "_" '{print $1"\t"$2"\t"$3}' > Common.CG.all.sampels.sort.bed &
cat Common.CG.all.sampels.sort.bed | sort -k 1,1 -k2,2n  > Common.CG.all.sampels.sort.as.bed &

##############################################
## awk -v OFS="\t" to avoid tab in every step#
##############################################


#Extract
for i in *cov.major.egual.10; do bedtools intersect -a $i -b Common.CG.all.sampels.sort.as.bed > $i.common; done &

#Create master table
paste *common | cut -f 1,2,3,4,9,14,19 | awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$7"\t"$4"\t"$5"}' > Master.Table.wt1.wt2.tko2.tko3.txt &

#Compute Mean Values
cat Master.Table.wt1.wt2.tko2.tko3.txt | awk '{print $1"\t"$2"\t"$3"\t"($4+$5)/2"\t"($6+$7)/2}' > Master.Table.wt1.wt2.tko2.tko3.Mean.Values.txt  &

#Compute delta
cat Master.Table.wt1.wt2.tko2.tko3.Mean.Values.txt | awk '{print $1"\t"$2"\t"$3"\t"$5-$4}' > Master.Table.wt1.wt2.tko2.tko3.Mean.Values.and.Delta.TKO.minus.WT.txt  &


##Extract delta mean in bins of 100 bp
bedtools intersect -a Master.Table.wt1.wt2.tko2.tko3.Mean.Values.and.Delta.TKO.minus.WT.bedgraph -b mm10.binned.100.bp.ready.bed -wb | awk  '{print $8"\t"$4}' > mm10.binned.100.bp.ready.CpG.mean.values.2.TKO.minus.2.WT
datamash-1.3/datamash -g 1 mean 2 <mm10.binned.100.bp.ready.CpG.mean.values.2.TKO.minus.2.WT  > Master.Table.wt1.wt2.tko2.tko3.Mean.Values.and.Delta.TKO.minus.WT.bin100bp.txt
cat Master.Table.wt1.wt2.tko2.tko3.Mean.Values.and.Delta.TKO.minus.WT.bin100bp.txt | awk -F "_" '{print $1"\t"$2"\t"$3"\t"$4}' > cat Master.Table.wt1.wt2.tko2.tko3.Mean.Values.and.Delta.TKO.minus.WT.bin100bp.bed

#Create bigWig with bedgraphTobigWig tool from UCSC https://www.encodeproject.org/software/bedgraphtobigwig/


