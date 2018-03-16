control2 is ES wt from Neil experiment

###################################################
###Extract methylation values with coverage >= 10##
###then take the common CG shared among all the###
###replicates######################################



#From the methylation extracted files (all of them with 43735674 CG extracted ) create a coverage file with methylation values

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


#Extract the common CG with coverage 10 that are present in all the 5 samples

#First create the Id values
for i in *.cov.major.egual.10;
do awk -F "\t" '{print $1"_"$2"_"$3}' $i > $i.Id;
done &




#Take common Id
comm -12 Gadd45.tko1.CpG_report.sort.bed.cov.major.egual.10.Id  Gadd45.tko3.CpG_report.sort.bed.cov.major.egual.10.Id | comm -12 - mESC_CpG_report.fromNeil.sort.bed.cov.major.egual.10.Id |
 comm -12 - Gadd45.tko2.CpG_report.sort.bed.cov.major.egual.10.Id | comm -12 -  mESC.CpG_report.fromGadd45.sort.bed.cov.major.egual.10.Id  > Common.CG.all.sampels.txt

#Format the Coordinates and sort
cat Common.CG.all.sampels.sort.txt | awk -F "_" '{print $1"\t"$2"\t"$3}' > Common.CG.all.sampels.sort.bed &
cat Common.CG.all.sampels.sort.bed | sort -k 1,1 -k2,2n  > Common.CG.all.sampels.sort.as.bed &


#Extract
for i in *cov.major.egual.10; do bedtools intersect -a $i -b Common.CG.all.sampels.sort.as.bed > $i.common; done &

#Create master table
 paste * | cut -f 1,2,3,4,9,14,19,24,29 | awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$7"\t"$4"\t"$5"\t"$6}' > Master.Table.wtNeil.wtGadd.tko1.tko2.tko3.txt &

#Compute Mean Values
cat Master.Table.wtNeil.wtGadd.tko1.tko2.tko3.txt | awk '{print $1"\t"$2"\t"$3"\t"($4+$5)/2"\t"($6+$7+$8)/3}' > Master.Table.wtNeil.wtGadd.tko1.tko2.tko3.Mean.Values.txt  &

#Compute delta
cat Master.Table.wtNeil.wtGadd.tko1.tko2.tko3.Mean.Values.txt | awk '{print $1"\t"$2"\t"$3"\t"$5-$4}' > Master.Table.wtNeil.wtGadd.tko1.tko2.tko3.Mean.Values.and.Delta.TKO.minus.WT.txt  &

#Convert to bigWig
