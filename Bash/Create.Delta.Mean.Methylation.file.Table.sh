#!/bin/bash

#SBATCH --job-name=metaPlot 
#SBATCH --nodes=1                            # nb of nodes
#SBATCH -p short                             # queue -- reserves entire node 
#SBATCH -A imb
#SBATCH --cpus-per-task=2         # nb of cores/task
##SBATCH --ntasks=1     # variable set in main.sh script
#SBATCH --mem=24G        # large memory needed for sorting step (8G not enough )
#SBATCH --time=2:00:00                      # time  hh:mm:ss  

# Computes the DNA methylation landscape at a CEBP/hyperDMR site.
# Script by Tommaso Andreani
# adapted by David Fournier January 2019
# original script by TA to be found at the following URL:
# https://raw.githubusercontent.com/tAndreani/DNA-Demethylation-Gadd45/master/Bash/Create.Delta.Mean.Methylation.file.Table.sh
# to be launched using: sbatch metagene_plot.sh


# Module loading - IMB cluster version

module load bedtools
module load datamash

# Variables

wdir="/home/dafourni/folders/Lmna_aging/DMR"  # DMR analysis folder  -- IMB cluster
#dir="/home/dafourni/Desktop/Epigen/projects/Lmna_aging/WGBS_analysis/DMR"  # DMR analysis folder -- local computer
datadir="/fsimb/groups/imb-niehrsgr/imb_niehrs_2018_DKFZ_XI101_Lmna_WGBS/bismark_methyl_extractions"


###################################################
###Extract methylation values with coverage >= 10##
###then take the common CG shared among all the###
###replicates######################################

#From the methylation extracted files (all of them with 43735674 CG extracted ) create a coverage file with methylation values
# +filters CG with Coverage >= 10

for sample in AS-265808-LR-38379 AS-265809-LR-38380 AS-265810-LR-38381 AS-265811-LR-38382; do # sorting step (maybe not needed)
#    sortBed -i $i > ${i}_sorted 
    sort -V -k1,1 -k2,2 ${datadir}/${sample}/${sample}_bismark_bt2_pe.CpG_report.txt > ${sample}_sorted.txt
    #sort -k 1,1 -k2,2n ${datadir}/${sample}/${sample}_bismark_bt2_pe.CpG_report.txt > ${sample}_sorted.txt
done
    
for sample in AS-265808-LR-38379 AS-265809-LR-38380 AS-265810-LR-38381 AS-265811-LR-38382; do
    awk '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print $1"\t"$2-1"\t"$2"\t"$4/($4+$5), $4+$5;}' ${sample}_sorted.txt > tmp.bed    
    cat tmp.bed | awk '$5 >= 10' > ${sample}_cov10.bed 
done

rm tmp.bed

#Create Id values for each CpG
for i in *cov10.bed;
do
    awk -F "\t" '{print $1"_"$2"_"$3}' $i > $i.Id
done 

# DF: new sorting step necessary for command "comm" (otherwise comm won't work)
for i in *.Id; do 
sort $i -n > $i.sort
done


#Take common Id

comm -12 AS-265810-LR-38381_cov10.bed.Id.sort  AS-265811-LR-38382_cov10.bed.Id.sort | comm -12 - AS-265808-LR-38379_cov10.bed.Id.sort | comm -12 -  AS-265809-LR-38380_cov10.bed.Id.sort  > Common.CG.all.samples.txt

#Format the Coordinates and sort again
cat Common.CG.all.samples.txt | awk -F "_" '{print $1"\t"$2"\t"$3}' > Common.CG.all.samples.sort.bed
cat Common.CG.all.samples.sort.bed | sort -k 1,1 -k2,2n  > Common.CG.all.samples.sort.as.bed 

##############################################
## awk -v OFS="\t" to avoid tab in every step#
##############################################
    
#Extract and quality check all are good with the same number
for i in *cov10.bed; do bedtools intersect -a $i -b Common.CG.all.samples.sort.as.bed > $i.common; done

#Create master table
paste *common | cut -f 1,2,3,4,9,14,19 | awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$7"\t"$4"\t"$5}' > Master.Table.wt1.wt2.dhe1.dhe2.txt 

#Compute Mean Values
cat Master.Table.wt1.wt2.dhe1.dhe2.txt | awk '{print $1"\t"$2"\t"$3"\t"($4+$5)/2"\t"($6+$7)/2}' > Master.Table.wt.dhe.Mean.Values.txt  

#Compute delta
cat Master.Table.wt.dhe.Mean.Values.txt | awk '{print $1"\t"$2"\t"$3"\t"$5-$4}' > Master.Table.wt.dhe.Mean.Values.and.Delta.dhe.minus.WT.txt

#Binning the mouse genome into 100 bp bins

wget https://raw.githubusercontent.com/tAndreani/DNA-Demethylation-Gadd45/master/List.Files.Job.Array/mm10.chr.size
sort -n  mm10.chr.size > tmp.bed
mv tmp.bed mm10.chr.size
bedtools makewindows -g mm10.chr.size -w 100 > mm10.binned.100.bp.ready.bed

##Extract delta mean in bins of 100 bp
sed -i 's/\./,/g' Master.Table.wt.dhe.Mean.Values.and.Delta.dhe.minus.WT.txt # datamash does not support . as floating point
bedtools intersect -a Master.Table.wt.dhe.Mean.Values.and.Delta.dhe.minus.WT.txt -b mm10.binned.100.bp.ready.bed -wb > mm10.binned.100.bp.ready.CpG.mean.values.2.dhe.minus.2.WT
datamash -g 5,6,7 mean 4 < mm10.binned.100.bp.ready.CpG.mean.values.2.dhe.minus.2.WT > Master.Table.wt.dhe.Mean.Values.and.Delta.dhe.minus.WT.bin100bp.txt

sortBed -i Master.Table.wt.dhe.Mean.Values.and.Delta.dhe.minus.WT.bin100bp.txt > Master.Table.wt.dhe.Mean.Values.and.Delta.dhe.minus.WT.bin100bp.bedGraph


# bigWig file generation:

#wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig -P ~/folders/software/
sed -i 's/,/\./g' Master.Table.wt.dhe.Mean.Values.and.Delta.dhe.minus.WT.bin100bp.bedGraph # bedGraphToBigWig does not support , as floating point
~/folders/software/bedGraphToBigWig Master.Table.wt.dhe.Mean.Values.and.Delta.dhe.minus.WT.bin100bp.bedGraph mm10.chr.size Master.Table.wt.dhe.Mean.Values.and.Delta.dhe.minus.WT.bin100bp.bw
# This bigWig file will be used by the metagene plot tool (deepTools f.i.)

# Cleaning

if false;then
    rm *cov10.bed *.Id *Id.sort Common* *common mm10.binned.100.bp.ready.CpG.mean.values.2.dhe.minus.2.WT Master.Table.wt.dhe.Mean.Values.and.Delta.dhe.minus.WT.bin100bp.txt
fi    



