#!/bin/bash

#SBATCH --job-name=metaPlot 
#SBATCH --nodes=1                            # nb of nodes
#SBATCH -p short                             # queue -- reserves entire node 
#SBATCH -A imb
#SBATCH --cpus-per-task=2         # nb of cores/task
#SBATCH --ntasks=1    
#SBATCH --mem=24G        # large memory needed for sorting step (8G not enough )
#SBATCH --time=2:00:00           

# Computes metagene plot showing the DNA methylation changes around a genomic feature (here CEBPb peaks)
# Script by Tommaso Andreani
# with contribution by David Fournier
# to be launched using: sbatch Create.Delta.Mean.Methylation.file.Table.sh in SLURM cluster or bash ./Create.Delta.Mean.Methylation.file.Table.sh in your command line
# Two conditions: WT (wild-type) and Cond (condition)
# each have two replicates: WT replicates are named sample1 and sample2
# Cond replicates are names sample3 and sample4

# Module loading 

module load bedtools
module load datamash
module load deepTools

# Variables

datadir="/path/to/bismark_methyl_extractions" # directory where the CpG methylation calls are (files from Bismark with name format CpG_report.txt)
outdir="/path/to/plots" # output directory which will contain the metagene plot

###################################################
###Extract methylation values with coverage >= 10##
###then take the common CG shared among all the###
###replicates######################################

#From the methylation extracted files, creates a coverage file with methylation values
# +filters CG with Coverage >= 10

for sample in Sample1 Sample2 Sample3 Sample4; do # sorting step
    sort -V -k1,1 -k2,2 ${datadir}/${sample}/${sample}_bismark_bt2_pe.CpG_report.txt > ${sample}_sorted.txt
done
    
for sample in Sample1 Sample2 Sample3 Sample4; do
    awk '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print $1"\t"$2-1"\t"$2"\t"$4/($4+$5), $4+$5;}' ${sample}_sorted.txt > tmp.bed    
    cat tmp.bed | awk '$5 >= 10' > ${sample}_cov10.bed 
done


rm tmp.bed

#Create Id values for each CpG
for i in *cov10.bed;
do
    awk -F "\t" '{print $1"_"$2"_"$3}' $i > $i.Id
done 

# Sorting step necessary for command "comm" (otherwise comm won't work)
for i in *.Id; do 
sort $i -n > $i.sort
done

#Take common Id

comm -12 Sample3_cov10.bed.Id.sort  Sample4_cov10.bed.Id.sort | comm -12 - Sample1_cov10.bed.Id.sort | comm -12 -  Sample2_cov10.bed.Id.sort  > Common.CG.all.samples.txt

#Format the Coordinates and sort again
cat Common.CG.all.samples.txt | awk -F "_" '{print $1"\t"$2"\t"$3}' > Common.CG.all.samples.sort.bed
cat Common.CG.all.samples.sort.bed | sort -k 1,1 -k2,2n  > Common.CG.all.samples.sort.as.bed 

#Extract and quality check all are good with the same number
for i in *cov10.bed; do bedtools intersect -a $i -b Common.CG.all.samples.sort.as.bed > $i.common; done

#Create master table
paste *common | cut -f 1,2,3,4,9,14,19 | awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$7"\t"$4"\t"$5}' > Master.Table.wt1.wt2.cond1.cond2.txt 

#Compute Mean Values
cat Master.Table.wt1.wt2.cond1.cond2.txt | awk '{print $1"\t"$2"\t"$3"\t"($4+$5)/2"\t"($6+$7)/2}' > Master.Table.wt.cond.Mean.Values.txt  

#Compute delta
cat Master.Table.wt.cond.Mean.Values.txt | awk '{print $1"\t"$2"\t"$3"\t"$5-$4}' > Master.Table.wt.cond.Mean.Values.and.Delta.cond.minus.WT.txt

#Binning the mouse genome into 100 bp bins

wget https://raw.githubusercontent.com/tAndreani/DNA-Demethylation-Gadd45/master/List.Files.Job.Array/mm10.chr.size
sort -n  mm10.chr.size > tmp.bed
mv tmp.bed mm10.chr.size
bedtools makewindows -g mm10.chr.size -w 100 > mm10.binned.100.bp.ready.bed

##Extract delta mean in bins of 100 bp
sed -i 's/\./,/g' Master.Table.wt.cond.Mean.Values.and.Delta.cond.minus.WT.txt # datamash does not support . as floating point
bedtools intersect -a Master.Table.wt.cond.Mean.Values.and.Delta.cond.minus.WT.txt -b mm10.binned.100.bp.ready.bed -wb > mm10.binned.100.bp.ready.CpG.mean.values.2.cond.minus.2.WT
datamash -g 5,6,7 mean 4 < mm10.binned.100.bp.ready.CpG.mean.values.2.cond.minus.2.WT > Master.Table.wt.cond.Mean.Values.and.Delta.cond.minus.WT.bin100bp.txt

sortBed -i Master.Table.wt.cond.Mean.Values.and.Delta.cond.minus.WT.bin100bp.txt > Master.Table.wt.cond.Mean.Values.and.Delta.cond.minus.WT.bin100bp.bedGraph


# bigWig file generation:

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig -P ~/folders/software/
sed -i 's/,/\./g' Master.Table.wt.cond.Mean.Values.and.Delta.cond.minus.WT.bin100bp.bedGraph # bedGraphToBigWig does not support , as floating point
~/folders/software/bedGraphToBigWig Master.Table.wt.cond.Mean.Values.and.Delta.cond.minus.WT.bin100bp.bedGraph mm10.chr.size Master.Table.wt.cond.Mean.Values.and.Delta.cond.minus.WT.bin100bp.bw


# Cleaning

if false;then # toggle on to remove unecessary files
    rm *cov10.bed *.Id *Id.sort Common* *common mm10.binned.100.bp.ready.CpG.mean.values.2.cond.minus.2.WT Master.Table.wt.cond.Mean.Values.and.Delta.cond.minus.WT.bin100bp.txt
fi    

# metagene plot generation using deepTools

i="/path/to/ChIP-seq/data/CEBPB_peaks.bed"

j="Master.Table.wt.cond.Mean.Values.and.Delta.cond.minus.WT.bin100bp.bw"
computeMatrix reference-point --referencePoint center -S ${j} -R ${i} -a 5000 -b 5000 --samplesLabel DNAmet_CEBPBb -out ${outdir}/matrix_dt_CEBPb_vs_DNAmet_join.gz 

gunzip ${outdir}/matrix_dt_CEBPb_vs_DNAmet_join.gz 
sed -i 's/nan/0/g' ${outdir}/matrix_dt_CEBPb_vs_DNAmet_join
gzip ${outdir}/matrix_dt_CEBPb_vs_DNAmet_join
fi

plotHeatmap -m ${outdir}/matrix_dt_CEBPb_vs_DNAmet_join.gz -out ${outdir}/centre_dt_CEBPb_vs_DNAmet_join_5kb.png --colorList cornflowerblue,yellow,darkred -x CEBPb_peaks -y DNAmet_level --refPointLabel center --legendLocation none --boxArouncondatmaps no 


