##### Manuscript accepted in Genes and Development: "Gadd45 promotes locus specific DNA demethylation and 2C cycling in embryonic stem cells"


## Methylome data pre-processing and downstream analysis in Co mESCs and Gadd45-TKO
Collection of scripts to detect Differentially Methylated Regions and perform downstream enrichment analysis at several genomic features of mouse ES cells.

## Set up the files for Job Array and subsequent parallelization of the workflow
Before performing the quality control we want to split the fastq files in small pieces of 1mln reads in order to parallelize all the analysis:

`zcat example.fq.gz | split -4000000 &` 

we use 4mln lines because the information of every read within the fastq file takes 4 lines. 4 lines x 1mln reads= 4mln. Afterwards we create a folder and we put inside all the files for each sample and we create a list file with the name of all the files. This will be important for the job array in every step of the analysis. An example of list file can be found in the folder "list files job array". 

### Create variable to export the file names for each sample  
```
export fastq_file_G452  =`sed -n "$SLURM_ARRAY_TASK_ID"p G45.TKO1`  
export fastq_file_G453  =`sed -n "$SLURM_ARRAY_TASK_ID"p G45.TKO2`  
export fastq_file_mESC1 =`sed -n "$SLURM_ARRAY_TASK_ID"p Co.ESCs.1`  
export fastq_file_mESC2 =`sed -n "$SLURM_ARRAY_TASK_ID"p Co.ESCs.2`  
export OxFeatures       =`sed -n "$SLURM_ARRAY_TASK_ID"p list.files.heatMap`  
export Genomicfeatures  =`sed -n "$SLURM_ARRAY_TASK_ID"p list.files.frequency.plot`  
```

### Bash script template for Job array
```
#!/bin/sh

#SBATCH --job-name=Alignment     # job name
#SBATCH --array=1-494            # Number of Job (according to the number of splitted files in the exported variables)
#SBATCH --nodes=1                # nodes
#SBATCH -p andrade               # queue
#SBATCH -A jgu-cbdm
#SBATCH -c 1                     # cores
#SBATCH --mem=10000M             # memory
#SBATCH --time=72:00:00          # time
#SBATCH --error=Alignment.err    # error file name
#SBATCH --output=Alignment.out   # output file name
#SBATCH --mail-user=t.andreani@imb-mainz.de  # email
#SBATCH --mail-type=ALL                      # type notification
```
The bash commands above are used for all the steps of the pre processing and when required in the post processing (for heatMap and frquency plot) of the methylome data using Job Arrays in the Slurm queue system environment of Mogon.     

## Quality Control and Trim performed with trim_galore version 0.4.4: select reads with a quality value > 20 and remove illumina adapters    

#### Gadd45.TKO samples  
```
trim_galore --paired --trim1 $fastq_file_G452.R1.fastq.gz $fastq_file_G452.R2.fastq.gz    
trim_galore --paired --trim1 $fastq_file_G453.R1.fastq.gz $fastq_file_G453.R2.fastq.gz  
```
#### Co mESCs   

```
trim_galore --paired --trim1 $fastq_file_mESC1.R1.fastq.gz $fastq_file_mESC1.R2.gz   
trim_galore --paired --trim1 $fastq_file_mESC2.R1.fastq.gz $fastq_file_mESC2.R2.gz   
```
## Alignment, merging and extraction of CGH, CHG and CHH  with bismark_v0.18.0

### Alignment  
#### Gadd45.TKO samples  
```
bismark --bowtie2 -n 1 -I 0 -X 1000 --score_min L,0,-0.6 genome -1 $fastq_file_G452.R1.fastq.gz -2 $fastq_file_G452.R2.fastq.gz -o Alignment/   
bismark --bowtie2 -n 1 -I 0 -X 1000 --score_min L,0,-0.6 genome -1 $fastq_file_G453.R1.fastq.gz -2 $fastq_file_G453.R2.fastq.gz -o Alignment/   
```

#### Co mESCs 
```
bismark --bowtie2 -n 1 -I 0 -X 1000 --score_min L,0,-0.6 genome -1  $fastq_file_mESC1.R1.fastq.gz -2  $fastq_file_mESC1.R2.fastq.gz -o Alignment/   
bismark --bowtie2 -n 1 -I 0 -X 1000 --score_min L,0,-0.6 genome -1  $fastq_file_mESC2.R1.fastq.gz -2  $fastq_file_mESC2.R2.fastq.gz -o Alignment/   
```

### Merge the samples sorting for bismark by name always 

#### Gadd45.TKO samples  
```
samtools merge -n $gadd45Tko2/Alignment/Gadd45.tko2.bam $gadd45Tko2/Alignment/x*.bam    
samtools merge -n $gadd45Tko3/Alignment/Gadd45.tko3.bam $gadd45Tko3/Alignment/x*.bam  
```
#### Co mESCs  
```
samtools merge -n $mESC1/Alignment/mESC1.bam $mESC1/Alignment/x*.bam    
samtools merge -n $mESC2/Alignment/mESC2.bam $mESC2/Alignment/x*.bam    
```
## Extraction  
#### Gadd45.TKO samples    
```
bismark_methylation_extractor -p --ignore 5 --ignore_r2 5 --ample_memory --bedGraph --counts --cytosine_report --buffer_size 10G --genome_folder $genome/genome $extraction/Gadd45.tko2.bam  
bismark_methylation_extractor -p --ignore 5 --ignore_r2 5 --ample_memory --bedGraph --counts --cytosine_report --buffer_size 10G --genome_folder $genome/genome $extraction/Gadd45.tko3.bam  
```
#### Co mESCs     
```
bismark_methylation_extractor -p --ignore 5 --ignore_r2 5 --ample_memory --bedGraph --counts --cytosine_report --buffer_size 10G --genome_folder $genome/genome $extraction/mESC1.bam  
bismark_methylation_extractor -p --ignore 5 --ignore_r2 5 --ample_memory --bedGraph --counts --cytosine_report --buffer_size 10G --genome_folder $genome/genome $extraction/mESC2.bam  
```



## Call of Differentially Methylated Regions (DMRs) with MethylKit (contents of the script in the folder "R"):
`Rscript DMRs.Estimation.r`  

## Downstream post-processing analysis (enrichment at regulatory features from Table Browser UCSC)  

```
#Create Random set regions similar to the DMRs in terms of numbers (6904), length and CG composition to compute the expectation. This task is performed with bedtools shuffle:  

bedtools shuffle -incl background.file.with.all.DNA.regions.tested.bed -i Hyper.DMRs.G45.TKO.100bp.2CpG.Delta30.FDR.0.05.bed -g mm10.chr.sizes > RANDOM.bed

#The Intersection of the Hyper-DMRs in each feature is obtained:  

#Observed
for i in *UCSC.feature.bed;
do
  bedtools intersect -a Hyper.DMRs.G45.TKO.100bp.2CpG.Delta30.FDR.0.05.bed -b $i -wa | sort -u > Hyper.DMRs.in.$i; 
done   

#Expected
for i in *UCSC.feature.bed;
do
  bedtools intersect -a RANDOM.bed -b $i -wa | sort -u > RANDOM.regions.in.$i;
done  
```

## Motif Analysis with HOMER v3.12, 6-8-2012(fig. 3-A)
```
findMotifsGenome.pl Hyper.DMRs.G45.TKO.100bp.2CpG.Delta30.FDR.0.05.Promoter.bed mm10 motifDMR.Promoter -size 25 -len 8
findMotifsGenome.pl Hyper.DMRs.G45.TKO.100bp.2CpG.Delta30.FDR.0.05.Enhancer.bed mm10 motifDMR.Enhancer -size 25 -len 8
findMotifsGenome.pl Hyper.DMRs.G45.TKO.100bp.2CpG.Delta30.FDR.0.05.5hmC.bed mm10 motifDMR.5hmC -size 25 -len 8
findMotifsGenome.pl Hyper.DMRs.G45.TKO.100bp.2CpG.Delta30.FDR.0.05.5fC.bed mm10 motifDMR.5fC -size 25 -len 8
findMotifsGenome.pl Hyper.DMRs.G45.TKO.100bp.2CpG.Delta30.FDR.0.05.5caC.bed mm10 motifDMR.5caC -size 25 -len 8
findMotifsGenome.pl Hyper.DMRs.G45.TKO.100bp.2CpG.Delta30.FDR.0.05.Promoter.5hmC.bed mm10 motifDMR.Promoter.5hmC -size 25 -len 8
findMotifsGenome.pl Hyper.DMRs.G45.TKO.100bp.2CpG.Delta30.FDR.0.05.Promoter.5fC.bed mm10 motifDMR.Promoter.5fC -size 25 -len 8
findMotifsGenome.pl Hyper.DMRs.G45.TKO.100bp.2CpG.Delta30.FDR.0.05.Promoter.5caC.bed mm10 motifDMR.Promoter.5caC -size 25 -len 8
findMotifsGenome.pl Hyper.DMRs.G45.TKO.100bp.2CpG.Delta30.FDR.0.05.Enhancer.5hmC.bed mm10 motifDMR.Enhancer.5hmC -size 25 -len 8
findMotifsGenome.pl Hyper.DMRs.G45.TKO.100bp.2CpG.Delta30.FDR.0.05.Enhancer.5fC.bed mm10 motifDMR.Enhcaner.5fC -size 25 -len 8
findMotifsGenome.pl Hyper.DMRs.G45.TKO.100bp.2CpG.Delta30.FDR.0.05.Enhancer.5caC.bed mm10 motifDMR.Enhcnaer.5caC -size 25 -len 8
```

## Heatmap and frequency plot with deepTools 3.0.1  

#### For Heatmap of Hyper-DMRs at 5mC oxidative products in figure 2-D 
```
computeMatrix reference-point --referencePoint center -b 5000 -a 5000
 -R Hyper.DMRs.G45.TKO.100bp.2CpG.Delta30.FDR.0.05.bed
 -S $OxFeatures ## Oxidative feature from 'list.files.heatMap in the List.Files.Job.Array folder' 
 --missingDataAsZero
 -out Matrix.$OxFeatures.gz
 --outFileSortedRegions regions.$OxFeatures.bed

plotHeatmap
 -m Matrix.$OxFeatures.gz 
 -out Matrix.$OxFeatures.png --colorList --colorList cornflowerblue,yellow,red --missingDataColor white   
```
#### For frequency plot at several genomic features figure 2-F

#### First create the delta methylation differences in bins of 100 bp in TKO vs WT and the bigWig file

```bash ./Create.Delta.Mean.Methylation.file.Table.sh ```

#### Second create the matrix for the frequency plot

```
computeMatrix reference-point --referencePoint center -b 5000 -a 5000
 -R $Genomicfeatures ## feature from 'list.files.frequency.plot in the List.Files.Job.Array folder'
 -S G45.TKO.minus.ESCsWT.bin100bp.bigWig 
 --missingDataAsZero
 -out Matrix.$Genomicfeatures.gz
 --outFileSortedRegions regions.$Genomicfeatures.bed 
```

#### Third plot according to the matrix created: the plot represents the delta methylation differences signal in TKO vs WT at the center of each genomic feature

``` 
  plotHeatmap
 -m Matrix.$Genomicfeatures.gz 
 -out Matrix.$Genomicfeatures.png --colorList cornflowerblue,yellow,red --missingDataColor white --legendLocation none
 ```
  

