## Methylome Data Processing and Analysis in mESCs upon Gadd45a,b,g Proteins Knock Out
Collection of scripts to detect Differentially Methylated Regions and perform downstream enrichment analysis at several genomic features of mouse ES cells.

## Quality Control and Trim
Before performing the quality control we want to split the fastq files in small pieces of 1mln reads in order to parallelize all the analysis:

`zcat example.fq.gz | split -4000000 &` 

we use 4mln lines because the information of every read within the fastq file takes 4 lines. 4 lines x 1mln reads= 4mln. Afterwards we create a folder and we put inside all the files for each sample and we create a list file with the name of all the files. This will be important for the job array in every step of the analysis. An example of list file can be found in the folder "list files job array". After we can start to exclude the reads with a quality value lower than 20:

#Create variable to export the files' names for each sample  
`export fastq_file_G452`=`sed -n "$SLURM_ARRAY_TASK_ID"p list.files.gadd45.2`  
`export fastq_file_G453`=`sed -n "$SLURM_ARRAY_TASK_ID"p list.files.gadd45.3`  

`export fastq_file_mESC1`=`sed -n "$SLURM_ARRAY_TASK_ID"p list.files.mESC.1`  
`export fastq_file_mESC2`=`sed -n "$SLURM_ARRAY_TASK_ID"p list.files.mESC.2`  


#Select reads with a quality value > 20 All the Reads    
`trim_galore --paired --trim1 $fastq_file_G452.R1.fastq.gz $fastq_file_G452.R2.fastq.gz`  
`trim_galore --paired --trim1 $fastq_file_G453.R1.fastq.gz $fastq_file_G453.R2.fastq.gz`  

`trim_galore --paired --trim1 $fastq_file_mESC1.R1.fastq.gz $fastq_file_mESC1.R2.gz`  
`trim_galore --paired --trim1 $fastq_file_mESC2.R1.fastq.gz $fastq_file_mESC2.R2.gz`  

## Alignment (parallelization), merging and extraction of CGH, CHG and CHH  

#Alignment  
#Gadd45.TKO samples  
`bismark --bowtie2 -n 1 -I 0 -X 1000 --score_min L,0,-0.6 genome -1 $fastq_file_G452.R1.fastq.gz -2 $fastq_file_G452.R2.fastq.gz -o Alignment/`  
`bismark --bowtie2 -n 1 -I 0 -X 1000 --score_min L,0,-0.6 genome -1 $fastq_file_G453.R1.fastq.gz -2 $fastq_file_G453.R2.fastq.gz -o Alignment/`  

#mESC Wild Type samples  
`bismark --bowtie2 -n 1 -I 0 -X 1000 --score_min L,0,-0.6 genome -1  $fastq_file_mESC1.R1.fastq.gz -2  $fastq_file_mESC1.R2.fastq.gz -o Alignment/`  
`bismark --bowtie2 -n 1 -I 0 -X 1000 --score_min L,0,-0.6 genome -1  $fastq_file_mESC2.R1.fastq.gz -2  $fastq_file_mESC2.R2.fastq.gz -o Alignment/`  

#Merge the samples  
#sort for bismark by name always  
`samtools merge -n $gadd45Tko2/Alignment/Gadd45.tko2.bam $gadd45Tko2/Alignment/x*.bam  `  
`samtools merge -n $gadd45Tko3/Alignment/Gadd45.tko3.bam $gadd45Tko3/Alignment/x*.bam  `  
`samtools merge -n $mESC1/Alignment/mESC1.bam $mESC1/Alignment/x*.bam  `  
`samtools merge -n $mESC2/Alignment/mESC2.bam $mESC2/Alignment/x*.bam`  

#Extraction  
#Gadd45.TKO samples    
`$bismark/bismark_methylation_extractor -p --ignore 5 --ignore_r2 5 --ample_memory --bedGraph --counts --cytosine_report --buffer_size 10G --genome_folder $genome/genome $extraction/Gadd45.tko2.bam`
`$bismark/bismark_methylation_extractor -p --ignore 5 --ignore_r2 5 --ample_memory --bedGraph --counts --cytosine_report --buffer_size 10G --genome_folder $genome/genome $extraction/Gadd45.tko3.bam`

#mESC Wild Type samples    
`$bismark/bismark_methylation_extractor -p --ignore 5 --ignore_r2 5 --ample_memory --bedGraph --counts --cytosine_report --buffer_size 10G --genome_folder $genome/genome $extraction/Gadd45.tko3.bam`
`$bismark/bismark_methylation_extractor -p --ignore 5 --ignore_r2 5 --ample_memory --bedGraph --counts --cytosine_report --buffer_size 10G --genome_folder $genome/genome $extraction/mESC.bam`




## Call of Differentially Methylated Regions (DMRs) with MethylKit from the bash run:
`Rscript MethylKit.Analysis.R`  

## Downstream post-processing analysis (enrichment at regulatory features, Oxidative products DIP-seq peaks after Tet and Tdg knock out, Rloops etc..)

