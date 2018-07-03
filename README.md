## Methylome Data Processing and Analysis in mESCs upon Gadd45a,b,g Proteins Knock Out
Collection of scripts to detect Differentially Methylated Regions and perform downstream enrichment analysis at several genomic features of mouse ES cells.

## Quality Control and Trim
Before performing the quality control we want to split the fastq files in small pieces of 1mln reads in order to parallelize all the analysis:

`zcat example.fq.gz | split -4000000 &` 

we use 4mln lines because the information of every read within the fastq file takes 4 lines. 4 lines x 1mln reads= 4mln. Afterwards we create a folder and we put inside all the files for each sample and we create a list file with the name of all the files. This will be important for the job array in every step of the analysis. An example of list file can be found in the folder "list files job array". After we can start to exclude the reads with a quality value lower than 20 and trim those with bad quality at 5' and 3' ends. Here


`export fastq_file_G452`=`sed -n "$SLURM_ARRAY_TASK_ID"p list.files.gadd45.2`  
`export fastq_file_mESC`=`sed -n "$SLURM_ARRAY_TASK_ID"p list.files.mESC`  



#Select reads with a quality value > 20 All the Reads    
`trim_galore --paired --trim1 $gadd45Tko2/$fastq_file_G45.R1.fastq.gz $gadd45Tko1/$fastq_file_G45.R2.fastq.gz`  
`trim_galore --paired --trim1 $mESC/$fastq_file_mESC.R1.fastq.gz $mESC/$fastq_file_mESC.R2.gz`  

## Alignment (parallelization) and extraction of CGH, CHG and CHH  

#Gadd45.TKO2 494 files  
`$bismark/bismark --bowtie2 -n 1 -I 0 -X 1000 --score_min L,0,-0.6 $genome/genome -1 $gadd45Tko2/R1/$fastq_file -2 $gadd45Tko2/R2/$fastq_file -o $gadd45Tko2/Alignment/`  

#mESC Wild Type 502 files  
`$bismark/bismark --bowtie2 -n 1 -I 0 -X 1000 --score_min L,0,-0.6 $genome/genome -1 $mESC/R1/$fastq_file -2 $mESC/R2/$fastq_file -o $mESC/Alignment/`  

## Call of Differentially Methylated Regions (DMRs) with MethylKit and Methpipe


## Downstream post-processing analysis (enrichment at regulatory features, Oxidative products DIP-seq peaks after Tet and Tdg knock out, Rloops etc..)

