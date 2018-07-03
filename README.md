## Methylome Data Processing and Analysis in mESCs upon Gadd45a,b,g Proteins Knock Out
Collection of scripts to detect Differentially Methylated Regions and perform downstream enrichment analysis at several genomic features of mouse ES cells.

## Quality Control and Trim
Before performing the quality control we want to split the fastq files in small pieces of 1mln reads in order to parallelize all the analysis:

`zcat example.fq.gz | split -4000000 &` 

we use 4mln lines because the information of every read within the fastq file takes 4 lines. 4 lines x 1mln reads= 4mln. Afterwards we create a folder and we put inside all the files for each sample and we create a list file with the name of all the files. This will be important for the job array in every step of the analysis. An example of list file can be found in the folder "list files job array". After we can start to exclude the reads with a quality value lower than xxx and trim those with bad quality at 5' and 3' ends. Here


`export fastq_file=`sed -n "$SLURM_ARRAY_TASK_ID"p list.files.gadd45.1
`export fastq_file=`sed -n "$SLURM_ARRAY_TASK_ID"p list.files.mESC``



#Select reads with a quality value > 20 All the Reads
../../trim_galore --paired --trim1 $gadd45Tko1/AS-180617-LR-27523_R1.fastq.gz $gadd45Tko1/AS-180617-LR-27523_R2.fastq.gz
../../trim_galore --paired --trim1 $gadd45Tko2/AS-180618-LR-27524_R1.fastq.gz $gadd45Tko2/AS-180618-LR-27524_R2.fastq.gz
../../trim_galore --paired --trim1 $gadd45Tko3/AS-180619-LR-27525_R1.fastq.gz $gadd45Tko3/AS-180619-LR-27525_R2.fastq.gz
../../trim_galore --paired --trim1 $mESC/AS-180620-LR-28032_R1.fastq.gz $mESC/AS-180620-LR-28032_R2.fastq.gz

## Alignment (parallelization) and extraction of CGH, CHG and CHH


## Call of Differentially Methylated Regions (DMRs) with MethylKit and Methpipe


## Downstream post-processing analysis (enrichment at regulatory features, Oxidative products DIP-seq peaks after Tet and Tdg knock out, Rloops etc..)

