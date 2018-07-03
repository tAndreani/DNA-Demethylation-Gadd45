## Methylome Data Processing and Analysis in mESCs upon Gadd45a,b,g Proteins Knock Out
Collection of scripts to detect Differentially Methylated Regions and perform downstream enrichment analysis at several genomic features of mouse ES cells.

## Quality Control and Trim
Before performing the quality control we want to split the fastq files in small pieces of 1mln reads in order to parallelize all the analysis:

`zcat example.fq.gz | split -4000000 &` 

## Alignment (parallelization) and extraction of CGH, CHG and CHH


## Call of Differentially Methylated Regions (DMRs) with MethylKit and Methpipe


## Downstream post-processing analysis (enrichment at regulatory features, Oxidative products DIP-seq peaks after Tet and Tdg knock out, Rloops etc..)

