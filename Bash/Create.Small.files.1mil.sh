#Create files with 1 mil reads from fastq files

#G45 TKO1
mkdir Gadd45.TKO1
cd Gadd45.TKO1/
Gadd45.TKO1.R1.fq.gz R1
Gadd45.TKO1.R2.fq.gz R2
mkdir R1
cd R1/
zcat Gadd45.TKO1.R1.fq.gz | split -4000000 &
mkdir R2
cd R2/
zcat Gadd45.TKO1.R2.fq.gz | split -4000000 &
  
#G45 TKO2
mkdir Gadd45.TKO2
cd Gadd45.TKO2/
mkdir R1
cd R1/
zcat Gadd45.TKO2.R1.fq.gz | split -4000000 &
mkdir R2
cd R2/
zcat Gadd45.TKO2.R2.fq.gz | split -4000000 &
  
  
#G45 TKO3
mkdir Gadd45.TKO3
cd Gadd45.TKO3/
mkdir R1
cd R1/
zcat Gadd45.TKO3.R1.fq.gz | split -4000000 &
mkdir R2
cd R2/
zcat Gadd45.TKO3.R2.fq.gz | split -4000000 &
  
  
#ES cells WT
mkdir ES.wt
cd ES.wt/
mkdir R1
cd R1/
zcat ES.wt.R1.fq.gz | split -4000000 &
mkdir R2
cd R2/
zcat ES.wt.R2.fq.gz | split -4000000 &