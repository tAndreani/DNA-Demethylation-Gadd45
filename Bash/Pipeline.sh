#!/bin/sh

#SBATCH --job-name=DMRtko        # job name
##SBATCH --array=401-494
#SBATCH --nodes=1                            # nodes
#SBATCH -p andrade                              # queue
#SBATCH -A jgu-cbdm
#SBATCH -c 5
#SBATCH --mem=10000M                          # memory
#SBATCH --time=72:00:00                      # time
#SBATCH --error=DMR.tko.err                 # error file name
#SBATCH --output=DMR.tko.out                # output file name
#SBATCH --mail-user=t.andreani@imb-mainz.de  # email
#SBATCH --mail-type=ALL                      # type notification


module unuse /cluster/easybuild/modules/all
module use /cluster/easybuild/nehalem/modules/all 
module load bio/Bowtie2/2.3.2-intel-2017.02
module load bio/Bowtie/1.1.2-intel-2017.02
module load bio/SAMtools/1.5-foss-2017a
module load lang/R/3.4.1-foss-2017a
module load bio/MethPipe/3.4.3-intel-2017.02
module load lang/Perl/5.24.1-intel-2017.02

genome="/project/jgu-cbdm/andradeLab/scratch/tandrean/Genomes/Mus_musculus_mm10/UCSC/mm10/Sequence/WholeGenomeFasta"
files="/project/jgu-cbdm/andradeLab/scratch/tandrean/Data/WGBS/test/Fastq.Trimmed"
genome2="/project/jgu-cbdm/andradeLab/scratch/tandrean/Genomes/Mus_musculus_mm10.2/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta"
testSamples="/project/jgu-cbdm/andradeLab/scratch/tandrean/Data/WGBS/test/subset.samples.test/"
testQval="/project/jgu-cbdm/andradeLab/scratch/tandrean/Data/WGBS/test/subset.samples.test/q.val.40/"
tiles="/project/jgu-cbdm/andradeLab/scratch/tandrean/Data/WGBS/test/Fastq.Trimmed/Remove.Tiles"
clone6R1="/project/jgu-cbdm/andradeLab/scratch/tandrean/Data/WGBS/test/Fastq.Trimmed/Split.AS-165429-LR-25517.R1"
clone6R2="/project/jgu-cbdm/andradeLab/scratch/tandrean/Data/WGBS/test/Fastq.Trimmed/Split.AS-165429-LR-25517.R2"
seqtk="/project/jgu-cbdm/andradeLab/scratch/tandrean/Data/WGBS/test/seqtk/seqtk"
gadd45Tko1="/project/jgu-cbdm/andradeLab/scratch/tandrean/Data/WGBS/Gadd45.TKO/Gadd45.tko.Quality.Filter/Gadd45.TKO1/"
gadd45Tko2="/project/jgu-cbdm/andradeLab/scratch/tandrean/Data/WGBS/Gadd45.TKO/Gadd45.tko.Quality.Filter/Gadd45.TKO2/"
gadd45Tko3="/project/jgu-cbdm/andradeLab/scratch/tandrean/Data/WGBS/Gadd45.TKO/Gadd45.tko.Quality.Filter/Gadd45.TKO3/"
mESC="/project/jgu-cbdm/andradeLab/scratch/tandrean/Data/WGBS/Gadd45.TKO/Gadd45.tko.Quality.Filter/mESC/"
extraction="/project/jgu-cbdm/andradeLab/scratch/tandrean/Data/WGBS/Gadd45.TKO/Gadd45.tko.Quality.Filter/Extraction/"
bismark="/project/jgu-cbdm/andradeLab/scratch/tandrean/Tools_Script/bismark/bismark_v0.18.0/"
fasta="/project/jgu-cbdm/andradeLab/scratch/tandrean/Genomes/Mus_musculus_mm10.2/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome"

#export fastq_file=`sed -n "$SLURM_ARRAY_TASK_ID"p list.files.gadd45.1`
#export fastq_file=`sed -n "$SLURM_ARRAY_TASK_ID"p list.files.gadd45.2`
#export fastq_file=`sed -n "$SLURM_ARRAY_TASK_ID"p list.files.gadd45.3`
#export fastq_file=`sed -n "$SLURM_ARRAY_TASK_ID"p list.files.mESC`



#Select reads with a quality value > 20 All the Reads
../../trim_galore --paired --trim1 $gadd45Tko1/AS-180617-LR-27523_R1.fastq.gz $gadd45Tko1/AS-180617-LR-27523_R2.fastq.gz
../../trim_galore --paired --trim1 $gadd45Tko2/AS-180618-LR-27524_R1.fastq.gz $gadd45Tko2/AS-180618-LR-27524_R2.fastq.gz
../../trim_galore --paired --trim1 $gadd45Tko3/AS-180619-LR-27525_R1.fastq.gz $gadd45Tko3/AS-180619-LR-27525_R2.fastq.gz
../../trim_galore --paired --trim1 $mESC/AS-180620-LR-28032_R1.fastq.gz $mESC/AS-180620-LR-28032_R2.fastq.gz

#Alignment

#Gadd45.TKO1 496 files
$bismark/bismark --bowtie2 -n 1 -I 0 -X 1000 --score_min L,0,-0.6 $genome/genome -1 $gadd45Tko1/R1/$fastq_file -2 $gadd45Tko1/R2/$fastq_file -o $gadd45Tko1/Alignment/ 

#Gadd45.TKO2 494 files
$bismark/bismark --bowtie2 -n 1 -I 0 -X 1000 --score_min L,0,-0.6 $genome/genome -1 $gadd45Tko2/R1/$fastq_file -2 $gadd45Tko2/R2/$fastq_file -o $gadd45Tko2/Alignment/

#Gadd45.TKO3 489 files
$bismark/bismark --bowtie2 -n 1 -I 0 -X 1000 --score_min L,0,-0.6 $genome/genome -1 $gadd45Tko3/R1/$fastq_file -2 $gadd45Tko3/R2/$fastq_file -o  $gadd45Tko3/Alignment/

#mESC Wild Type 502 files
$bismark/bismark --bowtie2 -n 1 -I 0 -X 1000 --score_min L,0,-0.6 $genome/genome -1 $mESC/R1/$fastq_file -2 $mESC/R2/$fastq_file -o $mESC/Alignment/


#sort for bismark by name always
samtools merge -n $gadd45Tko1/Alignment/Gadd45.tko1.bam $gadd45Tko1/Alignment/x*.bam
samtools merge -n $gadd45Tko2/Alignment/Gadd45.tko2.bam $gadd45Tko2/Alignment/x*.bam
samtools merge -n $gadd45Tko3/Alignment/Gadd45.tko3.bam $gadd45Tko3/Alignment/x*.bam
samtools merge -n $mESC/Alignment/mESC.bam $mESC/Alignment/x*.bam


#Extract
$bismark/bismark_methylation_extractor -p --ignore 5 --ignore_r2 5 --ample_memory --bedGraph --counts --cytosine_report --buffer_size 10G --genome_folder $genome/genome $extraction/Gadd45.tko1.bam
$bismark/bismark_methylation_extractor -p --ignore 5 --ignore_r2 5 --ample_memory --bedGraph --counts --cytosine_report --buffer_size 10G --genome_folder $genome/genome $extraction/Gadd45.tko2.bam
$bismark/bismark_methylation_extractor -p --ignore 5 --ignore_r2 5 --ample_memory --bedGraph --counts --cytosine_report --buffer_size 10G --genome_folder $genome/genome $extraction/Gadd45.tko3.bam
$bismark/bismark_methylation_extractor -p --ignore 5 --ignore_r2 5 --ample_memory --bedGraph --counts --cytosine_report --buffer_size 10G --genome_folder $genome/genome $extraction/mESC.bam


#Methpipe conversion
##Convert files format
to-mr -o Gadd45.tko1.mr -m bismark Gadd45.tko1.bam
to-mr -o Gadd45.tko2.mr -m bismark Gadd45.tko2.bam
to-mr -o Gadd45.tko3.mr -m bismark Gadd45.tko3.bam
to-mr -o mESC.mr -m bismark mESC.bam


#Sort to remove duplicate
LC_ALL=C sort -T temporary.mESC  -k 1,1 -k 2,2n -k 3,3n -k 6,6 -o mESC.mr.sort.for.duplicate.mr mESC.mr
LC_ALL=C sort -T temporary.tko1  -k 1,1 -k 2,2n -k 3,3n -k 6,6 -o Gadd45.tko1.mr.sort.for.duplicate.mr Gadd45.tko1.mr
LC_ALL=C sort -T temporary.tko2  -k 1,1 -k 2,2n -k 3,3n -k 6,6 -o Gadd45.tko2.mr.sort.for.duplicate.mr Gadd45.tko2.mr
LC_ALL=C sort -T temporary.tko3  -k 1,1 -k 2,2n -k 3,3n -k 6,6 -o Gadd45.tko3.mr.sort.for.duplicate.mr Gadd45.tko3.mr

#Remove duplicate
duplicate-remover -S mESC.mr.sort.for.duplicate.mr.remove.stat -o mESC.mr.sort.for.duplicate.mr.remove mESC.mr.sort.for.duplicate.mr
duplicate-remover -S Gadd45.tko1.mr.sort.for.duplicate.mr.remove.stat -o Gadd45.tko1.mr.sort.for.duplicate.mr.remove Gadd45.tko1.mr.sort.for.duplicate.mr
duplicate-remover -S Gadd45.tko2.mr.sort.for.duplicate.mr.remove.stat -o Gadd45.tko2.mr.sort.for.duplicate.mr.remove Gadd45.tko2.mr.sort.for.duplicate.mr
duplicate-remover -S Gadd45.tko3.mr.sort.for.duplicate.mr.remove.stat -o Gadd45.tko3.mr.sort.for.duplicate.mr.remove Gadd45.tko3.mr.sort.for.duplicate.mr


#Bisulphite Conversion
bsrate -c $genome/genome/genome.fa -o mESC.mr.sort.for.duplicate.mr.remove.bsrate mESC.mr.sort.for.duplicate.mr.remove
bsrate -c $genome/genome/genome.fa -o Gadd45.tko1.mr.sort.for.duplicate.mr.remove.bsrate Gadd45.tko1.mr.sort.for.duplicate.mr.remove
bsrate -c $genome/genome/genome.fa -o Gadd45.tko2.mr.sort.for.duplicate.mr.remove.bsrate Gadd45.tko2.mr.sort.for.duplicate.mr.remove
bsrate -c $genome/genome/genome.fa -o Gadd45.tko3.mr.sort.for.duplicate.mr.remove.bsrate Gadd45.tko3.mr.sort.for.duplicate.mr.remove

#Sort for methylation count
LC_ALL=C sort -T temporary.mESC -k 1,1 -k 3,3n -k 2,2n -k 6,6 -o mESC.WT.sort.for.duplicate.mr.remove.sorted_end_first.mr mESC.mr.sort.for.duplicate.mr.remove 
LC_ALL=C sort -T temporary.tko1 -k 1,1 -k 3,3n -k 2,2n -k 6,6 -o Gadd45.tko1.sort.for.duplicate.mr.remove.sorted_end_first.mr Gadd45.tko1.mr.sort.for.duplicate.mr.remove
LC_ALL=C sort -T temporary.tko2 -k 1,1 -k 3,3n -k 2,2n -k 6,6 -o Gadd45.tko2.sort.for.duplicate.mr.remove.sorted_end_first.mr Gadd45.tko2.mr.sort.for.duplicate.mr.remove
LC_ALL=C sort -T temporary.tko3 -k 1,1 -k 3,3n -k 2,2n -k 6,6 -o Gadd45.tko3.sort.for.duplicate.mr.remove.sorted_end_first.mr Gadd45.tko3.mr.sort.for.duplicate.mr.remove

#Count
methcounts -c $fasta/genome.fa -o mESC.WT.meth mESC.WT.sort.for.duplicate.mr.remove.sorted_end_first.mr
methcounts -c $fasta/genome.fa -o Gadd45.tko1.meth Gadd45.tko1.sort.for.duplicate.mr.remove.sorted_end_first.mr
methcounts -c $fasta/genome.fa -o Gadd45.tko2.meth Gadd45.tko2.sort.for.duplicate.mr.remove.sorted_end_first.mr
methcounts -c $fasta/genome.fa -o Gadd45.tko3.meth Gadd45.tko3.sort.for.duplicate.mr.remove.sorted_end_first.mr


###Here you should merge the files of the tko1, tko2, tko3
merge-methcounts Gadd45.tko1.meth Gadd45.tko2.meth Gadd45.tko3.meth -o Gadd45.tko.meth

#Statistic on all the methylated cytosines
levels -o mESC.WT.meth.stat  mESC.WT.meth
levels -o Gadd45.tko1.meth.stat Gadd45.tko1.meth
levels -o Gadd45.tko2.meth.stat Gadd45.tko2.meth
levels -o Gadd45.tko3.meth.stat Gadd45.tko3.meth
levels -o Gadd45.tko.meth.stat Gadd45.tko.meth


#Take only the CpG
grep CpG mESC.WT.meth > mESC.WT.CpG.meth
grep CpG Gadd45.tko1.meth > Gadd45.tko1.CpG.meth
grep CpG Gadd45.tko2.meth > Gadd45.tko2.CpG.meth
grep CpG Gadd45.tko3.meth > Gadd45.tko3.CpG.meth
grep CpG Gadd45.tko.meth > Gadd45.tko.CpG.meth


#Keep symmetric CpG including the mutated sites
symmetric-cpgs -m -o mESC.WT.CpG.symmetric.and.mutations.meth mESC.WT.CpG.meth
symmetric-cpgs -m -o Gadd45.tko1.CpG.symmetric.and.mutations.meth Gadd45.tko1.CpG.meth
symmetric-cpgs -m -o Gadd45.tko2.CpG.symmetric.and.mutations.meth Gadd45.tko2.CpG.meth
symmetric-cpgs -m -o Gadd45.tko3.CpG.symmetric.and.mutations.meth Gadd45.tko3.CpG.meth
symmetric-cpgs -m -o Gadd45.tko.CpG.symmetric.and.mutations.meth Gadd45.tko.CpG.meth


#Call Hypo Methylated Cystosine
hmr -o mESC.WT.CpG.symmetric.and.mutations.hypo.hmr mESC.WT.CpG.symmetric.and.mutations.meth
hmr -o Gadd45.tko1.CpG.symmetric.and.mutations.hypo.hmr Gadd45.tko1.CpG.symmetric.and.mutations.meth
hmr -o Gadd45.tko2.CpG.symmetric.and.mutations.hypo.hmr Gadd45.tko2.CpG.symmetric.and.mutations.meth
hmr -o Gadd45.tko3.CpG.symmetric.and.mutations.hypo.hmr Gadd45.tko3.CpG.symmetric.and.mutations.meth
hmr -o Gadd45.tko.CpG.symmetric.and.mutations.hypo.hmr Gadd45.tko.CpG.symmetric.and.mutations.meth


#Invert the methylome to call Hyper Methylated Regions
awk '{$5=1-$5; print $0}' mESC.WT.CpG.symmetric.and.mutations.meth > mESC.WT.CpG.symmetric.and.mutations.inverted.meth
awk '{$5=1-$5; print $0}' Gadd45.tko1.CpG.symmetric.and.mutations.meth > Gadd45.tko1.CpG.symmetric.and.mutations.inverted.meth
awk '{$5=1-$5; print $0}' Gadd45.tko2.CpG.symmetric.and.mutations.meth > Gadd45.tko2.CpG.symmetric.and.mutations.inverted.meth
awk '{$5=1-$5; print $0}' Gadd45.tko3.CpG.symmetric.and.mutations.meth > Gadd45.tko3.CpG.symmetric.and.mutations.inverted.meth
awk '{$5=1-$5; print $0}' Gadd45.tko.CpG.symmetric.and.mutations.meth > Gadd45.tko.CpG.symmetric.and.mutations.inverted.meth


#Call Hyper Methylated Cystosine
hmr -o mESC.WT.CpG.symmetric.and.mutations.hyper.hmr mESC.WT.CpG.symmetric.and.mutations.inverted.meth
hmr -o Gadd45.tko1.CpG.symmetric.and.mutations.hyper.hmr Gadd45.tko1.CpG.symmetric.and.mutations.inverted.meth
hmr -o Gadd45.tko2.CpG.symmetric.and.mutations.hyper.hmr Gadd45.tko2.CpG.symmetric.and.mutations.inverted.meth
hmr -o Gadd45.tko3.CpG.symmetric.and.mutations.hyper.hmr Gadd45.tko3.CpG.symmetric.and.mutations.inverted.meth
hmr -o Gadd45.tko.CpG.symmetric.and.mutations.hyper.hmr Gadd45.tko.CpG.symmetric.and.mutations.inverted.meth


#Call Hyper Methylated Cystosines second approach
hypermr -o mESC.WT.hypemr mESC.WT.CpG.symmetric.and.mutations.meth
hypermr -o Gadd45.tko1.hypermr Gadd45.tko1.CpG.symmetric.and.mutations.meth
hypermr -o Gadd45.tko2.hypermr Gadd45.tko2.CpG.symmetric.and.mutations.meth
hypermr -o Gadd45.tko3.hypermr Gadd45.tko3.CpG.symmetric.and.mutations.meth
hypermr -o Gadd45.tko.hypermr Gadd45.tko.CpG.symmetric.and.mutations.meth


#Call Methylation Differences wt vs tko (merged files)
methdiff -o mESC.Vs.Gadd45.tko.hypo.methdiff mESC.WT.CpG.symmetric.and.mutations.meth Gadd45.tko.CpG.symmetric.and.mutations.meth 
methdiff -o mESC.Vs.Gadd45.tko.hyper.methdiff mESC.WT.CpG.symmetric.and.mutations.inverted.meth Gadd45.tko.CpG.symmetric.and.mutations.inverted.meth


#Differential Methylated Regions wt vs tko (merged tko samples) 
#Hyper
dmr mESC.Vs.Gadd45.tko.hyper.methdiff mESC.WT.CpG.symmetric.and.mutations.hyper.hmr Gadd45.tko.CpG.symmetric.and.mutations.hyper.hmr mESC.higher.than.Gadd45.tko Gadd45.tko.higher.than.mESC

#Hypo
dmr mESC.Vs.Gadd45.tko.hypo.methdiff mESC.WT.CpG.symmetric.and.mutations.hypo.hmr Gadd45.tko.CpG.symmetric.and.mutations.hypo.hmr mESC.lower.than.Gadd45.tko Gadd45.tko.lower.than.mESC 

#Hyper in Gadd45 TKO second approach wt vs tko
dmr mESC.Vs.Gadd45.tko.hyper.methdiff mESC.WT.hypemr Gadd45.tko.hypermr mESC.higher.than.Gadd45.tko.version2 Gadd45.tko.higher.than.mESC.version2
