#!/bin/sh

#SBATCH --job-name=heatMap        # job name
#SBATCH --array=1-22
#SBATCH --nodes=1                            # nodes
#SBATCH -p andrade                              # queue
#SBATCH -A jgu-cbdm
#SBATCH -c 5
#SBATCH --mem=16000M                          # memory
#SBATCH --time=120:00:00                      # time
#SBATCH --error=heatMap.err                 # error file name
#SBATCH --output=heatMap.out                # output file name
#SBATCH --mail-user=t.andreani@imb-mainz.de  # email
#SBATCH --mail-type=ALL                      # type notification

module load bio/Bowtie2/2.3.2-intel-2017.02
module load bio/Bowtie/1.1.2-intel-2017.02
module load bio/SAMtools/1.5-foss-2017a
module unuse /cluster/easybuild/modules/all
module use /cluster/easybuild/nehalem/modules/all
module load lang/R/3.4.1-foss-2017a
module load lang/Python/2.7.12-foss-2017a

export mark=`sed -n "$SLURM_ARRAY_TASK_ID"p list.markers`

computeMatrix reference-point -S GSM1863079_ESC_DKO_hMeDIP.fastq_monoExt300.mm10.No.Overlap.bw -R Hyper.DMRs.bin.100.delta.30.min.2.CpG.mm10.merged.consecutive.bed -out GSM1863079_ESC_DKO_hMeDIP.in.HyperMe.DMRs.mm10.gz -a 5000 -b 5000 --outFileSortedRegions GSM1863079_ESC_DKO_hMeDIP.in.HyperMe.DMRs.mm10.bed--missingDataAsZero --skipZeros


plotHeatmap -m Matrix.$mark.gz -out Matrix.$mark.png ---colorList cornflowerblue,yellow,red --missingDataColor white
