#!/bin/bash
## Job Name
#SBATCH --job-name=20190903_bismark_yaamini
## Allocation Definition
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=3-00:00:00
## Memory per node
#SBATCH --mem=120G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yaaminiv@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/yaamini/analyses/Gigas-WGBS/2019-09-03-Gigas-Bismark

#Genome preparation
#Sam White prepared the bisulfite genome: https://robertslab.github.io/sams-notebook/2019/02/21/Data-Management-Create-C.gigas-Bisulfite-Genome-with-Bismark-on-Mox.html

#Alignment

/gscratch/srlab/programs/Bismark-0.21.0/bismark \
--path_to_bowtie2 /gscratch/srlab/programs/bowtie2-2.3.4.1-linux-x86_64/ \
--samtools_path /gscratch/srlab/programs/samtools-1.9/ \
--non_directional \
-p 4 \
-score_min L,0,-0.9 \
--genome /gscratch/scrubbed/yaamini/data/Gigas-WGBS/2019-09-03-Bismark-Inputs/Crassostrea_gigas.oyster_v9.dna_sm.toplevel/ \
-1 /gscratch/scrubbed/yaamini/data/Gigas-WGBS/2019-09-03-Trimmed-Files/YRVL_R1_001.fastq.gz \
-2 /gscratch/scrubbed/yaamini/data/Gigas-WGBS/2019-09-03-Trimmed-Files/YRVL_R2_001.fastq.gz