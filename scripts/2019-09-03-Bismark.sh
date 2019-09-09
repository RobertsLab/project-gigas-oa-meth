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
#SBATCH --time=00-09:00:00
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

%%bash

find /gscratch/scrubbed/yaamini/data/Gigas-WGBS/2019-09-03-Trimmed-Files/YRV*_R1_001.fastq.gz \
| xargs basename -s _R1_001.fastq.gz | xargs -I{} /gscratch/srlab/programs/Bismark-0.21.0/bismark \
--path_to_bowtie2 /gscratch/srlab/programs/bowtie2-2.3.4.1-linux-x86_64/ \
--samtools_path /gscratch/srlab/programs/samtools-1.9/ \
--non_directional \
-p 4 \
-score_min L,0,-0.9 \
--genome /gscratch/scrubbed/yaamini/data/Gigas-WGBS/2019-09-03-Bismark-Inputs/ \
-1 /gscratch/scrubbed/yaamini/data/Gigas-WGBS/2019-09-03-Trimmed-Files/{}_R1_001.fastq.gz \
-2 /gscratch/scrubbed/yaamini/data/Gigas-WGBS/2019-09-03-Trimmed-Files/{}_R2_001.fastq.gz

#Deduplication

%%bash

/gscratch/srlab/programs/Bismark-0.21.0/deduplicate_bismark \
--samtools_path /gscratch/srlab/programs/samtools-1.9/ \
-p \
--bam \
YRV*_R1_001_bismark_bt2_pe.bam

#Sorting for Downstream Applications

%%bash

find *deduplicated.bam \
| xargs basename -s _R1_001_bismark_bt2_pe.deduplicated.bam | xargs -I{} /gscratch/srlab/programs/samtools-1.9/samtools \
sort {}_R1_001_bismark_bt2_pe.deduplicated.bam \
-o {}_dedup.sorted.bam

#Indexing for Downstream Applications

%%bash

find *dedup.sorted.bam \
| xargs basename -s _dedup.sorted.bam | xargs -I{} /gscratch/srlab/programs/samtools-1.9/samtools \
index {}_dedup.sorted.bam

#Methylation Extraction

%%bash

/gscratch/srlab/programs/Bismark-0.21.0/bismark_methylation_extractor \
--samtools_path /gscratch/srlab/programs/samtools-1.9/ \
-p \
--bedGraph \
--counts \
--scaffolds \
--multicore 4 \
*deduplicated.bam

#HTML Processing Report

%%bash

/gscratch/srlab/programs/Bismark-0.21.0/bismark2report

#Summary Report

%%bash

/gscratch/srlab/programs/Bismark-0.21.0/bismark2summary
