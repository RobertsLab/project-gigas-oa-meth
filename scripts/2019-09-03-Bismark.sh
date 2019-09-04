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
#SBATCH --time=30-00:00:00
## Memory per node
#SBATCH --mem=120G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yaaminiv@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/yaamini/analyses/Gigas-WGBS/2019-09-03-Gigas-Bismark

#Add Directories and Programs to Path

bismark_dir="/gscratch/srlab/programs/Bismark-0.21.0"
bowtie2_dir="/gscratch/srlab/programs/bowtie2-2.3.4.1-linux-x86_64/"
trimmed_files="/gscratch/scrubbed/yaamini/data/Gigas-WGBS/2019-09-03-Trimmed-Files/"
genome="/gscratch/scrubbed/yaamini/data/Gigas-WGBS/2019-09-03-Bismark-Inputs/"
samtools_dir="/gscratch/srlab/programs/samtools-1.9/"

#Genome preparation
#Sam White prepared the bisulfite genome: https://robertslab.github.io/sams-notebook/2019/02/21/Data-Management-Create-C.gigas-Bisulfite-Genome-with-Bismark-on-Mox.html

#Alignment

find ${trimmed_files}YRV*_R1_001.fastq.gz \
| xargs basename -s _R1_001.fastq.gz | xargs -I{} ${bismark_dir}/bismark \
--path_to_bowtie ${bowtie2_dir} \
--samtools_path ${samtools_dir} \
--non_directional \
-p 4 \
-score_min L,0,-0.9 \
--genome ${genome} \
-1 ${trimmed_files}/{}_R1_001.fastq.gz \
-2 ${trimmed_files}/{}_R2_001.fastq.gz \

#Deduplication

${bismark_dir}/deduplicate_bismark \
--samtools_path ${samtools_dir} \
-p \
--bam \
YRV*_R1_001_bismark_bt2_pe.bam

#Sorting for Downstream Applications

find *deduplicated.bam \
| xargs basename -s _R1_001_bismark_bt2_pe.deduplicated.bam | xargs -I{} ${samtools_dir}/samtools \
sort {}_R1_001_bismark_bt2_pe.deduplicated.bam \
-o {}_dedup.sorted.bam

#Indexing for Downstream Applications

find *dedup.sorted.bam \
| xargs basename -s _dedup.sorted.bam | xargs -I{} ${samtools_dir}/samtools \
index {}_dedup.sorted.bam

#Methylation Extraction

${bismark_dir}/bismark_methylation_extractor \
--samtools_path ${samtools_dir} \
-p \
--bedGraph \
--counts \
--scaffolds \
--multicore 4 \
*deduplicated.bam

#HTML Processing Report

${bismark_dir}/bismark2report

#Summary Report

${bismark_dir}/bismark2summary
