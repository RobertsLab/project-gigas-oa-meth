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

#Summary Report

/gscratch/srlab/programs/Bismark-0.21.0/bismark2summary \
YRVA_R1_001_bismark_bt2_pe.bam \
YRVL_R1_001_bismark_bt2_pe.bam
