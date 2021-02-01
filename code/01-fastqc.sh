#!/bin/bash
## Job Name
#SBATCH --job-name=manchester-fastqc
## Allocation Definition
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=3-00:00:00
## Memory per node
#SBATCH --mem=200G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yaaminiv@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/yaaminiv/Manchester/analyses/fastqc

# Directories and programs
reads_dir="/gscratch/scrubbed/yaaminiv/Manchester/data/"

source /gscratch/srlab/programs/scripts/paths.sh

#FastQC: Assess raw sequence quality
find ${reads_dir}* \
| fastqc
