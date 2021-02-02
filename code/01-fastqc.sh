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

#Exit script if any command fails

set -e

# Directories and programs
reads_dir="/gscratch/scrubbed/yaaminiv/Manchester/data/"
fastqc_dir="/gscratch/srlab/programs/fastqc_v0.11.9/"
source /gscratch/srlab/programs/scripts/paths.sh

# Load Python Mox module for Python module availability
module load intel-python3_2017

#FastQC: Assess raw sequence quality
find ${reads_dir}* \
| ${fastqc_dir}fastqc \
--outdir /gscratch/scrubbed/yaaminiv/Manchester/analyses/fastqc \
--threads 28

#MultiQC: Combine reports
/gscratch/srlab/programs/anaconda3/bin/multiqc \
/gscratch/scrubbed/yaaminiv/Manchester/analyses/fastqc/.
