#!/bin/bash
## Job Name
#SBATCH --job-name=methylKit-manchester
## Allocation Definition
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=10-00:00:00
## Memory per node
#SBATCH --mem=100G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yaaminiv@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/yaaminiv/Manchester/analyses/methylKit

# Load R 3.6.0

module load r_3.6.0

Rscript > output.txt 2>&1 /gscratch/home/yaaminiv/06-methylKit.R
