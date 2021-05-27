#!/bin/bash
## Job Name
#SBATCH --job-name=roslin-diamond
## Allocation Definition
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=5-00:00:00
## Memory per node
#SBATCH --mem=120G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yaaminiv@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/yaaminiv/data

#Part of functional annotation workflow for C. gigas DML GO term enrichment

# Exit script if any command fails
set -e

# Program paths
diamond=/gscratch/srlab/programs/diamond-2.0.4/diamond

# Run DIAMOND with blastx
# Output format 6 produces a standard BLAST tab-delimited file
${diamond} blastx \
--db /gscratch/srlab/blastdbs/uniprot_sprot_20200123/uniprot_sprot.dmnd \
--query /gscratch/scrubbed/yaaminiv/data/cgigas_uk_roslin_v1_genomic-mito.fa \
--out 20210526-cgigas-roslin-blastx.outfmt6 \
--outfmt 6 \
--evalue 1e-4 \
--max-target-seqs 1 \
--block-size 15.0 \
--index-chunks 4
