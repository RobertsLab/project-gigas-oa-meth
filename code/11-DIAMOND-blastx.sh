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
#SBATCH --mem=500G
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

# Create database from Uniprot-SwissProt information
# Use 27 threads to create database
# Save to srlab/yaamini/blastdbs
${diamond} makedb \
--in /gscratch/srlab/yaaminiv/blastdbs/20210601-uniprot_sprot.fasta \
--threads 27 \
-d /gscratch/srlab/yaaminiv/blastdbs/20210601-uniprot-sprot.dmnd

# Run DIAMOND with blastx
# Output format 6 produces a standard BLAST tab-delimited file
${diamond} blastx \
--db /gscratch/srlab/yaaminiv/blastdbs/20210601-uniprot-sprot.dmnd \
--query /gscratch/scrubbed/yaaminiv/data/cgigas_uk_roslin_v1_genomic-mito.fa \
--out 20210601-cgigas-roslin-mito-blastx.outfmt6 \
--long-reads \
--outfmt 6 \
--evalue 1e-4 \
--max-target-seqs 1 \
--block-size 15.0 \
--index-chunks 4
