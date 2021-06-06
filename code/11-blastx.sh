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
#SBATCH --time=3-00:00:00
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
blast=/gscratch/srlab/programs/ncbi-blast-2.10.1+/bin/

# Create database from Uniprot-SwissProt information
# Speciy protein database
# Save to srlab/yaamini/blastdbs

${blast}makeblastdb \
-in /gscratch/srlab/yaaminiv/blastdbs/20210601-uniprot_sprot.fasta \
-dbtype prot \
-out /gscratch/srlab/yaaminiv/blastdbs/20210601-uniprot-sprot.blastdb

# Run blastx
# Output format 6 produces a standard BLAST tab-delimited file
${blast}blastx \
-db /gscratch/srlab/yaaminiv/blastdbs/20210601-uniprot-sprot.blastdb \
-query /gscratch/scrubbed/yaaminiv/data/cgigas_uk_roslin_v1_genomic-mito.fa \
-out 20210605-cgigas-roslin-mito-blastx.outfmt6 \
-outfmt 6 \
-evalue 1e-4 \
-max_target_seqs 1 \
-num_threads 28
