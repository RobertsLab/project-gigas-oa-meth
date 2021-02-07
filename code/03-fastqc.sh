#!/bin/bash
## Job Name
#SBATCH --job-name=manchester-fastqc2
## Allocation Definition
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=1-00:00:00
## Memory per node
#SBATCH --mem=200G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yaaminiv@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore

### Prepare script ###

# Exit script if any command fails
set -e

# Paths to programs
fastqc=/gscratch/srlab/programs/fastqc_v0.11.9/fastqc
multiqc=/gscratch/srlab/programs/anaconda3/bin/multiqc

# Set CPU threads to use
threads=28

### FastQC for first trim ###

# Populate array with FastQ files
fastq_array=(/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore/*.fq.gz)

# Pass array contents to new variable
fastqc_list=$(echo "${fastq_array[*]}")

# Run FastQC
# NOTE: Do NOT quote ${fastqc_list}
${fastqc} \
--threads ${threads} \
--outdir /gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore \
${fastqc_list}

echo FastQC 1 done

#MultiQC
${multiqc} \
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore/.

echo MultiQC 1 complete

### FastQC for second trim ###

# Populate array with FastQ files
fastq_array=(/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/*.fq.gz)

# Pass array contents to new variable
fastqc_list=$(echo "${fastq_array[*]}")

# Run FastQC
# NOTE: Do NOT quote ${fastqc_list}
${fastqc} \
--threads ${threads} \
--outdir /gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2 \
${fastqc_list}

#MultiQC: Files after second round of adapter trimming
${multiqc} \
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/.

echo MultiQC2 complete

### FastQC for poly-G trim ###

# Populate array with FastQ files
fastq_array=(/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/poly-G/*.fq.gz)

# Pass array contents to new variable
fastqc_list=$(echo "${fastq_array[*]}")

# Run FastQC
# NOTE: Do NOT quote ${fastqc_list}
${fastqc} \
--threads ${threads} \
--outdir /gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/poly-G \
${fastqc_list}

#MultiQC: Files after poly-G tail trimming
/${multiqc}\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/poly-G/.

echo MultiQC 3 complete
