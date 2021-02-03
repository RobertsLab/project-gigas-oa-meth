#!/bin/bash
## Job Name
#SBATCH --job-name=manchester-trimgalore
## Allocation Definition
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=5-00:00:00
## Memory per node
#SBATCH --mem=200G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yaaminiv@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore

#Exit script if any command fails
set -e

#TrimGalore: Remove most abundant sequences and hard trim 10 bp from each end
/gscratch/srlab/programs/TrimGalore-0.6.6/trim_galore \
--output_dir /gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore \
--paired \
--fastqc_args \
"--outdir /gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore \
--threads 28" \
--illumina \
--clip_R1 10 \
--clip_R2 10 \
--three_prime_clip_R1 10 \
--three_prime_clip_R2 10 \
--path_to_cutadapt /usr/lusers/yaaminiv/.local/bin/cutadapt \
/gscratch/scrubbed/yaaminiv/Manchester/data/zr3616_1_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Manchester/data/zr3616_1_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Manchester/data/zr3616_2_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Manchester/data/zr3616_2_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Manchester/data/zr3616_3_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Manchester/data/zr3616_3_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Manchester/data/zr3616_4_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Manchester/data/zr3616_4_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Manchester/data/zr3616_5_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Manchester/data/zr3616_5_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Manchester/data/zr3616_6_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Manchester/data/zr3616_6_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Manchester/data/zr3616_7_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Manchester/data/zr3616_7_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Manchester/data/zr3616_8_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Manchester/data/zr3616_8_R2.fq.gz

echo TrimGalore 1 complete

#MultiQC
/gscratch/srlab/programs/anaconda3/bin/multiqc \
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore/.

echo MultiQC 1 complete

#TrimGalore: Remove adapter sequences

/gscratch/srlab/programs/TrimGalore-0.6.6/trim_galore \
--output_dir /gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2 \
--paired \
--fastqc_args \
"--outdir /gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2 \
--threads 28" \
--illumina \
--path_to_cutadapt /usr/lusers/yaaminiv/.local/bin/cutadapt \
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore/zr3616_1_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore/zr3616_1_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore/zr3616_2_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore/zr3616_2_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore/zr3616_3_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore/zr3616_3_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore/zr3616_4_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore/zr3616_4_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore/zr3616_5_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore/zr3616_5_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore/zr3616_6_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore/zr3616_6_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore/zr3616_7_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore/zr3616_7_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore/zr3616_8_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore/zr3616_8_R2_val_2.fq.gz

echo TrimGalore 2 complete

#MultiQC: Files after second round of adapter trimming
/gscratch/srlab/programs/anaconda3/bin/multiqc \
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/.

echo MultiQC2 complete

#TrimGalore: Remove poly-G tails. Sequence obtained from FastQC reports

/gscratch/srlab/programs/TrimGalore-0.6.6/trim_galore \
--output_dir /gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/poly-G \
--paired \
--fastqc_args \
"--outdir /gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/poly-G \
--threads 28" \
--adapter GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG \
--adapter2 GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG \
--path_to_cutadapt /usr/lusers/yaaminiv/.local/bin/cutadapt \
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/zr3616_1_R1_val_1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/zr3616_1_R2_val_2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/zr3616_2_R1_val_1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/zr3616_2_R2_val_2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/zr3616_3_R1_val_1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/zr3616_3_R2_val_2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/zr3616_4_R1_val_1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/zr3616_4_R2_val_2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/zr3616_5_R1_val_1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/zr3616_5_R2_val_2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/zr3616_6_R1_val_1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/zr3616_6_R2_val_2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/zr3616_7_R1_val_1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/zr3616_7_R2_val_2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/zr3616_8_R1_val_1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/zr3616_8_R2_val_2_val_2.fq.gz

echo TrimGalore 3 complete

#MultiQC: Files after poly-G tail trimming
/gscratch/srlab/programs/anaconda3/bin/multiqc \
/gscratch/scrubbed/yaaminiv/Manchester/analyses/trimgalore-2/poly-G/.

echo MultiQC 3 complete
