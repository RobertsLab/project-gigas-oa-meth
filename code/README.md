# code

Bash and R Markdown scripts used to process and analyze data. Output from these scripts can be found [here](https://github.com/RobertsLab/project-gigas-oa-meth/tree/master/output).

- [01-fastqc.sh](https://github.com/RobertsLab/project-gigas-oa-meth/blob/master/code/01-fastqc.sh): Mox script used to assess raw data quality
- [02-trimgalore.sh](https://github.com/RobertsLab/project-gigas-oa-meth/blob/master/code/02-trimgalore.sh): Mox script used to trim adapter sequences out of data
- [03-fastqc.sh](https://github.com/RobertsLab/project-gigas-oa-meth/blob/master/code/03-fastqc.sh): Mox script used to assess trimmed data quality after [`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) did not run successfully within `trimgalore`
- [04-bismark.sh](https://github.com/RobertsLab/project-gigas-oa-meth/blob/master/code/04-bismark.sh): Mox script to align sequences to bisulfite-converted genome and extract methylation information
