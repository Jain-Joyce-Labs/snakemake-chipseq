#!/bin/sh
########
# This script can be used to launch snakemake
# on the cluster without an interactive session.
# Use a command something like this:
# bsub -o snakemake.log -e snakemake.err bash submit.sh
########

module load python/3.9.1
python -m venv venv
source venv/bin/activate

snakemake --use-singularity --jobs 45 --singularity-args "--bind /home/rabdill/snakemake-chipseq/resources" --cluster 'bsub -n 16 -o {log} -M {resources.memory}000 -R "span[hosts=1]"'
