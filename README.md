# Snakemake workflow: ChIP-Seq

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/Jain_Joyce-Labs/snakemake-chipseq/workflows/Tests/badge.svg?branch=main)](https://github.com/Jain_Joyce-Labs/snakemake-chipseq/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for processing ChIP-seq experiments. Input is demultiplexed FASTQ sequencing files, output is bedgraph files.


## Usage

Running on the HPC cluster:

```sh
git clone https://github.com/Jain-Joyce-Labs/snakemake-chipseq.git NAME_OF_PROJECT
cd NAME_OF_PROJECT

# (Here is where you move the input files into the `resources` directory)

python3 -m venv venv
source venv/bin/activate
pip install snakemake

snakemake --use-singularity --jobs 45 --cluster 'bsub -n 16 -o {log} -M 12000 -R "span[hosts=1]"'
```

The usage of this workflow may eventually be described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=Jain_Joyce-Labs%2Fsnakemake-chipseq).
