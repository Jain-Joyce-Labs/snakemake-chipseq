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

# we can't run snakemake from the head node, so we
# jump into an interactive session to start it:
bsub -Is bash

# wait for the session to start...

source venv/bin/activate
pip install snakemake # this line only needs to be run the first time the pipeline is installed

snakemake --use-singularity --jobs 45 --singularity-args "--bind /home/rabdill/snakemake-chipseq/resources" --cluster 'bsub -n 16 -o {log} -M {resources.memory}000 -R "span[hosts=1]"'
```

If you want to run snakemake in headless mode (or you can't sit around staring at the interactive job until the whole pipeline finishes), you can put the above code, starting after the line "wait for the session to start..." in a shell script called `submit.sh` and run this instead:

```sh
bsub -o snakemake.log -e snakemake.err bash submit.sh
```

The snakemake output you would usually monitor in the terminal will instead be written to a log file at `.snakemake/log/`.

The usage of this workflow may eventually be described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=Jain_Joyce-Labs%2Fsnakemake-chipseq).
