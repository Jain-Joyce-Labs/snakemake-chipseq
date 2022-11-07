# Snakemake workflow: ChIP-Seq

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/Jain_Joyce-Labs/snakemake-chipseq/workflows/Tests/badge.svg?branch=main)](https://github.com/Jain_Joyce-Labs/snakemake-chipseq/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for processing ChIP-seq experiments. Input is demultiplexed FASTQ sequencing files, output is bedgraph files.

## Steps

This is a general summary of the pipeline with example commands, for running these steps without snakemake. See "Usage" below for instructions on running the pipeline automatically.

**1. Quality check.** Uses [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to generate reports about the sequencing quality in a single sample. Saves results to a directory called `qc/`.

```sh
fastqc sample1.fastq.gz --noextract --threads 4 --outdir qc/
```

**2. Trimming.** This uses [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) to remove low-quality reads and trim uncertain bases from the ends of reads in a single sample and generates a new FASTQ file in a directory called `trimmed/`. These settings are very flexible and will probably change depending on what you see in the fastqc report. Note that this command is for paired-end reads, and expects to be given both forward and reverse FASTQ files in a single command. You will probably want to repeat step 1 with the trimmed FASTQ files to confirm the trimming was effective.

```sh
java -jar /Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 sample1.R1.fastq.gz sample1.R2.fastq.gz trimmed/sample1.R1.fastq.gz trimmed/sample1.R2.fastq.gz SLIDINGWINDOW:4:20 MINLEN:20
```

**3. Alignment.** This command uses [Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to map paired-end reads to a reference genome (in this example, human assembly hg38) and create a SAM file. Reference genomes ("indices") can be downloaded from [the Bowtie homepage](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml).

```sh
bowtie2 -p 4 -q --local -x GRCh38/GRCh38_noalt_as -1 trimmed/sample1.R1.fastq.gz -2 trimmed/sample1.R2.fastq.gz -X 300 -S sample1.sam
```

**4. SAM conversion.** SAM files are inefficient to access and take up a lot of space, so we use [samtools](http://www.htslib.org/) to convert the output to BAM files.

```sh
samtools sort -l 9 -O bam -@ 4 -o sample1.bam sample1.sam
```

**5. Remove duplicate reads.** There is some debate about whether this step is required. In any case, it uses [Picard](https://broadinstitute.github.io/picard/) to remove duplicate reads found in the BAM file, with the assumption that they are actually PCR artifacts and should not be counted multiple times.

```sh
java -jar /usr/picard/picard.jar MarkDuplicates I=sample1.bam O=sample1.nodupes.bam M=sample1.txt REMOVE_DUPLICATES=true
```

**6. Strip blacklisted regions.** This uses [bedtools](https://bedtools.readthedocs.io/en/latest/) to remove regions of the genome determined to be "problematic" by the ENCODE Project because of "anomalous, unstructured, or high signal in next-generation sequencing experiments independent of cell line or experiment." Different versions and organisms can be downloaded from places [such as GitHub](https://github.com/Boyle-Lab/Blacklist/tree/master/lists).

```sh
bedtools subtract -a sample1.nodupes.bam -b hg38-blacklist.v2.bed > sample1.final.bam
```

**7. Index BAM file.** Analyzing the BAM files require they be indexed. This uses samtools again, this time to create an index file in the same directory as the BAM.

```sh
samtools index -@ 4 sample1.bam
```

**8. Generate bedgraph file.** This step breaks the genome into equally sized bins (in the example, 10 kilobases) and counts how many reads fall into those bins. Critically, **this step also does input normalization of the ChIP-seq data**: In the version of the command below, each bin is assigned two numbers: The count of reads in the IP sample, and the number of reads in the "input" sample, in which the immunoprecipitation step has been skipped. The final output is the log2-fold change between these two numbers.

```sh
bamCompare -b1 sample1.bam -b2 sample1.INPUT.bam --minMappingQuality 15 --binSize 10000 --numberOfProcessors 4 --outFileFormat bedgraph --outFileName sample1.bedgraph
```

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
