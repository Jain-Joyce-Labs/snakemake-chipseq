# Manual pipeline

This is a general summary of the pipeline with example commands, for running these steps without snakemake. See "Usage" in the readme file for instructions on running the pipeline automatically.

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
