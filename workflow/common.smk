
def which_fastq(wildcards):
    # Determines the path to fastq files depending on whether they've
    # been trimmed
    if config['trim']:
        return(f"results/trimmed/{wildcards['sample']}.fastq.gz")
    else:
        return(f"resources/reads/{wildcards['sample']}.fastq.gz")

def which_fastq_r1(wildcards):
    # Determines the path to fastq files depending on whether they've
    # been trimmed
    if config['trim']:
        return(f"results/trimmed/{wildcards['sample']}.R1.fastq.gz")
    else:
        return(f"resources/reads/{wildcards['sample']}.R1.fastq.gz")

def which_fastq_r2(wildcards):
    # Determines the path to fastq files depending on whether they've
    # been trimmed
    if config['trim']:
        return(f"results/trimmed/{wildcards['sample']}.R2.fastq.gz")
    else:
        return(f"resources/reads/{wildcards['sample']}.R2.fastq.gz")
