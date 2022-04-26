
def which_fastq(wildcards):
    # Determines the path to fastq files depending on whether they've
    # been trimmed
    if config['trim']:
        return(f"results/trimmed/{wildcards['sample']}.fastq.gz")
    else:
        return(f"resources/reads/{wildcards['sample']}.fastq.gz")