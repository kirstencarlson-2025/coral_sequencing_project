# ------------------------------------------------ #
# 2bRAD SRA Raw Fastq Download Snakefile
# Kirsten Carlson
# Updated 11/2025
# ------------------------------------------------ #

# This Snakefile downloads raw fastq files from NCBI.

# ------------------------------------------------ #
# CONFIGURATION
# ------------------------------------------------ #

configfile: "config.yaml"
# PRJNA accession number
prjna = config["prjna"]
# Conda environment
env = config["env"]
# Resource directory (SRR numbers, run info)
resource_dir = config["resource_dir"]
# Raw fastq file directory
rawfq_dir = config["rawfq_dir"] 
# Script directory
scripts_dir = config["scripts_dir"]

SAMPLES = {}
with open(f"{resource_dir}/srr_2brad_list.txt") as f:
    for line in f:
        srr, sampleid = line.strip().split(",")
        SAMPLES[sampleid] = srr

# Define final target(s)
rule all:
    input:
        expand(f"{rawfq_dir}/{{sample}}.fastq", sample=SAMPLES.keys())

# ------------------------------------------------ #
# RULES
# ------------------------------------------------ #

# Download raw fastq files from SRA
# ------------------------------------------------
rule download_sra:
    output:
        fq = f"{rawfq_dir}/{{sample}}.fastq"
    params:
        srr = lambda wildcards: SAMPLES[wildcards.sample]
    conda:
        config["env"]
    log:
        f"{resource_dir}/download.log"
    threads: 4   
    shell:
        """
        fasterq-dump {params.srr} \
            --threads {threads} \
            -o {output.fq} \
            &> {log}
        """
