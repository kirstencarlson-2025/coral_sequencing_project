# ------------------------------------------------ #
# 2bRAD Trim and QC Snakefile
# Kirsten Carlson
# Updated 11/2025
# ------------------------------------------------ #

# This Snakefile runs the pipeline for processing 2bRAD data. Processing produces trimmed, deduplicated reads.

# Scripts needed for this pipeline:
# From https://github.com/z0on/2bRAD_denovo
#  - trim2bRAD_2barcodes_dedup.pl

# ------------------------------------------------ #
# CONFIGURATION
# ------------------------------------------------ #

configfile: "config.yaml"
# Conda environment
env: config["env"]
# Raw fastq file directory
rawfq_dir = config["rawfq_dir"] 
# Raw QC file output directory
rawqc_dir = config["rawqc_dir"]  
# Trimmed fastq file directory
trimfq_dir = config["trimfq_dir"] 
# Trimmed QC file output directory
trimqc_dir = config["trimqc_dir"] 
# Script directory
scripts_dir = config["scripts_dir"]
# Number of threads to use
threads = config["threads"] 
quality = config["quality"]
minlen = config["minlen"] 

# Get all FASTQ files and extract sample names
SAMPLES=glob_wildcards(f"{rawfq_dir}/{{sample}}.fastq").sample
print("Found samples:", SAMPLES)


# ------------------------------------------------ #
# RULES
# ------------------------------------------------ #

# Initial fastqc on raw reads
# ------------------------------------------------ 
rule fastqc_raw:
    input:
        f"{rawfq_dir}/{{sample}}.fastq"
    output:
        f"{rawqc_dir}/{{sample}}_fastqc.zip",
        f"{rawqc_dir}/{{sample}}_fastqc.html"
    shell:
        "fastqc {input} -o {rawqcq_dir}"

# Multi-QC on raw reads
# ------------------------------------------------ 
rule multiqc_raw:
    input:
        expand(f"{rawqc_dir}/{{sample}}_fastqc.zip", sample=SAMPLES)
    output:
        f"{rawqc_dir}/multiqc_report.html"
    shell:
        "multiqc {input} -o {rawqc_dir}"


# Trim using custom perl script
# This step appears to have already been done before uploading to NCBI SRA, so this is just for reference.
# ------------------------------------------------
rule trim_dedup:
    input:
        f"{rawfq_dir}/{{sample}}.fastq"
    output:
        f"{trimfq_dir}/{{sample}}.tr0"
    params:
        site= lambda wc: ".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}",
        adaptor="AGATC?",
        sampleID=100,
        barcode= lambda wc: "[ATGC]{4}"
    shell:
        """
        perl {scripts_dir}/trim2bRAD_2barcodes_dedup.pl \
        input={input} \
        site="{params.site}" \
        adaptor="{params.adaptor}" \
        sampleID={params.sampleID} \
        deduplicate=1 \
        bc="{params.barcode}"
        """

# Quality filter with cutadapt
# If you need to run trim2bRAD_2barcodes_dedup.pl, run this rule after that step, changing input to the trimmed output file (.tr0).
# This step appears to have already been done before uploading to NCBI SRA, so this is just for reference.
# ------------------------------------------------
rule quality_filter:
    input:
        f"{rawfq_dir}/{{sample}}.fastq"
    output:
        f"{trimfq_dir}/{{sample}}.fastq"
    params:
        quality=config["quality"],
        minlen=config["minlen"]
    shell:
        """
        cutadapt -q {params.quality} -m {params.minlen} -o {output} {input}
        """

# Run fastqc on trimmed reads
# Above step did not trim or remove any reads, so this is just for reference.
# ------------------------------------------------
rule fastqc_trimmed:
    input:
        f"{trimfq_dir}/{{sample}}.fastq"
    output:
        f"{trimqc_dir}/{{sample}}_fastqc.zip",
        f"{trimqc_dir}/{{sample}}_fastqc.html"
    shell:
        "fastqc {input} -o {trimqc_dir}"

# Multi-QC on trimmed reads
# Above step did not trim or remove any reads, so this is just for reference.
# ------------------------------------------------
rule multiqc_trimmed:   
    input:
        expand(f"{trimqc_dir}/{{sample}}_fastqc.zip", sample=SAMPLES)
    output:
        f"{trimqc_dir}/multiqc_report.html"
    shell:
        "multiqc {trimqc_dir} -o {trimqc_dir}"
