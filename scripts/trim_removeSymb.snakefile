# IN PROGRESS

# This Snakefile runs the pipeline for processing 2bRAD data from Stephanocoenia intersepta samples.
# Processing produces trimmed, deduplicated reads with symbiont reads removed.

# Scripts needed for this pipeline:
# From https://github.com/z0on/2bRAD_denovo
#  - 2bRAD_trim_launch_dedup.pl
#  - trim2bRAD_2barcodes_dedup.pl

#####################################################
# Snakemake rule pipeline: 
# fastq_raw -> make_trim_commands -> execute_trims -> fastqc_trimmed -> create_symb_reference -> align_symbiont -> exclude_symbiont_reads
#####################################################

# CONFIGURATION
configfile: "config.yaml"
# Raw fastq file directory
rawfq_dir = config["rawfq_dir"] 
# Raw QC file output directory
rawqc_dir = config["rawqc_dir"]  
# Trimmed fastq file directory
trimfq_dir = config["trimfq_dir"] 
# Trimmed QC file output directory
trimqc_dir = config["trimqc_dir"] 
# Excluded symbiont trimmed fastq directory
trimfq_noSymb_dir = config["trimfq_noSymb_dir"] 
# Merged unique reads directory
merge_dir = config["merge_dir"] 
# Aligned symbiont bam directory
symb_align_dir = config["symb_align_dir"]
# Kraken database directory
kraken_db = config["kraken_db"]
# De novo reference genome directory
denovo_ref_dir = config["denovo_ref_dir"]
# Script directory
scripts_dir = config["scripts_dir"] 
# Symbiont genomes
symbiont_genomes = config["symbiont_genomes"]
# Symbiont concatenated reference genome
reference_symbiont = config["reference_symbiont"]   
# Number of threads to use
threads: config["threads"]  

# Get all FASTQ files and extract sample names
SAMPLES=glob_wildcards(f"{rawfq_dir}/{{sample}}.fastq").sample

# Define final target(s)
# Temprary rule_all for debugging
rule all:
    input:
        expand(f"{rawqc_dir}/{{sample}}_fastqc.html", sample=SAMPLES),
        expand(f"{trimqc_dir}/{{sample}}_fastqc.html", sample=SAMPLES),
        expand(f"{trimfq_noSymb_dir}/{{sample}}.fastq", sample=SAMPLES)

#############################################
# RULES
#############################################

##############################
# Quality control and trimming
##############################

# Initial fastqc on raw reads
rule fastqc_raw:
    input:
        f"{rawfq_dir}/{{sample}}.fastq"
    output:
        f"{rawqc_dir}/{{sample}}_fastqc.html"
    shell:
        "fastqc {input} -o {output}"

# Multi-QC on raw reads
rule multiqc_raw:
    input:
        expand(f"{rawqc_dir}/{{sample}}_fastqc.html", sample=SAMPLES)
    output:
        f"{rawqc_dir}/multiqc_report.html"
    shell:
        "multiqc {rawqc_dir} -o {rawqc_dir}"

# Generate trimming commands using custom perl script
rule make_trim_commands:
    input:
        f"{rawfq_dir}/{{sample}}.fastq"
    output:
        f"{scripts_dir}/trims.sh"
    params:
        site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}",
        adaptor="AGATC?",
        sampleID=100, # This tells the script to name files from whole original file name
        barcode="[ATGC]{4}"
    shell:
        """
        perl {scripts_dir}/2bRAD_trim_launch_dedup.pl {input} \
            site='{params.site}' \
            adaptor='{params.adaptor}' \
            sampleID='{params.sampleID}' \
            barcode2='{params.barcode}' > {output}
        """

# Execute trimming commands
rule execute_trims:
    input:
        f"{scripts_dir}/trims.sh"
    output:
        f"{trimfq_dir}/{{sample}}.tr0"
    shell:
        "bash {input}"

# Quality filter with cutadapt
rule quality_filter:
    input:
        f"{trimfq_dir}/{{sample}}.tr0"
    output:
        f"{trimfq_dir}/{{sample}}.fastq"
    params:
        quality="15,15",
        minlen=36
    shell:
        """
        cutadapt -q {params.quality} -m {params.minlen} -o {output} {input}
        """

# Run fastqc on trimmed reads
rule fastqc_trimmed:
    input:
        f"{trimfq_dir}/{{sample}}.fastq"
    output:
        f"{trimqc_dir}/{{sample}}_fastqc.html"
    shell:
        "fastqc {input} -o {output}"

# Multi-QC on trimmed reads
rule multiqc_trimmed:   
    input:
        expand(f"{trimqc_dir}/{{sample}}_fastqc.html", sample=SAMPLES)
    output:
        f"{trimqc_dir}/multiqc_report.html"
    shell:
        "multiqc {trimqc_dir} -o {trimqc_dir}"

##############################
# Symbiont read removal
##############################

# Create concatenated symbiont reference genome from individual symbiont genomes
rule create_symb_reference:
    input:
        config["symbiont_genomes"]
    output:
        config["reference_symbiont"]
    shell:
        """
        cat {input} > {output}
        """

# Align reads to symbiont reference genome
rule align_symbiont:
    input:
        fq = f"{trimfq_dir}/{{sample}}.fastq",
        ref = config["reference_symbiont"]
    output:
        bam = f"{symb_align_dir}/{{sample}}_symb.bam",
        bai = f"{symb_align_dir}/{{sample}}_symb.bam.bai"
    params:
        threads = config["threads"]
    shell:
        """
        bwa mem -t {params.threads} {input.ref} {input.fq} | \
        samtools sort -@ {params.threads} -o {output.bam}
        samtools index {output.bam}
        """

# Exclude reads that align to symbiont genomes
rule exclude_symbiont_reads:
    input:
        fq = f"{trimfq_dir}/{{sample}}.fastq",
        symb_bam = f"{symb_align_dir}/{{sample}}_symb.bam"
    output:
        f"{trimfq_noSymb_dir}/{{sample}}.fastq"
    shell:
        """
        samtools view -b -f 4 {input.symb_bam} | \
        samtools fastq - > {output}
        """
