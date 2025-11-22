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
# Environment
bowtie2_env = config["bowtie2_env"]
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
print("Found samples:", SAMPLES)

# Define final target(s)
# Temprary rule_all for debugging
rule all:
    input:
        expand(f"{trimfq_dir}/{{sample}}.fastq", sample=SAMPLES)

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
        f"{rawqc_dir}/{{sample}}_fastqc.zip",
        f"{rawqc_dir}/{{sample}}_fastqc.html"
    shell:
        "fastqc {input} -o {rawqc_dir}"

# Multi-QC on raw reads
rule multiqc_raw:
    input:
        expand(f"{rawqc_dir}/{{sample}}_fastqc.zip", sample=SAMPLES)
    output:
        f"{rawqc_dir}/multiqc_report.html"
    shell:
        "multiqc {input} -o {rawqc_dir}"

# Trim using custom perl script
# This step was already done before uploading to NCBI SRA, so this is just for reference.
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
rule quality_filter:
    input:
        f"{rawfq_dir}/{{sample}}.fastq"
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
        f"{config["reference_symbiont"]}.fasta"
    shell:
        """
        cat {input} > {output}
        """

# Build bowtie2 index for symbiont reference genome
rule index_symbiont_reference:
    input:
        fasta = f"{config['reference_symbiont']}.fasta",
        index = f"{config['reference_symbiont']}"
    output:
        f"{config['reference_symbiont']}.1.bt2"
    shell:
        """
        bowtie2-build {input.fasta} {input.index}
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
