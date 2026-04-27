# ------------------------------------------------ #
# 2bRAD SRA Raw Fastq Download Snakefile
# Kirsten Carlson
# Updated 1/2026
# ------------------------------------------------ #

# This Snakefile downloads raw fastq files from NCBI.

# ------------------------------------------------ #
# CONFIGURATION
# ------------------------------------------------ #

configfile: "../config/config.yaml"
# Conda environment
env = config["env"]
# List of 2bRAD srr numbers and sample IDs
srr_2brad_list = config["srr_2brad_list"]   
# Resource directory (SRR numbers, run info)
resource_dir = config["resource_dir"]
# Raw fastq file directory
rawfq_dir = config["rawfq_dir"] 

# ------------------------------------------------ #
# RULES
# ------------------------------------------------ #

# Download raw fastq files from SRA
# ------------------------------------------------
rule download_sra:
    output:
        f"{rawfq_dir}/{{sample}}.fastq"
    params:
        SRR = lambda wildcards: SRR_MAP[wildcards.sample]
    conda:
        config["env"]
    threads: 4
    shell:
        """
        fasterq-dump {params.SRR} -O {rawfq_dir} --threads {threads}
        mv {rawfq_dir}/{params.SRR}.fastq {output}
        """
