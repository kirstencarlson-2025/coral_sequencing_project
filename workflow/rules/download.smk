# ------------------------------------------------ #
# 2bRAD SRA Raw Fastq Download Snakefile
# Kirsten Carlson
# Updated 1/2026
# ------------------------------------------------ #

# This Snakefile downloads raw fastq files from NCBI.

# ------------------------------------------------ #
# CONFIGURATION
# ------------------------------------------------ #

configfile: "../../config/config.yaml"
# Conda environment
env = config["env"]
# List of 2bRAD srr numbers and sample IDs
srr_2brad_list = config["srr_2brad_list"]
# Resource directory (SRR numbers, run info)
resource_dir = config["resource_dir"]
# Raw fastq file directory
rawfq_dir = config["rawfq_dir"] 

import pandas as pd

# Read SRR list
srr_df = pd.read_csv(config["srr_2brad_list"], header=None, names=["SRR", "SAMPLEID"])

SRR = srr_df["SRR"].tolist()
SAMPLEID = srr_df["SAMPLEID"].tolist()


# ------------------------------------------------ #
# RULES
# ------------------------------------------------ #

# Download raw fastq files from SRA
# ------------------------------------------------
rule download_sra:
    output:
        fastq = f"{rawfq_dir}/{{sample}}.fastq"
    params:
        srr = lambda wildcards: srr_df.loc[srr_df["SAMPLEID"] == wildcards.sample, "SRR"].values[0]
    conda:
        config["env"]
    threads: 4   
    shell:
        """
        fasterq-dump {params.srr} \
            --threads {threads} \
            -o {output.fastq}
        """
