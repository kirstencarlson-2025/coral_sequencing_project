# ------------------------------------------------ #
# 2bRAD Get Metadata Snakefile
# Kirsten Carlson
# Updated 11/2025
# ------------------------------------------------ #

# This Snakefile gets SRR numbers from PRJNA accession number.

# ------------------------------------------------ #
# CONFIGURATION
# ------------------------------------------------ #

configfile: "../../config/config.yaml"
# PRJNA accession number
prjna = config["prjna"]
# Conda environment
env = config["env"]
# Resource directory (SRR numbers, run info)
resource_dir = config["resource_dir"]

rule all_metadata:
    input:
        f"{resource_dir}/srr_2brad_list.txt",
        f"{resource_dir}/srr_its2_list.txt"

# ------------------------------------------------ #
# RULES
# ------------------------------------------------ #

# Fetch run info from PRJNA accession number
# ------------------------------------------------
rule fetch_srr:
    output:
        f"{resource_dir}/srr_runinfo.txt"
    params:
        prjna = config["prjna"]
    shell:
        """
        esearch -db sra -query {params.prjna} | efetch -format runinfo > {output}
        """

# Save SRR numbers to a list
# ------------------------------------------------
rule srr_2brad_list:
    input:
        f"{resource_dir}/srr_runinfo.txt"
    output:
        f"{resource_dir}/srr_2brad_list.txt"
    shell:
        """
        grep 2bRAD {input} | cut -f1,30 -d "," | grep SRR > {output}
        """

# Save SRR ITS2 numbers to a list (optional)
# ------------------------------------------------
rule srr_its2_list:
    input:
        f"{resource_dir}/srr_runinfo.txt"
    output:
        f"{resource_dir}/srr_its2_list.txt"
    shell:
        """
        grep ITS2 {input} | cut -f1,30 -d "," | grep SRR > {output}
        """
