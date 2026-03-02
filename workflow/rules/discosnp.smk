# IN PROGRESS
# ------------------------------------------------ #
# 2bRAD DiscoSNP-RAD Parameter Sweep Snakefile
# Kirsten Carlson
# Updated 1/2026
# ------------------------------------------------ #

# This Snakefile runs the pipeline for processing filtered 2bRAD fastq files from Stephanocoenia intersepta samples with DiscoSNP to call variants.
# It also filters and processes vcf files to create a final set of variants for downstream analyses.


# ------------------------------------------------ #
# CONFIGURATION
# ------------------------------------------------ #

# CONFIGURATION
configfile: "../../config/config.yaml"
# Conda environment
env = config["env"]
# Raw fastq file directory
rawfq_dir = config["rawfq_dir"]
# Excluded symbiont trimmed fastq directory
trimfq_noSymb_dir = config["trimfq_noSymb_dir"]
# Sint alignment directory
sint_align_dir = config["sint_align_dir"]
# Resource directory (SRR numbers, run info)
resource_dir = config["resource_dir"]
# De novo reference directory
denovo_ref_dir = config["denovo_ref_dir"]
# De novo reference base name
denovo_ref_basename = config["denovo_ref_basename"]
# Scripts directory
scripts_dir = config["scripts_dir"]


# Get all FASTQ files and extract sample names
SAMPLES=glob_wildcards(f"{rawfq_dir}/{{sample}}.fastq").sample

# Parameter sweep values
KMERS = [25] # kmer lengths to test
DELS = [5] # max deletion size

# Define final target(s)
# Temporary rule_all for debugging
rule all:
    input:
        expand(
            f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_mapped.vcf.gz",
            k=KMERS,
            D=DELS
        )

# ------------------------------------------------ #
# Rules
# ------------------------------------------------ #

# Create file-of-files
# ------------------------------------------------
rule create_file_of_files:
    input:
        expand(f"{trimfq_noSymb_dir}/{{sample}}.fastq", sample=SAMPLES)
    output:
        f"{sint_align_dir}/discosnp/fof.txt"
    shell:
        "ls {input} > {output}"

# DiscoSNP align and call variants
# ------------------------------------------------
rule run_discosnpRad:
    input:
        fof = f"{sint_align_dir}/discosnp/fof.txt"
    output:
        vcf = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_clustered.vcf.gz"
    params:
        env = config["env"]
    threads: 12
    conda:
        config["env"]
    resources:
        mem_mb=50000,
        runtime=1440
    shell:
        """
        # Create output directory
        mkdir -p {sint_align_dir}/discosnp/k{wildcards.k}_D{wildcards.D}

        cd {sint_align_dir}/discosnp/k{wildcards.k}_D{wildcards.D}

        # Run DiscoSnp-RAD
        run_discoSnpRad.sh \
        -r "{input.fof}" \
        -S "$CONDA_PREFIX/bin" \
        -k {wildcards.k} \
        -D {wildcards.D} \

        # Compress and index
        bgzip -c discoRad_k_{wildcards.k}_c_3_D_{wildcards.D}_P_5_m_5_clustered.vcf > discoRad_k_{wildcards.k}_c_3_D_{wildcards.D}_P_5_m_5_clustered.vcf.gz
        tabix -p vcf discoRad_k_{wildcards.k}_c_3_D_{wildcards.D}_P_5_m_5_clustered.vcf.gz
        """

# Create VCF
# ------------------------------------------------
rule create_vcf:
    input:
        fa = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_raw_filtered.fa",
        ref = f"{denovo_ref_dir}/{denovo_ref_basename}_cc.fasta"
    output:
        vcf = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/temp.vcf"
    conda:
        config["env"]
    resources:
        mem_mb=50000,
        runtime=600
    shell:
        """
        $CONDA_PREFIX/scripts/run_VCF_creator.sh \
        -G {input.ref} \
        -p {input.fa} \
        -e \
        -o {output.vcf}
        """

# Convert VCF from 0 to 1 format
# ------------------------------------------------
rule convert_vcf_format:
    input:
        vcf = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/temp.vcf"
    output:
        vcf = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/temp_1.vcf"
    conda:
        config["env"]
    shell:
        """
        python {scripts_dir}/zero2one.py -i {input.vcf} -o {output.vcf}
        """


# Add cluster info to mapped VCF
# ------------------------------------------------
rule add_cluster_info:
    input:
        vcf = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/temp_1.vcf",
        vcf_clustered = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_clustered.vcf"
    output:
        vcf = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_mapped.vcf"
    resources:
        mem_mb=50000,
        runtime=1440
    conda:
        config["env"]
    shell:
        """
        python $CONDA_PREFIX/discoSnpRAD/post-processing_scripts/add_cluster_info_to_mapped_vcf.py \
        -m {input.vcf} \
        -u {input.vcf_clustered} \
        -o {output.vcf}
        """

# Compress and index final mapped VCF
# ------------------------------------------------
rule compress_index_vcf:
    input:
        vcf = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_mapped.vcf"
    output:
        vcf = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_mapped.vcf.gz"
    conda:
        config["env"]
    shell:
        """
        bcftools sort {input.vcf} -Oz -o {output.vcf}
        tabix -p vcf {output.vcf}
        """
