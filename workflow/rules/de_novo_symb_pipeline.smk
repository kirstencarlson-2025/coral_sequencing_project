# ------------------------------------------------ #
# 2bRAD De Novo Reference Construction Snakefile
# Kirsten Carlson
# Updated 1/2026
# ------------------------------------------------ #

# This Snakefile runs the pipeline for processing 2bRAD data from Stephanocoenia intersepta samples.
# This script aligns reads to concatenated symbiont genomes to remove symbiont reads. A standard kraken2 database is used 
# to remove any remaining contaminant reads. Finally, reads are constructed into a de novo reference genome.

# Scripts needed for this pipeline:
# From https://github.com/z0on/2bRAD_denovo
#  - 2bRAD_trim_launch_dedup.pl
#  - trim2bRAD_2barcodes_dedup.pl
#  - uniquerOne.pl
#  - mergeUniq.pl
#  - concatFasta.pl

# ------------------------------------------------ #
# CONFIGURATION
# ------------------------------------------------ #

configfile: "../config/config.yaml"
# Conda environment
env = config["env"]
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
kraken_db_std = config["kraken_db_std"]
# De novo reference genome directory
denovo_ref_dir = config["denovo_ref_dir"]
# Script directory
scripts_dir = config["scripts_dir"] 
# Symbiont genomes
symbiont_genomes = config["symbiont_genomes"]
# Symbiont concatenated reference genome
reference_symbiont = config["reference_symbiont"] 
# De novo reference genome
denovo_ref_basename = config["denovo_ref_basename"]   

# Get all FASTQ files and extract sample names
SAMPLES=glob_wildcards(f"{rawfq_dir}/{{sample}}.fastq").sample

# Which rules to run locally
local_rules = ["create_symb_reference",
                "filter_tags",
                "tab_to_fasta",
                "rename_denovo_ref",
                "construct_denovo_ref",
                "cleanup_denovo_intermediate"]

# ------------------------------------------------ #
# Rules
# ------------------------------------------------ #

# Create concatenated symbiont reference genome from individual symbiont genomes
# ------------------------------------------------ 
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
# ------------------------------------------------ 
rule index_symbiont_reference:
    input:
        fasta=f"{config['reference_symbiont']}.fasta"
    output:
        f"{config['reference_symbiont']}.1.bt2",
        f"{config['reference_symbiont']}.2.bt2",
        f"{config['reference_symbiont']}.3.bt2",
        f"{config['reference_symbiont']}.4.bt2",
        f"{config['reference_symbiont']}.rev.1.bt2",
        f"{config['reference_symbiont']}.rev.2.bt2"
    threads: 4
    resources:
        mem_mb=32000
    params:
        basename=config['reference_symbiont']
    conda:
        config["env"]
    shell:
        """
        bowtie2-build --threads {threads} {input.fasta} {params.basename}
        """


# Align reads to symbiont reference genome
# ------------------------------------------------ 
rule map_symbiont:
    input:
        fq = f"{trimfq_dir}/{{sample}}.fastq",
        index = f"{config["reference_symbiont"]}.1.bt2"
    output:
        sam = f"{symb_align_dir}/{{sample}}_symb.sam",
        aligned = f"{symb_align_dir}/{{sample}}_symb.fastq",
        unaligned = f"{trimfq_noSymb_dir}/{{sample}}.fastq"
    params:
        index_prefix = {config['reference_symbiont']}
    threads: 4
    resources:
        mem_mb=16000
    conda:
        config["env"]
    shell:
        """
        bowtie2 \
          --score-min L,16,1 --local -L 16 \
          -x {params.index_prefix} \
            -U {input.fq} \
            --no-unal \
            --al {output.aligned} \
            --un {output.unaligned} \
            -S {output.sam} \
            -p {threads}
        """

# Uniquing reads
# ------------------------------------------------ 
rule unique_reads:
    input:
        f"{trimfq_noSymb_dir}/{{sample}}.fastq"
    output:
       f"{merge_dir}/{{sample}}.uni"
    conda:
        config["env"]
    threads: 1
    resources:
        mem_mb=16000
    shell:
        """
       {scripts_dir}/uniquerOne.pl {input} > {output}
        """

# Merge unique reads
# ------------------------------------------------ 
rule merge_unique_reads:
    input:
        expand(f"{merge_dir}/{{sample}}.uni", sample=SAMPLES)
    output:
        f"{merge_dir}/all.uniq"
    params:
        suffix= "uni",
        minInd = 30
    threads: 1
    resources:
        mem_mb=200000
    shell:
        """
        cd {merge_dir}
        {scripts_dir}/mergeUniq.pl {params.suffix} minInd={params.minInd}  > all.uniq
        """

# Filter tags based on quality and count
# ------------------------------------------------ 
rule filter_tags:
    input:
        f"{merge_dir}/all.uniq"
    output:
        f"{merge_dir}/all.tab"
    shell:
        """
        awk '!($3>7 && $4==0) && $2!="seq"' {input} > {output}
        """

# Convert filtered tags to fasta format
# ------------------------------------------------ 
rule tab_to_fasta:
    input:
        f"{merge_dir}/all.tab"
    output:
        f"{merge_dir}/all.fasta"
    shell:
        """
        awk '{{print ">"$1"\\n"$2}}' {input} > {output}
        """

# Cluster filtered tags with CD-HIT
# ------------------------------------------------ 
rule cluster_tags: 
    input:
        f"{merge_dir}/all.fasta"
    output:
        f"{merge_dir}/cdh_alltags.fas"
    params:
        identity=0.91
    resources:
        mem_mb=200000
    shell:
        """
        cd-hit-est -i {input} -o {output} -c {params.identity} -aL 1 -aS 1 -g 1 -T 0 -M 0
        """

# Remove contamination with Kraken2 standard database
# ------------------------------------------------ 
rule kraken2_filter:
    input:
        f"{merge_dir}/cdh_alltags.fas"
    output:
        unclass = f"{merge_dir}/cdh_alltags.unclass.fa",
        report = f"{merge_dir}/kraken_report.txt"
    params:
        kraken_db_std = config["kraken_db_std"]
    threads: 8
    resources:
        mem_mb=200000
    conda:
        config["env"]
    shell:
        """
        kraken2 --threads {threads} --db {params.kraken_db_std} --output /dev/null --report {output.report} --unclassified-out {output.unclass} {input} 
        """

# Copy and rename de novo reference genome
# ------------------------------------------------
rule rename_denovo_ref:
    input:
        f"{merge_dir}/cdh_alltags.unclass.fa"
    output:
        f"{denovo_ref_dir}/{denovo_ref_basename}.fa"
    shell:
        """
        mv {input} {output}
        """ 

# Construct de novo reference with 30 psuedo chromosomes
# ------------------------------------------------ 
rule construct_denovo_ref:
    input:
        f"{denovo_ref_dir}/{denovo_ref_basename}.fa"
    output:
        fasta = f"{denovo_ref_dir}/{denovo_ref_basename}_cc.fasta",
        tab = f"{denovo_ref_dir}/{denovo_ref_basename}_cc.tab"
    params:
        num_chr=30
    shell:
        """
        perl {scripts_dir}/concatFasta.pl fasta={input} num={params.num_chr}
        """

# Format and index de novo reference genome
# ------------------------------------------------ 
rule index_denovo_ref:
    input:
        f"{denovo_ref_dir}/{denovo_ref_basename}_cc.fasta"
    output:
        f"{denovo_ref_dir}/{denovo_ref_basename}_cc.1.bt2",
        f"{denovo_ref_dir}/{denovo_ref_basename}_cc.2.bt2",
        f"{denovo_ref_dir}/{denovo_ref_basename}_cc.3.bt2",
        f"{denovo_ref_dir}/{denovo_ref_basename}_cc.4.bt2",
        f"{denovo_ref_dir}/{denovo_ref_basename}_cc.rev.1.bt2",
        f"{denovo_ref_dir}/{denovo_ref_basename}_cc.rev.2.bt2",
        f"{denovo_ref_dir}/{denovo_ref_basename}_cc.fasta.fai"
    threads: 4
    resources:
        mem_mb=32000
    conda:
        config["env"]
    params:
        basename=f"{denovo_ref_dir}/{denovo_ref_basename}_cc"
    shell:
        """
        bowtie2-build --threads {threads} {input} {params.basename}
        samtools faidx {input}
        """
        

# Final cleanup of de novo reference intermediate files
# ------------------------------------------------ 
rule cleanup_denovo_intermediate:   
    input:
        fasta = f"{denovo_ref_dir}/{denovo_ref_basename}_cc.fasta",
        tab = f"{denovo_ref_dir}/{denovo_ref_basename}_cc.tab",
        index = f"{denovo_ref_dir}/{denovo_ref_basename}_cc.1.bt2"
    output:
        touch(f"{denovo_ref_dir}/cleanup.done")
    shell:
        """
        rm {merge_dir}/*.uni
        rm {merge_dir}/all.uniq
        rm {merge_dir}/all.tab
        rm {merge_dir}/all.fasta
        rm {input.tab}
        touch {output}
     """
