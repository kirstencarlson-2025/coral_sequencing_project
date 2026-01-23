# IN PROGRESS

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

# Snakemake rule pipeline: 
# create_symb_reference -> index_symbiont -> align_symbiont -> exclude_symbiont_reads -> unique_reads -> merge_unique_reads -> 
# filter_tags -> tab_to_fasta -> cluster_tags -> kraken2_filter -> copy_denovo_ref -> construct_denovo_ref -> index_denovo_ref -> cleanup_denovo_intermediate

#####################################################
# CONFIGURATION
######################################################

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
kraken_db_std = config["kraken_db_std"]
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

rule all:
    input:
        f"{denovo_ref_dir}/cleanup.done"



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

#######THESE DON'T WORK###############
# Build bowtie2 index for symbiont reference genome
rule index_symbiont_reference:
    input:
        fasta=f"{config['reference_symbiont']}.fasta",
        index=f"{config['reference_symbiont']}",
        conda="/home/kcarls36/.conda/envs/get_data"
    output:
        f"{config['reference_symbiont']}.1.bt2"
    shell:
        """
        bowtie2-build {input.fasta} {input.index}
        """
     
# Build bowtie2 index for symbiont reference genome
#rule index_symbiont:
#    input:
#        fasta=config['reference_symbiont'] + ".fasta"
#    output:
#        expand(f"{config['reference_symbiont']}.{{ext}}",
#               ext=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"])
#    conda:
#        config["bowtie2_env"]
#    run:
#        # strip .fasta from index prefix
#        index_prefix = input.fasta.rsplit('.',1)[0]
#        shell(f"module load mamba/latest \ source activate {conda} \ bowtie2-build {input.fasta} {index_prefix}")
#rule index_symbiont:
#    input:
#        fasta=f"{config['reference_symbiont']}.fasta"
#    output:
#        expand(f"{config['reference_symbiont']}.{{ext}}",
#               ext=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"])
#    shell:
#        f"bowtie2-build {input.fasta} {config['reference_symbiont']}"
####################################### 


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

##############################
# Create de novo reference
##############################
# Uniquing reads
rule unique_reads:
    input:
        f"{trimfq_noSymb_dir}/{{sample}}.fastq"
    output:
       f"{merge_dir}/{{sample}}.uni"
    threads: config["threads"]
    shell:
        """
        perl {scripts_dir}/uniquerOne.pl {input} > {output}
        """
# Merge unique reads
rule merge_unique_reads:
    input:
        expand(f"{merge_dir}/{{sample}}.uni", sample=SAMPLES)
    output:
        f"{merge_dir}/all.uniq"
    params:
        minInd = 30
    shell:
        """
        perl {scripts_dir}/mergeUniq.pl {input} minInd={params.minInd}  > {output}
        """
# Filter tags based on quality and count
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
rule cluster_tags: 
    input:
        f"{merge_dir}/all.fasta"
    output:
        f"{merge_dir}/cdh_alltags.fas"
    params:
        identity = 0.91
    shell:
        """
        cd-hit-est -i {input} -o {output} -c {params.identity} -aL 1 -aS 1 -g 1 -T 0 -M 0
        """

# Remove contamination with Kraken2 standard database
rule kraken2_filter:
    input:
        f"{merge_dir}/cdh_alltags.fas"
    output:
        f"{merge_dir}/cdh_alltags.unclass.fa"
    params:
        kraken_db_std = config["kraken_db_std"]
    shell:
        """
        kraken2 --db {params.kraken_db_std} --output /dev/null --report {merge_dir}/kraken_report.txt --unclassified-out {output} {input}
        """
    
# Copy and rename unclassified cdh tags as de novo reference genome
rule copy_denovo_ref:
    input:
        f"{merge_dir}/cdh_alltags.unclass.fa"
    output:
        f"{denovo_ref_dir}/sint_denovo.fa"
    shell:
        """
        cp {input} {output}
        """
# Construct de novo refernece with 30 psuedo chromosomes
rule construct_denovo_ref:
    input:
        f"{denovo_ref_dir}/sint_denovo.fa"
    output:
        fasta = f"{denovo_ref_dir}/sint_denovo_cc.fa",
        tab = f"{denovo_ref_dir}/sint_denovo_cc.tab"
    params:
        num_chr = 30,
        script = f"{scripts_dir}/concatFasta.pl"
    shell:
        """
        perl {params.script} fasta={input} num_chr={params.num_chr}
        """

# Format and index de novo reference genome with bwa and samtools
rule index_denovo_ref:
    input:
        f"{denovo_ref_dir}/sint_denovo_cc.fa"
    output:
        touch(f"{denovo_ref_dir}/sint_denovo_cc.fa.bwt")
    shell:
        """
        bwa index {input}
        samtools faidx {input}
        touch {output}
        """
# Final cleanup of de novo reference intermediate files
rule cleanup_denovo_intermediate:   
    input:
        fasta = f"{denovo_ref_dir}/sint_denovo_cc.fa",
        tab = f"{denovo_ref_dir}/sint_denovo_cc.tab",
        fq = f"{merge_dir}/cdh_alltags.fas"
    output:
        touch(f"{denovo_ref_dir}/cleanup.done")
    shell:
        """
        rm {merge_dir}/*.uni
        rm {input.fq}
        rm {merge_dir}/all.uniq
        rm {merge_dir}/all.tab
        rm {merge_dir}/all.fasta
        rm {input.tab}
        touch {output}
     """
