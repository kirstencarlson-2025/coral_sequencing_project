# ------------------------------------------------ #
# 2bRAD DiscoSNP-RAD Snakefile
# Kirsten Carlson
# Updated 5/2026
# ------------------------------------------------ #

# This Snakefile runs the pipeline for processing filtered 2bRAD fastq files from Stephanocoenia intersepta samples with DiscoSNP to call variants.
# It can be used to perform a parameter sweep of kmer length and max deletion size.
# It also creates a slimmed vcf for exploratory QC with SeqArray in R, then filters variants by cluster size, rank, coverage, missing data, minor allele frequency, and paralogs.
# Variant reports are created for parameter sweep of kmer length and max deletion size, as well as for filtering steps.


# ------------------------------------------------ #
# CONFIGURATION
# ------------------------------------------------ #

# CONFIGURATION
configfile: "../config/config.yaml"
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
# Results directory
results_dir = config["results_dir"]
# De novo reference directory
denovo_ref_dir = config["denovo_ref_dir"]
# De novo reference base name
denovo_ref_basename = config["denovo_ref_basename"]
# Scripts directory
scripts_dir = config["scripts_dir"]
# Parameter sweep for discoSnp_Rad
KMERS = config["discosnp_kmers"]
DELS = config["discosnp_dels"]
disco_percent_heterozygotes = config["disco_percent_heterozygotes"]
disco_percent_variants = config["disco_percent_variants"]


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
        vcf = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_clustered.vcf.gz",
        tbi = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_clustered.vcf.gz.tbi"
    threads: 12
    conda:
        config["env"]
    resources:
        mem_mb=50000,
        runtime=1440,
        cpus_per_task=12
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
        -t {threads}

        # Compress and index
        bgzip -c discoRad_k_{wildcards.k}_c_3_D_{wildcards.D}_P_5_m_5_clustered.vcf > discoRad_k_{wildcards.k}_c_3_D_{wildcards.D}_P_5_m_5_clustered.vcf.gz
        tabix -p vcf discoRad_k_{wildcards.k}_c_3_D_{wildcards.D}_P_5_m_5_clustered.vcf.gz
        """

### Steps to map discoSnp_Rad output with a reference
# Create VCF
# ------------------------------------------------
rule create_vcf:
    input:
        fa = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_raw_filtered.fa",
        ref = f"{denovo_ref_dir}/{denovo_ref_basename}_cc.fasta"
    output:
        temp_vcf = temp(f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/add_cluster_info_temp_1.vcf")
    conda:
        config["env"]
    resources:
        mem_mb=128000,
        runtime=1440,
        cpus_per_task=12
    shell:
        """
        outdir=$(dirname {output.temp_vcf})
        mkdir -p $outdir

        cd $outdir

        $CONDA_PREFIX/scripts/run_VCF_creator.sh \
        -G {input.ref} \
        -p {input.fa} \
        -e \
        -o {output.temp_vcf}
        """

# Convert VCF from 0 to 1 format
# ------------------------------------------------
rule convert_vcf_format:
    input:
        temp_vcf = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/add_cluster_info_temp_1.vcf"
    output:
        temp_vcf = temp(f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/add_cluster_info_temp_2.vcf")
    conda:
        config["env"]
    shell:
        """
        python {scripts_dir}/zero2one.py -i {input.temp_vcf} -o {output.temp_vcf}
        """

# Add cluster info to mapped VCF
# ------------------------------------------------
rule add_cluster_info:
    input:
        temp_vcf = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/add_cluster_info_temp_2.vcf",
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
        -m {input.temp_vcf} \
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
        bgzip -f {input.vcf}
        tabix -f -p vcf {output.vcf}
        """
###

# Create variant report before filtering
# ------------------------------------------------
rule create_variant_report_before_filtering:
    input:
        clustered = expand(f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_clustered.vcf.gz", k=KMERS, D=DELS),
        mapped = expand(f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_mapped.vcf.gz", k=KMERS, D=DELS)
    output:
        report = f"{results_dir}/variant_report_before_filtering.txt"
    conda:
        config["env"]
    threads: 1
    run:
        import re   
        
        with open(output.report, "w") as out:
            out.write("k\tD\tall_variants\tsnps\tindels\n")


            # Map labels to files listed
            file_groups = {
                "clustered": input.clustered,
                "mapped": input.mapped
            }
            
            for ftype, vcfs in file_groups.items():
                for vcf in vcfs:
                    # Extract k and D
                    m = re.search(r"k(\d+)_D(\d+)", vcf)
                    k_val, D_val = m.groups()

                # Count variants
                all_count = int(shell(f"bcftools view -H {vcf} | wc -l", read=True).strip())
                snp_count = int(shell(f"bcftools view -H -v snps {vcf} | wc -l", read=True).strip())
                indel_count = int(shell(f"bcftools view -H -v indels {vcf} | wc -l", read=True).strip())

                # Write output
                out.write(f"{k_val}\t{D_val}\t{all_count}\t{snp_count}\t{indel_count}\n")


# Prepare VCF for exploratory QC with SeqArray in R
# Sort the clustered VCF
# ------------------------------------------------
rule sort_clustered_vcf:
    input:
        clustered = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_clustered.vcf.gz",
        mapped = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_mapped.vcf.gz"
    output:
        clustered = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_sorted_clustered.vcf.gz",
        mapped = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_sorted_mapped.vcf.gz"
    conda:
        config["env"]
    shell:
        """
        bcftools sort {input.clustered} -Oz -o {output.clustered}
        bcftools sort {input.mapped} -Oz -o {output.mapped}
        """

# Reheader the VCF with correct sampleID (discoSnp_Rad ouputs sampleID as G1-Gx in the order of the fof.txt)
# Create sample map from fof.txt
# ------------------------------------------------
rule create_sample_map:
    input:
        fof = f"{sint_align_dir}/discosnp/fof.txt"
    output:
        sample_map = f"{sint_align_dir}/discosnp/sample_map.txt"
    conda:
        config["env"]
    run:
        """
        awk '{split($1,a,"/"); sub(".fastq$","",a[length(a)]); print "G"NR"\t"a[length(a)]}' {input.fof} > {output.sample_map}
        """

# Reheader VCF
# ------------------------------------------------
rule reheader_vcf:
    input:
        clustered = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_sorted_clustered.vcf.gz",
        mapped = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_sorted_mapped.vcf.gz",
        sample_map = f"{sint_align_dir}/discosnp/sample_map.txt"
    output:
        clustered = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_sorted_reheader_clustered.vcf.gz",
        mapped = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_sorted_reheader_mapped.vcf.gz"
    conda:
        config["env"]
    shell:
        """
        bcftools reheader -s {input.sample_map} {input.clustered} -o {output.clustered}
        bcftools reheader -s {input.sample_map} {input.mapped} -o {output.mapped}
        """

# Slim up the sorted and reheadered clustered VCF for exploratory QC with SeqArray in R
# ------------------------------------------------
rule slim_clustered_vcf:
    input:
        vcf = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_sorted_reheader_clustered.vcf.gz"
    output:
        vcf = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_D_{{D}}_slim.vcf.gz",
        tbi = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_D_{{D}}_slim.vcf.gz.tbi"
    conda:
        config["env"]
    shell:
        """
        bcftools annotate -x INFO,FORMAT,CONTIG {input.vcf} -Oz -o {output.vcf}
        tabix -p vcf {output.vcf}
        """

# Great! You can now use the Rmd file in the workflow/code directory to perform exploratory QC with SeqArray in R.

# ------------------------------------------------ #
# Filtering steps
# ------------------------------------------------ #

# Download discoSnp_Rad filtering scripts from GitHub
# ------------------------------------------------
rule download_discosnp_filtering_scripts:
    output:
        filter_csr = f"{scripts_dir}/filter_by_cluster_size_and_rank.py",
        filter_cmismaf = f"{scripts_dir}/filter_vcf_by_indiv_cov_max_missing_and_maf.py",
        filter_paralogs = f"{scripts_dir}/filter_paralogs.py"
    shell:
        """
        wget -O {output.filter_csr} https://github.com/GATB/DiscoSnp/blob/master/discoSnpRAD/post-processing_scripts/filter_by_cluster_size_and_rank.py
        wget -O {output.filter_cmismaf} https://github.com/GATB/DiscoSnp/blob/master/discoSnpRAD/post-processing_scripts/filter_vcf_by_indiv_cov_max_missing_and_maf.py
        wget -O {output.filter_paralogs} https://github.com/GATB/DiscoSnp/blob/master/discoSnpRAD/post-processing_scripts/filter_paralogs.py
        chmod +x {output.filter_csr} {output.filter_cmismaf} {output.filter_paralogs}
        """

# Filter by cluster size and rank
# ------------------------------------------------
rule filter_cluster_size_rank:
    input:
        vcf = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_sorted_reheader_mapped.vcf.gz",
        script = f"{scripts_dir}/filter_by_cluster_size_and_rank.py"
    output:
        temp_vcf = temp(f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/temp.vcf"),
        vcf = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_filter_csr.vcf",
        zip = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_filter_csr.vcf.gz"
    params:
        env = config["env"],
        min_cluster_size = 3,
        Max_cluster_size = 100,
        rank = 0.4
    threads: 4
    resources:
        mem_mb=200000
    shell:
        """
        gunzip -c {input.vcf} > {output.temp_vcf}

        python {input.script} -i {output.temp_vcf} -o {output.vcf} \
        -m {params.min_cluster_size} \
        -M {params.Max_cluster_size} \
        -r {params.rank}

        bgzip -c {output.vcf} > {output.zip}
        tabix -p vcf {output.zip}
        """

# Filter by coverage, missing genotypes, and minor allele frequency
# ------------------------------------------------
rule filter_coverage_missing_maf:
    input:
        vcf = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_filter_csr.vcf.gz",
        script = f"{scripts_dir}/filter_vcf_by_indiv_cov_max_missing_and_maf.py"
    output:
        temp_vcf = temp(f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/filter_temp_1.vcf"),
        vcf = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_filter_cmismaf.vcf",
        zip = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_filter_cmismaf.vcf.gz"
    params:
        env = config["env"],
        min_cov = 3, # From Eckert et al. 2024, ANGSD -setMinDepthInd 3, which is the minimum coverage per individual to call a genotype. 
                     # This is not the same as the minimum depth to call a variant, but it is a similar concept and provides a starting point for filtering. 
                     # We can adjust this threshold based on the distribution of coverage in our data and the results of downstream analyses.
        max_missing = 0.25, # From Eckert et al. 2024, ANGSD -minInd 165, which is the minimum number of individuals that must have a called genotype for a variant to be retained. 
                             # This corresponds to a maximum missing data threshold of ((220-165)/220) = 0.25. 
                             # We can adjust this threshold based on the distribution of missing data in our variants and the results of downstream analyses.
        min_maf = 0.05 # From Eckert et al. 2024, ANGSD -setMinMaf 0.05, which is the minimum minor allele frequency to retain a variant.
    threads: 4
    resources:
        mem_mb=200000
    shell:
        """
        gunzip -c {input.vcf} > {output.temp_vcf}
        
        python {input.script} -i {output.temp_vcf} -o {output.vcf} \
        -c {params.min_cov} \
        -m {params.max_missing} \
        -f {params.min_maf} 

        bgzip -c {output.vcf} > {output.zip}
        tabix -p vcf {output.zip}
        """
    

# Create variant report for first two filtering steps
# ------------------------------------------------
rule create_variant_by_filter_report:
    input:
        original = expand(f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_sorted_reheader_mapped.vcf.gz", k=KMERS, D=DELS),
        csr = expand(f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_filter_csr.vcf.gz", k=KMERS, D=DELS),
        cmismaf = expand(f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_filter_cmismaf.vcf.gz", k=KMERS, D=DELS)
    output:
        report = f"{results_dir}/variant_by_filter_report.txt"
    conda:
        config["env"]
    threads: 1
    run:
        import re

        # Map labels to file lists
        file_groups = {
            "original": input.original,
            "csr": input.csr,
            "cmismaf": input.cmismaf
        }

        with open(output.report, "w") as out:
            out.write("file\tk\tD\tall_variants\tsnps\tindels\n")

            for ftype, vcfs in file_groups.items():
                for vcf in vcfs:
                    # Extract k and D
                    m = re.search(r"k(\d+)_D(\d+)", vcf)
                    k_val, D_val = m.groups()

                    # Count variants
                    all_count = int(shell(f"bcftools view -H {vcf} | wc -l", read=True).strip())
                    snp_count = int(shell(f"bcftools view -H -v snps {vcf} | wc -l", read=True).strip())
                    indel_count = int(shell(f"bcftools view -H -v indels {vcf} | wc -l", read=True).strip())

                    # Write output
                    out.write(f"{ftype}\t{k_val}\t{D_val}\t{all_count}\t{snp_count}\t{indel_count}\n")


# Filter by paralogs
# ------------------------------------------------
rule filter_paralogs:
    input:
        vcf = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_filter_cmismaf.vcf.gz",
        script = f"{scripts_dir}/filter_paralogs.py"
    output:
        temp_vcf = temp(f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/filter_temp_2_hetero{{hetero}}_variants{{variants}}.vcf"),
        vcf = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_filtered_hetero{{hetero}}_variants{{variants}}.vcf",
        zip = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_filtered_hetero{{hetero}}_variants{{variants}}.vcf.gz",
        tbi = f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_filtered_hetero{{hetero}}_variants{{variants}}.vcf.gz.tbi"
    conda:
        config["env"]
        # Filter so that more than 50% of variants have more than 10% of heterozygous genotypes.
        # Example parameters given in DiscoSnp_Rad COOKBOOK. We can adjust these thresholds based on the distribution of heterozygous genotypes.
    threads: 4
    resources:
        mem_mb=200000
    shell:
        """
        gunzip -c {input.vcf} > {output.temp_vcf}
        
        python {input.script} -i {output.temp_vcf} -o {output.vcf} \
        -x {wildcards.hetero} \
        -y {wildcards.variants}
        
        bgzip -c {output.vcf} > {output.zip}
        tabix -p vcf {output.zip}
        """

# Create variant report after all filtering steps
# ------------------------------------------------
rule create_paralog_variant_report:
    input:
        expand(f"{sint_align_dir}/discosnp/k{{k}}_D{{D}}/discoRad_k_{{k}}_c_3_D_{{D}}_P_5_m_5_filtered_hetero{{hetero}}_variants{{variants}}.vcf.gz", 
        k=KMERS, 
        D=DELS, 
        hetero=disco_percent_heterozygotes, 
        variants=disco_percent_variants)
    output:     
        report = f"{results_dir}/paralog_variant_report.txt"
    conda:
        config["env"]
    threads: 1
    run:
        import re 
        import os  

        with open(output.report, "w") as out:
            out.write("k\tD\thetero\tvariants\tsnps\tindels\n")

            for vcf in input:
                fname = os.path.basename(vcf)
                parts = fname.replace(".vcf.gz", "").split("_")

                k_val = parts[parts.index("k") + 1]
                D_val = parts[parts.index("D") + 1]
                hetero_val = [p for p in parts if p.startswith("hetero")][0].replace("hetero", "")
                variants_val = [p for p in parts if p.startswith("variants")][0].replace("variants", "")

                # Count variants
                all_count = int(shell(f"bcftools view -H {vcf} | wc -l", read=True).strip())
                snp_count = int(shell(f"bcftools view -H -v snps {vcf} | wc -l", read=True).strip())
                indel_count = int(shell(f"bcftools view -H -v indels {vcf} | wc -l", read=True).strip())

                # Write output
                out.write(f"{k_val}\t{D_val}\t{hetero_val}\t{variants_val}\t{snp_count}\t{indel_count}\n")


