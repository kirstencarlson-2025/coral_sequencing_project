# 2b-RAD population genetics: *Stephanocoenia intersepta*
Kirsten Carlson
<br>Updated: 1/2026
<br>
<br>This workflow uses Snakemake to process 2bRAD sequencing reads from raw data to analysis-ready files. It includes downloading FASTQs from NCBI, trimming and quality control of reads, filtering symbiont sequences, and preparing data for population genetic analyses using INDELs as markers.
<br>
<details>
  <summary><strong>Setup</strong></summary>  
First, be sure snakemake is installed. I created a mamba environment to install/load tools to.
  
```bash
module load mamba/latest
mamba create -n snakemake
source activate snakemake
conda install bioconda::snakemake
```

The tools needed for the workflow are included in "/workflow/envs/env.yaml", which will install the first time you run snakemake.

Lastly, update the config file ("/config/config.yaml") to your naming and directory structure preferences.
</details>

<details>
  <summary><strong>Downloading raw fastqs </strong></summary>
In "/workflow/rules", this step includes "metadata.smk" and "download.smk". To begin, "metadata.smk" creates a text file including the SRR number and the SampleID number and saves it to "/resources". However, in this dataset there were three field triplicates that do not carry over the SampleID designation. I manually edited the text file using the NCBI database to update the SampleID numbers, and the file can be found in "resources/sample_rename.csv". You can use this to rename samples after downloading.
<br>
<br> If you'd like to estimate download size required, run:

```bash
bash download_size.sh srr_2brad_list.txt
```

The "download.smk" rule will download raw fastq files using the SRR numbers from the list, and name the files "[SampleID].fastq".

Check how many reads with ```countreads.sh``` by navigating to raw fastq directory and running in command line.

```bash
cd /scratch/user/data/raw_fastq
bash countreads.sh

Total reads: 455539135
```
</details>

<details>
  <summary><strong>Trimming and quality filtering</strong></summary>

## Quality control with fastqc and multiqc
To run quality control, use the rules ```fastq_raw``` and ```multiqc_raw``` in the Snakefile ``` trim_qc.smk```

```bash
snakemake -s trim_qc.smk multiqc_raw --cores 4
```
  
2bRAD sequencing has 2 barcodes, and can be trimmed with the Matz lab perl script ```trim2bRAD_2barcodes_dedup.pl```. However, this appears to have already been done. It is included in ``` trim_qc.smk``` if needed, though.
<br>
<br>To run quality filtering with cutadapt use the rule ```quality_filter```. The default is quality="15,15" minlen=36, but you can edit with ```--config```.
  
```bash
snakemake -s trim_qc.smk quality_filter --cores 4
```

After quality filtering, run ```countreads.sh``` in the trimmed fastq directory.

```bash
cd /scratch/user/data/trimmed_fastq --cores 4
bash countreads.sh

Total reads: 455539135
```

If needed, run fastqc and multiqc on trimmed, quality-filtered reads.

```bash
snakemake -s trim_qc.smk multiqc_trimmed --cores 4
```
</details>

<details>
  <summary><strong>De novo reference</strong></summary>

### Remove contamination
First, you'll need to download Symbiodiniaceae genomes:

```bash
cd /scratch/user/data/symbgenomes
wget http://symbs.reefgenomics.org/download/SymbC1.Genome.Scaffolds.fasta.gz
wget http://smic.reefgenomics.org/download/Smic.genome.scaffold.final.fa.gz
wget https://marinegenomics.oist.jp/symbd/download/102_symbd_genome_scaffold.fa.gz
esearch -db assembly -query "GCA_000507305.1" | elink -target nuccore | efetch -format fasta > Breviolum_minutum.v1.0.genome.fa

# Rename 2 for clarity:
mv SymbC1.Genome.Scaffolds.fasta.gz Cladocopium_goreaui_Genome.Scaffolds.fasta.gz
mv Smic.genome.scaffold.final.fa.gz Symbiodinium_microadriacticum_genome.scaffold.fasta.gz
```

Next, create a concatenated genome from all references:

```bash
snakemake -s trim_qc.smk create_symb_reference --cores 4
```

Index the concatenated reference:

```bash
snakemake -s  trim_qc.smk index_symbiont --cores 4
```

Align reads to reference:

```bash
snakemake -s  trim_qc.smk align_symbiont --cores 4
```
Exclude symbiont reads from trimmed fastq files:
```bash
snakemake -s  trim_qc.smk exclude_symbiont_reads --cores 4
```

