# 2b-RAD population genetics: *Stephanocoenia intersepta*
Kirsten Carlson
<br>Updated: 11/19/2025
<br>
<br>This document walks through the process of downloading raw fastq files from NCBI and processing 2b-RADseq reads. 
<br>
<details>
  <summary><strong>Setup</strong></summary>  
I created a mamba environment to install/load tools to.
  
```bash
module load mamba/latest
mamba create -n sint
source activate sint
```

<br> These are the tools you'll need:
```bash
# This needs to be edited
python
angsd
cd-hit
entrez-direct
htslib
multiqc
ngsrelate
parallel-fastq-dump
perl?
pysradb
snakemake
bowtie2
```
Lastly, most of the workflow is is processed by the Snakefile. Update the config to your preferences.
</details>

<details>
  <summary><strong>Downloading raw fastqs </strong></summary>
First, this script uses the PRJNA accession number to get all SRR numbers and sampleIDs and save them to a csv file (one per line). SRR numbers are split into two files, one for 2bRAD fastqs and one for ITS2 fastqs. SampleIDs will be used for naming to match metadata.
  
```bash
bash get_SRR.sh $PRJNA# $OutputDirectory $OutputfileName
```

Example:

```bash
bash get_SRR.sh PRJNA884416 /scratch/user/data/srr_numbers SRR_sampleID.csv
```

Next, you can estimate download size required with vdb-dump

```bash
bash download_size.sh 2bRAD_SRR_sampleID.csv
```

Finally, use the sbatch script to schedule the download via SLURM scheduler. By default, it loads mamba, but if you are using conda, be sure to edit the script. The script ```download_fastq.sh``` uses SRR numbers to download fastq files. Raw fastq files are named "[sampleID].fastq".

```bash
sbatch download_fastq.sbatch $mamba/conda_EnvironmentName $OutputDirectory $SRR_list
```

Example:

```bash
sbatch --mail-user=user@example.com download_fastq.sbatch download /scratch/user/data/raw_fastq 2bRAD_SRR_sampleID.csv
```

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

