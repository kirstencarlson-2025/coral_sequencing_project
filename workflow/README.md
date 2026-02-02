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
In "/workflow/rules", this step includes "metadata.smk" and "download.smk". To begin, "metadata.smk" creates a text file including the SRR number and the SampleID number and saves it to "/resources". However, in this dataset there were three field triplicates that do not carry over the SampleID designation. I manually edited the text file using the NCBI database to update the SampleID numbers, and the file can be found in "resources/sample_rename.csv". You can edit the snakefile rule to use this list to download and name files.
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
The rule ```trim_qc.smk``` will run quality control. First, ```fastq_raw``` and ```multiqc_raw``` will create a quality control summary of the raw fastqs.
<br>
<br>2bRAD sequencing has 2 barcodes, and can be trimmed with the Matz lab perl script ```trim2bRAD_2barcodes_dedup.pl``` using the snakemake rule ```trim_dedup```. However, this appears to have already been done. It is included if needed, though.
<br>
<br>The rule ```quality_filter``` runs quality filtering with ```cutadapt```. The default is quality="15,15" minlen=36, which you can edit in the config file.
<br>
<br>After quality filtering, run ```countreads.sh``` in the trimmed fastq directory.

```bash
cd /scratch/user/data/trimmed_fastq
bash countreads.sh

Total reads: 455539135
```
Lastly, ```trim_qc.smk``` runs quality control again on the trimmed and quality filtered fastqs.
</details>

<details>
  <summary><strong>De novo reference</strong></summary>

### Remove contamination
First, you'll need to download Symbiodiniaceae genomes:

```bash
cd /scratch/user/data/symbgenomes
wget http://symbs.reefgenomics.org/download/SymbC1.Genome.Scaffolds.fasta.gz > Cladocopium_goreaui_Genome.Scaffolds.fasta.gz
wget http://smic.reefgenomics.org/download/Smic.genome.scaffold.final.fa.gz > Symbiodinium_microadriacticum_genome.scaffold.fasta.gz
wget https://marinegenomics.oist.jp/symbd/download/102_symbd_genome_scaffold.fa.gz
esearch -db assembly -query "GCA_000507305.1" | elink -target nuccore | efetch -format fasta > Breviolum_minutum.v1.0.genome.fa
```

The rule ```de_novo_symb_pipeline.smk``` aligns reads to the concatenated symbiont genomes to remove any symbiont reads, then a standard ```kraken2``` database removes any remaining contamination. Finally, the reads are constructed into a de novo reference genome.
<br>
<br>Rules in this snakefile include:
<br>```create_symb_reference``` (concatenates symbiont genomes)
<br>```index_symbiont_reference``` (indexes the concatenated genome)
<br>```map_symbiont``` (aligns read to symbiont reference)
<br>```exlude_symbiont_reads``` (excludes any reads that align to the symbiont reference)
<br>```unique_reads``` (uses uniquerOne.pl to unique reads)
<br>```merge_unique_reads``` (merges the uniqued reads)
<br>```filter_tags``` (filters tags based on quality and count)
<br>```tab_to_fasta``` (converts filtered tags to fasta format)
<br>```cluster_tags``` (cluster the filtered tags with ```CD-HIT```)
<br>```kraken2_filter``` (removes contamination with ```kraken2``` standard database)
<br>```construct_denovo_ref``` (construct de novo reference with 30 psuedo chromosomes)
<br>```index_denovo_ref``` (index de novo reference)
<br>```cleanup_denovo_intermediates``` (final cleanup of de novo intermediate files)

