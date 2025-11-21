# 2b-RAD population genetics: *Stephanocoenia intersepta*
Kirsten Carlson
<br>Version: 11/19/2025
<br>
<br>This document walks through the process of downloading raw fastq files from NCBI and processing 2b-RADseq reads. 
<br>
## Setup
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
```
Lastly, most of the workflow is is processed by the Snakefile. Update the config to your preferences.
## Downloading raw fastqs from NCBI
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
## Quality control with fastqc and multiqc
To run quality control, use the rules ```fastq_raw``` and ```multiqc_raw``` in the Snakefile ```trim_removeSymb.snakefile```
```bash
snakemake -s trim_removeSymb.snakefile multiqc_raw
```
## Trimming and quality filtering
2bRAD sequencing has 2 barcodes, and is trimmed with the Matz lab perl script ```trim2bRAD_2barcodes_dedup.pl```. However, this appears to have already been done. It is included in ```trim_removeSymb.snakefile``` if needed, though.
<br>
<br>To run quality filtering with cutadapt use the rule ```quality_filter```. The default is quality="15,15" minlen=36, but you can edit with ```--config```.
```bash
snakemake -s trim_removeSymb.snakefile quality_filter
```
After quality filtering, run ```countreads.sh``` in the trimmed fastq directory.
```bash
cd /scratch/user/data/trimmed_fastq
bash countreads.sh

Total reads: 455539135
```
I got the same number of reads as the raw fastqs, so it appears trimmed fastqs may have been what was uploaded to NCBI. Since the number of reads are the same, skip the trimmed fastqc/multiqc step. However, these rules are included in trim_removeSymb.snakefile if needed.
## Remove contamination

