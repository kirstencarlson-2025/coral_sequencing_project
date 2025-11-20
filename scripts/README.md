# 2b-RAD population genetics: *Stephanocoenia intersepta*
Kirsten Carlson
<br>Version: 11/19/2025
<br>
<br>This document walks through the process of downloading raw fastq files from NCBI and processing
<br>processing 2b-RADseq reads. 
<br>
## Tools needed
I created a mamba environment to install/load these tools on.
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
## Downloading raw fastqs from NCBI
First, this script uses the PRJNA accession number to get all SRR numbers and sampleIDs and save them to a csv file (one per line). SRR numbers are split into two files, one for 2bRAD fastqs and one for ITS2 fastqs. SampleIDs will be used for naming to match metadata.
```bash
bash get_SRR.sh $PRJNA# $OutputDirectory $OutputfileName
```
Example:
```bash
bash get_SRR.sh PRJNA884416 ~/scratch/user/data/srr_numbers SRR_sampleID.csv
```
Next, you can estimate download size required with vdb-dump
```bash
bash download_size.sh 2bRAD_SRR_sampleID.csv
```
Finally, download raw fastqs with the script ```download_fastq.sh```, the SRR_sampleID list, the output directory, and the logfile name to check successful download. Raw fastq files are named sampleID.fastq.
```bash
bash download_fastq.sh SRR/sampleIDlist outputDirectory
```
Example:
```bash
bash download_fastq.sh 2bRAD_SRR_sampleID.csv /scratch/user/data/raw_fastq
```
Alternatively, use the sbatch script to schedule the download via SLURM scheduler. 
<br>By default, it loads mamba, but if you are using conda, be sure to edit the script.
```bash
sbatch download_fastq.sbatch $mamba/conda_EnvironmentName $OutputDirectory $SRR_list
```
Example:
```bash
sbatch --mail-user=user@example.com download_fastq.sbatch download /sctach/user/data/raw_fastq 2bRAD_SRR_sampleID.csv
```
