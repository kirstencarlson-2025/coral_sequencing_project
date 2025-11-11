# What are these scripts?

### **get_SRR.sh** - this uses a PRJNA accession number to create csv files listing all project SRR numbers in 2bRAD and ITS2 files
### **download_size.sh** - this uses the SRR csv files to calculate download size for fastq file download
### **download_fastq.sh** - this uses SRR csv files to download raw fastq files
### **download_fastq.sbatch** - sbatch script for SLURM scheduler
### **sample_rename** - this renames raw fastq files from SRR____.fastq to sampleID.fastq
### **Snakefile** - this runs the pipeline to call variants from raw fastq files
