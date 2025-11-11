#!/bin/bash
#SBATCH --job-name=download_fastq	# Job name
#SBATCH --output=/home/kcarls36/scratch/get_data/download_fastq_%j.log	# Standard output and error log (%j = job ID)
#SBATCH --time=04:00:00			# Maximum runtime (HH:MM:SS)
#SBATCH --ntasks=1			# Number of tasks (processes)
#SBATCH --cpus-per-task=2		# Number of CPU cores per task
#SBATCH --mem=8G			# Memory per node
#SBATCH --partition=public		# Partition/queue name
#SBATCH --qos=public			# QOS
#SBATCH --mail-type=END,FAIL		# Send email on job end or fail
#SBATCH --mail-user=kcarls36@asu.edu
#SBATCH --export=NONE			# Purge the job-submitting shell environment

# Load modules if necessary
module load mamba/latest
source activate get_data

# Path to your SRR list and output
SRR_LIST=/home/kcarls36/get_data/2bRAD_SRR.csv
LOGFILE=/scratch/kcarls36/get_data/downloaded_files.txt

# Run your script
bash /scratch/kcarls36/get_data/fastqer.sh "$SRR_LIST" fastq_raw "$LOGFILE"

