#!/bin/bash
#SBATCH --job-name=download_fastq	# Job name
#SBATCH --time=04:00:00			# Maximum runtime (HH:MM:SS)
#SBATCH --ntasks=1			# Number of tasks (processes)
#SBATCH --cpus-per-task=2		# Number of CPU cores per task
#SBATCH --mem=8G			# Memory per node
#SBATCH --partition=public		# Partition/queue name
#SBATCH --qos=public			# QOS
#SBATCH --mail-type=END,FAIL		# Send email on job end or fail
#SBATCH --export=NONE			# Purge the job-submitting shell environment


# Path to your SRR list and output
ENV=$1
OUTDIR=$2
SRR_LIST=$3
LOGFILE=$OUTDIR/log.txt

# Load modules if necessary
module load mamba/latest
source activate download

# Run your script
bash download_fastq.sh $SRR_LIST $OUTDIR $LOGFILE
