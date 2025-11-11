#!/bin/bash

# This script creates a csv file with SRR number and Sample ID
# Then uses the csv file to rename fastq files from SRR.fastq to SampleID.fastq

# Example: bash sample_rename.sh

# Set variables
FASTQDIR=/scratch/kcarls36/projects/data/fastq_raw
SAMPLELIST=~/projects/stephanocoenia_FKNMS/data/SRR_numbers/sample_rename.csv
#######################
# Create a csv file with SRR numbers and sampleID
#######################

esearch -db sra -query PRJNA884416 | efetch -format runinfo > PRJNA884416_runinfo.csv
grep 2bRAD PRJNA884416_runinfo.csv | cut -d "," -f1,30 > $SAMPLELIST

# Field triplicates didn't retain IDs, so I had to go through NCBI to rename them in the sample_rename.csv file manually.
# They are:
# SFK066.1 = SRR21715609
# SFK066.2 = SRR21715608
# SFK066.3 = SRR21715607
# SFK162.1 = SRR21715592
# SFK162.2 = SRR21715591
# SFK162.3 = SRR21715590
# SFK205.1 = SRR21715481
# SKF205.2 = SRR21715480
# SFK205.3 = SRR21715478


########################
# Rename fastq files
########################

# Loop through CSV
while IFS=',' read -r srr sampleid;
do
    old_file="$FASTQDIR/$srr.fastq"
    new_file="$FASTQDIR/$sampleid.fastq"

   if [ -f "$old_file" ]; then
       mv "$old_file" "$new_file"
    else
        echo "WARNING: $old_file does not exist, skipping"
    fi
done < $SAMPLELIST
