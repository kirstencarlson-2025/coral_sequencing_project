#!/bin/bash
# WHAT: uses entrez-direct package to search for SRA runs using PRJNA accession number, then saves in csv file. 
# Searches for SRR numbers and saves in seperate 2bRAD and ITS2 files
# Must have entrez-direct package installed. For Kirsten, use: source activate get_data

# Set variables
PRJNA=$1     #PRJNA accession number
OUTPUT=$2    #output file name .csv

# Fetch run info from PRJNA accession number, save to output
esearch -db sra -query $PRJNA | efetch -format runinfo > $OUTPUT


# Save only SRR numbers to a csv file, to use for later download
#grep 2bRAD $OUTPUT | cut -f1 -d "," | grep SRR > 2bRAD_$OUTPUT
#grep ITS2 $OUTPUT | cut -f1 -d "," | grep SRR > ITS2_$OUTPUT

# Save SRR number and SampleID for later fastq file rename
cut -d "," -f1,30 $OUTPUT | grep 2bRAD > 2bRAD_SampleID_$OUTPUT
