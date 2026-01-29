#!/bin/bash
# WHAT: uses entrez-direct package to search for SRA runs using PRJNA accession number, then saves in csv file. 
# Searches for SRR numbers and saves in seperate 2bRAD and ITS2 files
# Must have entrez-direct package installed.

# Set variables
PRJNA=$1     #PRJNA accession number
OUTPUTDIR=$2 #ouput directory
OUTPUT=$3    #output file name .csv

# Fetch run info from PRJNA accession number, save to output
esearch -db sra -query $PRJNA | efetch -format runinfo > $OUTPUTDIR/$OUTPUT


# Save only SRR numbers to a csv file, to use for later download
grep 2bRAD $OUTPUTDIR/$OUTPUT | cut -f1,30 -d "," | grep SRR > 2bRAD_$OUTPUT
grep ITS2 $OUTPUTDIR/$OUTPUT | cut -f1,30 -d "," | grep SRR > ITS2_$OUTPUT
