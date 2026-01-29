#!/bin/bash

SRR_LIST="$1"	# Input file of SRR numbers, one per line
OUTDIR="$2"	# Specify the name of the output directory to save files to
LOGFILE=$OUTDIR/log.txt	# The file to record SRR numbers successfully downloaded. This should match your original file.

# Create directory to save output files in
mkdir -p "$OUTDIR"

# Clear the log file to start
> "$LOGFILE"

# Loop through each accession
while IFS=',' read -r SRR SAMPLEID; do
	OUTFILE="$OUTDIR/${SAMPLEID}.fastq"
	# Downloads file into the output directory
	if fasterq-dump "$SRR" -o $OUTFILE; then
	# If successful, prints the SRR name to the log file
		echo "$SRR $SAMPLEID" | tee -a "$LOGFILE"
	else
		echo "Failed to download $SRR($SAMPLEID)" | tee -a "$LOGFILE"
	fi
done < "$SRR_LIST"

# Compare number of successful downloads to input SRR list
SUCCESS=$(wc -l < "$LOGFILE")
TOTAL=$(wc -l < "$SRR_LIST")

echo "Downloaded $SUCCESS out of $TOTAL SRRs."

if [ "$SUCCESS" -eq "$TOTAL" ]; then
    echo "All SRRs downloaded successfully!"
else
    echo "Some SRRs failed. Check the log."
fi
