#!/bin/bash
# This script calculates size of each SRR file and returns total size of all files
# This will tell you how much space you need for downloading your data
# Packages needed:
# - ncbi-vdb

SRR_LIST="$1"  # your CSV or TXT file with SRRs
TOTAL_SIZE=0   # initialize total size in bytes

while read -r SRR; do
    # Get size of SRR in bytes using vdb-dump
    SIZE=$(vdb-dump "$SRR" --info | grep -i "size" | awk '{print $3}' | tr -cd '0-9')

    if [ -n "$SIZE" ]; then
        TOTAL_SIZE=$((TOTAL_SIZE + SIZE))
        echo "$SRR: $SIZE bytes"
    else
        echo "Warning: Could not get size for $SRR"
    fi
done < "$SRR_LIST"

# Convert total bytes to GB
TOTAL_GB=$(echo "scale=2; $TOTAL_SIZE/1024/1024/1024" | bc)
echo "Estimated total download size: $TOTAL_GB GB"
