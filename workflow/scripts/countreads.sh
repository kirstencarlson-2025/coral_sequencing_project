#!/bin/bash

# This script count number of reads for each fastq file and save the file name and number of reads in a text file

# Set variables
extension=".fastq"
output="ReadCounts.txt"

> "$output"

for file in *"$extension"; do
	lines=$(wc -l < "$file")
	reads=$((lines/4))
	echo "${file}	${reads}" >> "$output"
done

sum=$(awk '{sum += $2} END {print sum}' "$output")

echo "Total reads: ${sum}"
