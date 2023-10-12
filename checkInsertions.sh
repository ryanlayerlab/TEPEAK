#!/bin/bash

# Parse command line arguments
while getopts ":s:f:" opt; do
  case $opt in
    s) species="$OPTARG"
    ;;
    f) filename="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
    :) echo "Option -$OPTARG requires an argument." >&2
    exit 1
    ;;
  esac
done

# Check if species and filename have been provided
if [ -z "$species" ] || [ -z "$filename" ]; then
    echo "Both -s (species) and -f (filename) must be provided."
    exit 1
fi

# Process each line of the file
output_file="count_${species}.txt"
echo -e "Sample\tINS Count" > "$output_file"  # Initialize output file with headers

while read -r line; do
    vcf_file="output/$species/${line}/out.pass.vcf.gz"
    echo $vcf_file    
    # Unzip the file
    gzip -d "$vcf_file"
    decompressed_file="output/$species/${line}/out.pass.vcf"

    # Count occurrences of "SVLEN"
    count=$(grep -c "SVLEN" "$decompressed_file")

    # Write to output file
    echo -e "${line}\t${count}" >> "$output_file"

done < "$filename"

echo "Processing complete. Check $output_file for results."

