#!/bin/bash
set -eu
set -o pipefail 

# Parse command line arguments
while getopts "s:o:f:" opt; do
  case $opt in
    s) species="$OPTARG";;
    o) output_dir="$OPTARG";;
    f) sample_file="$OPTARG";;
    *) echo "Usage: $0 -s species -o output_dir -t threads -f sample_file"
    exit 1;;
  esac
done

output_file="$output_dir/${species}_count.txt"
echo -e "Sample\tINS Count" > "$output_file"  # Initialize output file with headers

while read -r line; do
  vcf_file="$output_dir/${line}/out.pass.vcf.gz"

  # Unzip the file -- not using gzip -k to maintain compatibility with gzip<=1.5
  decompressed_file="$output_dir/${line}/out.pass.vcf"
  gzip -dc "$vcf_file" > "$decompressed_file"

  # Count occurrences of "SVLEN"
  count=$(grep -c "SVLEN" "$decompressed_file")

  # Write to output file
  echo -e "${line}\t${count}" >> "$output_file"
done < "$sample_file"

echo "Processing complete. Check $output_file for results."