#!/bin/bash
set -eu
set -o pipefail 

# Parse command line arguments
while getopts ":s:" opt; do
  case $opt in
    s) species="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
    :) echo "Option -$OPTARG requires an argument." >&2
    exit 1
    ;;
  esac
done


data_dir=$(grep 'data_directory:' configs/config_${species}.yaml | awk '{print $2}')
filename="$data_dir/$species/${species}_samples.txt"
output_file="${data_dir}/${species}/count_${species}.txt"
echo -e "Sample\tINS Count" > "$output_file"  # Initialize output file with headers

while read -r line; do
  vcf_file="output/$species/${line}/out.pass.vcf.gz"

  # Unzip the file -- not using gzip -k to maintain compatibility with gzip<=1.5
  decompressed_file="output/$species/${line}/out.pass.vcf"
  gzip -dc "$vcf_file" > "$decompressed_file"

  # Count occurrences of "SVLEN"
  count=$(grep -c "SVLEN" "$decompressed_file")

  # Write to output file
  echo -e "${line}\t${count}" >> "$output_file"

done < "$filename"

echo "Processing complete. Check $output_file for results."