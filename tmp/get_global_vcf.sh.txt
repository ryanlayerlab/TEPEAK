#!/bin/bash
set -eu
set -o pipefail 

if [ $# -eq 0 ]; then
    >&2 echo "No arguments provided"
    exit 1
fi

while getopts s:o:f: flag; do
    case "${flag}" in
        s) species="$OPTARG";;
        o) output_dir="$OPTARG";;
        f) sample_file="$OPTARG";;
        *) "Usage: $0 -s species -o output_dir -f sample_file"
    esac
done

while read -r line; do
    decompressed_file="$output_dir/$line/out.pass.vcf"
	bcftools query -f "%CHROM\t%POS\t%INFO/END\t%SVLEN\t%SVINSSEQ\n" "$decompressed_file" >> "${species}_global_vcf.txt"
done < "$sample_file"
mv  "${species}_global_vcf.txt" "$output_dir/${species}_global_vcf.txt"

# python3 src/build_histogram.py -f output/"$species"/"$species"_global_vcf.txt
