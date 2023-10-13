#!/bin/bash

if [ $# -eq 0 ]; then
    >&2 echo "No arguments provided"
    exit 1
fi

while getopts s:f: flag
do
    case "${flag}" in
        s) species=${OPTARG};;
        f) file=${OPTARG};;
    esac
done


while read -r line; do
	vcf_file="output/$species/${line}/out.pass.vcf.gz"
	echo $vcf_file
	gzip -d "$vcf_file"
	decompressed_file="output/$species/${line}/out.pass.vcf"
	bcftools query -f "%CHROM\t%POS\t%INFO/END\t%SVLEN\t%SVINSSEQ\n" "$decompressed_file" >> "$species"_global_vcf.txt

done < "$file"

python3 buildHistogram.py -f "$species"_global_vcf.txt

