#!/bin/bash
set -eu
set -o pipefail 

if [ $# -eq 0 ]; then
    >&2 echo "No arguments provided"
    exit 1
fi

while getopts s:f:o:l:h: flag
do
    case "${flag}" in
        s) species=${OPTARG};;
        f) sample_file=${OPTARG};;
        o) output_dir=${OPTARG};;
        l) low=${OPTARG};;
        h) high=${OPTARG};;
        *) "Usage: $0 -s species -f sample_file -o output_dir -l low -h high"
    esac
done

mkdir -p "$output_dir"
mkdir -p "$output_dir/peak_$low-$high/"

rm -f "$output_dir/peak_$low-$high/output_seqs.vcf"

while read -r line; do
    sra_example=$line

    bcftools view -i  'SVLEN>='${low}' && SVLEN<='${high}'' "$output_dir/${sra_example}/out.pass.vcf" -o  \
    "$output_dir/peak_$low-$high/output_seqs.vcf"

    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%SVINSSEQ\t'${sra_example}'\n' "$output_dir/peak_$low-$high/output_seqs.vcf" >> \
            "$output_dir/peak_$low-$high/${species}_$low-${high}_pop_vcf.txt"

    rm -f "$output_dir/peak_$low-$high/output_seqs.vcf"
done < "$sample_file"



