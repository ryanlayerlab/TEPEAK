#!/bin/bash
set -eu
set -o pipefail 

if [ $# -eq 0 ]; then
    >&2 echo "No arguments provided"
    exit 1
fi

while getopts s: flag
do
    case "${flag}" in
        s) species=${OPTARG};;
    esac
done

data_dir=$(grep 'data_directory:' configs/config_${species}.yaml | awk '{print $2}')
echo $data_dir
data_path="$(pwd)/$data_dir/${species}"
# threads=$(grep 'threads:' configs/config_${species}.yaml | awk '{print $2}')
picard_path="$(pwd)"

sra_file=$data_dir/$species/${species}_samples.txt


while read -r line; do
	vcf_file="output/$species/${line}/out.pass.vcf.gz"
	echo $vcf_file
	gzip -d "$vcf_file"
	decompressed_file="output/$species/${line}/out.pass.vcf"
	bcftools query -f "%CHROM\t%POS\t%INFO/END\t%SVLEN\t%SVINSSEQ\n" "$decompressed_file" >> "$species"_global_vcf.txt

done < "$sra_file"
mv  "$species"_global_vcf.txt output/"$species"/"$species"_global_vcf.txt

python3 src/build_histogram.py -f output/"$species"/"$species"_global_vcf.txt

python3 src/dfam_annotate.py -s "$species"
