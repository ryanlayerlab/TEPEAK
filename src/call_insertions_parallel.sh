#!/bin/bash
set -eu 
set -o pipefail

if [ $# -eq 0 ]; then
    >&2 echo "No arguments provided"
    exit 1
fi

while getopts s:n:d:l: flag; do
    case "${flag}" in
        l) sra_examples=${OPTARG};;
        s) species=${OPTARG};;
        d) data_dir=${OPTARG};;
        n) threads=${OPTARG};;
        *) echo "Usage: $0 -s species_name -d data_directory -n number_of_threads -l sra_file"
           exit 1;;
    esac
done

IFS=' ' read -ra items <<< "$sra_examples"
data_path="$(pwd)/$data_dir"
mkdir -p output/$species

for item in "${items[@]}"; do
    echo "Processing: $item"
    line=$(echo "$item" | tr -d '[:space:]')
    
    sra_example=$line
    mkdir -p output/$species/"${sra_example}"
    echo insurveyor.py --threads $threads data/"${sra_example}".bam output/$species/"${sra_example}" $data_path/$species.fa | cat -v
    insurveyor.py --threads $threads data/"${sra_example}".bam output/$species/"${sra_example}" $data_path/$species.fa
done
echo "new job"