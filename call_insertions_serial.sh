#!/bin/bash
if [ $# -eq 0 ]; then
    >&2 echo "No arguments provided"
    exit 1
fi

while getopts s:f:d:n: flag
do
    case "${flag}" in
        s) species=${OPTARG};;
        f) file=${OPTARG};;
        d) data_dir=${OPTARG};;
        n) threads=${OPTARG};;
    esac
done

data_path="$(pwd)/$data_dir"
file_path="$data_dir/$file"


mkdir -p output


if [ ! -f "$file_path" ]; then
    echo "Error: File '$file_path' does not exist."
    exit 3
fi

mkdir -p output/$species

while IFS= read -r line; do
    line=$(echo "$line" | tr -d '[:space:]')
    sra_example=$line

    mkdir -p output/$species/"${sra_example}"
    
    echo insurveyor.py --threads $threads data/"${sra_example}".bam output/$species/"${sra_example}" $data_path/$species.fa | cat -v
   
    insurveyor.py --threads $threads data/"${sra_example}".bam output/$species/"${sra_example}" $data_path/$species.fa

done < "$file_path"

