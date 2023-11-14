#!/bin/bash
set -eu
set -o pipefail 

while getopts "s:d:f:" flag; do
    case "${flag}" in 
        s) species=${OPTARG};;
        d) species_dir=${OPTARG};;
        f) zipped_gtf_dataset=${OPTARG};;
        *) echo "Usage: $0 -s species_name -d species_dir -f zipped_gtf_dataset"
            exit 1;;
    esac
done 

# data_dir=$(grep 'data_directory:' configs/config_${species}.yaml | awk '{print $2}')
echo "n" | unzip $zipped_gtf_dataset
# extracting the genomic.gtf refseq file
gcf_folder=$(ls ncbi_dataset/data | grep GCF)
gtf_file=$(ls "ncbi_dataset/data/$gcf_folder/")

# moving the gtf file to data_dir/species_dir and renaming as species.gtf
mv "ncbi_dataset/data/$gcf_folder/$gtf_file" "$species_dir/${species}.gtf"
rm -rf ncbi_dataset