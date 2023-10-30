#!/bin/bash
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



mkdir -p output
mkdir -p output/$species


BAM_DIR_VALUE=$(grep '^bam_dir:' configs/config_${species}.yaml | awk '{print $2}')

if [ -n "$BAM_DIR_VALUE" ]; then
    bam_dir=$BAM_DIR_VALUE
    data_path="$(pwd)/$bam_dir"

    # You can now use the $bam_dir variable in your script
else
     data_dir=$(grep 'data_directory:' configs/config_${species}.yaml | awk '{print $2}')
     data_path="$(pwd)/$data_dir/${species}"
fi

data_dir=$(grep 'data_directory:' configs/config_${species}.yaml | awk '{print $2}')

#data_path="$(pwd)/$data_dir/${species}"

threads=$(grep 'threads:' configs/config_${species}.yaml | awk '{print $2}')

sra_file=$data_dir/$species/${species}_samples.txt


#mkdir -p output/$species
while IFS= read -r line; do
    line=$(echo "$line" | tr -d '[:space:]')
    sra_example=$line

    mkdir -p output/$species/"${sra_example}"
    
    echo insurveyor.py --threads $threads "$data_path"/"${sra_example}".bam output/$species/"${sra_example}" "$data_path"/$species.fa | cat -v
   
    insurveyor.py --threads $threads "$data_path"/"${sra_example}".bam output/$species/"${sra_example}" "$data_path"/$species.fa

done < "$sra_file"

