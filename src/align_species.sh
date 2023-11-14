#!/bin/bash
set -eu
set -o pipefail

if [ $# -eq 0 ]; then
    >&2 echo "No arguments provided"
    exit 1
fi

while getopts s:d:t:f: flag; do
    case "${flag}" in
        s) species=${OPTARG};;
        d) data_dir=${OPTARG};;
        t) threads=${OPTARG};;
        f) sra_file=${OPTARG};;
        *) echo "Usage: $0 -s species -d species_dir -t threads -f sra_file"
            exit 1;;
    esac
done

data_path="$(pwd)/$data_dir"
picard_path="$(pwd)"

# sra_file=$data_dir/$species/${species}_samples.txt
while IFS= read -r line; do

    line=$(echo "$line" | tr -d '[:space:]')
    sra_example=$line

    mkdir -p "$data_path/$sra_example/bwa_errors"
    cd "$data_path/$sra_example"

    prefetch "${sra_example}" --max-size 200G
    fastq-dump --split-files "${sra_example}"

    bwa mem -M -t $threads -R "@RG\tID:1\tSM:""${sra_example}" \
         "$data_path/$species" \
         "${sra_example}"_1.fastq "${sra_example}"_2.fastq  \
         2> bwa_errors/bwa_"${sra_example}".err \
        > "${sra_example}".bam

    samtools sort "${sra_example}".bam -o "${sra_example}".sorted.bam #-@$threads
    samtools index "${sra_example}".sorted.bam #-@$threads

    java -jar "$picard_path"/picard/build/libs/picard.jar FixMateInformation \
   	I="${sra_example}".sorted.bam \
   	ADD_MATE_CIGAR=true \
    	O="${sra_example}".fixed.bam  \
    	TMP_DIR="$(pwd)/tmp"

    samtools index "${sra_example}".fixed.bam

    mv "${sra_example}".fixed.bam ../"${sra_example}".bam
    mv "${sra_example}".fixed.bam.bai ../"${sra_example}".bam.bai
    
    cd ../..
    rm -r "$data_path/$sra_example"
  
done < "$sra_file"


