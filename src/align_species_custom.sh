#!/bin/bash
set -x


if [ $# -eq 0 ]; then
    >&2 echo "No arguments provided"
    exit 1
fi

while getopts s:d:t:f:q: flag; do
    case "${flag}" in
        s) species=${OPTARG};;
        d) data_dir=${OPTARG};;
        t) threads=${OPTARG};;
        f) sample_file=${OPTARG};;
        q) fastq_dir=${OPTARG};;
        *) echo "Usage: $0 -s species -d species_dir -t threads -f sample_file -q fastq_dir"
            exit 1;;
    esac
done

data_path="$(pwd)/$data_dir"
picard_path="$(pwd)"
fastq_path="$(pwd)/$fastq_dir"

echo "DEBUG: sample_file = $sample_file"
echo "DEBUG: sample_file exists = $(test -f "$sample_file" && echo yes || echo no)"
echo "DEBUG: sample_file contents:"
cat "$sample_file" || echo "Failed to read sample file"

while IFS= read -r line; do

    line=$(echo "$line" | tr -d '[:space:]')
    sample_name=$line
    
    echo "DEBUG: Processing sample: $sample_name"
    echo "DEBUG: Creating directory: $data_path/$sample_name/bwa_errors"
    mkdir -p "$data_path/$sample_name/bwa_errors"
    echo "DEBUG: Changing to directory: $data_path/$sample_name"
    cd "$data_path/$sample_name"

    # Copy and decompress FASTQ files instead of downloading from SRA
    echo "Processing custom FASTQ files for sample: $sample_name"
    
    if [ -f "$fastq_path/${sample_name}_R1.fastq.gz" ]; then
        gunzip -c "$fastq_path/${sample_name}_R1.fastq.gz" > "${sample_name}_1.fastq"
    else
        echo "Error: ${sample_name}_R1.fastq.gz not found in $fastq_path"
        exit 1
    fi
    
    if [ -f "$fastq_path/${sample_name}_R2.fastq.gz" ]; then
        gunzip -c "$fastq_path/${sample_name}_R2.fastq.gz" > "${sample_name}_2.fastq"
    else
        echo "Error: ${sample_name}_R2.fastq.gz not found in $fastq_path"
        exit 1
    fi

 
    bwa mem -M -t $threads -R "@RG\tID:1\tSM:""${sample_name}" \
             "$data_path/${species}.fa" \
         "${sample_name}"_1.fastq "${sample_name}"_2.fastq  \
         2> bwa_errors/bwa_"${sample_name}".err \
        > "${sample_name}".bam

    samtools sort "${sample_name}".bam -o "${sample_name}".sorted.bam #-@$threads
    samtools index "${sample_name}".sorted.bam #-@$threads

    java -jar "$picard_path"/picard/build/libs/picard.jar FixMateInformation \
   	I="${sample_name}".sorted.bam \
   	ADD_MATE_CIGAR=true \
    	O="${sample_name}".fixed.bam  \
    	TMP_DIR="$(pwd)/tmp"

    samtools index "${sample_name}".fixed.bam

    mv "${sample_name}".fixed.bam ../"${sample_name}".bam
    mv "${sample_name}".fixed.bam.bai ../"${sample_name}".bam.bai
    
    cd ../..
    rm -r "$data_path/$sample_name"
  
done < "$sample_file"