#!/bin/bash
set -eu
set -o pipefail

# SRA-based alignment script for per-sample parallelization
# This handles downloading from SRA and aligning a single sample

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

# Get script directory for relative paths to Picard
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEPEAK_DIR="$(dirname "$SCRIPT_DIR")"

# Use local Picard build as per README instructions
PICARD_JAR="$TEPEAK_DIR/picard/build/libs/picard.jar"

# Check if Picard jar exists
if [ ! -f "$PICARD_JAR" ]; then
    echo "Error: Picard not found at $PICARD_JAR"
    echo "Please install Picard following the README instructions:"
    echo "  git clone https://github.com/broadinstitute/picard.git"
    echo "  cd picard/"
    echo "  ./gradlew shadowJar"
    exit 1
fi

echo "Using Picard: $PICARD_JAR"

data_path="$(pwd)/$data_dir"

# Check if picard is available
if ! command -v picard &> /dev/null; then
    echo "Error: picard not found. Please install via conda: conda install -c bioconda picard"
    exit 1
fi

while IFS= read -r line; do
    line=$(echo "$line" | tr -d '[:space:]')
    sra_example=$line

    echo "Processing SRA sample: $sra_example"
    
    mkdir -p "$data_path/$sra_example/bwa_errors"
    cd "$data_path/$sra_example"

    # Download and extract FASTQ
    echo "Downloading $sra_example..."
    prefetch "${sra_example}" --max-size 200G
    fastq-dump --split-files "${sra_example}"

    # Align with BWA
    echo "Aligning $sra_example..."
    bwa mem -M -t $threads -R "@RG\tID:1\tSM:${sra_example}" \
         "$data_path/${species}.fa" \
         "${sra_example}"_1.fastq "${sra_example}"_2.fastq  \
         2> bwa_errors/bwa_"${sra_example}".err \
        > "${sra_example}".bam

    # Sort and index
    echo "Sorting and indexing..."
    samtools sort "${sra_example}".bam -o "${sra_example}".sorted.bam
    samtools index "${sra_example}".sorted.bam

    # Fix mate information using local Picard build
    echo "Fixing mate information..."
    mkdir -p tmp
    java -jar "$PICARD_JAR" FixMateInformation \
        I="${sra_example}".sorted.bam \
        ADD_MATE_CIGAR=true \
        O="${sra_example}".fixed.bam  \
        TMP_DIR="$(pwd)/tmp"

    samtools index "${sra_example}".fixed.bam

    # Move final files to species directory
    mv "${sra_example}".fixed.bam ../"${sra_example}".bam
    mv "${sra_example}".fixed.bam.bai ../"${sra_example}".bam.bai
    
    echo "Completed $sra_example"
    
    # Cleanup
    cd ../..
    rm -r "$data_path/$sra_example"
  
done < "$sra_file"
