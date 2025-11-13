#!/bin/bash
set -e

# FASTQ-based alignment script for per-sample parallelization
# This handles local FASTQ files for a single sample

if [ $# -eq 0 ]; then
    >&2 echo "No arguments provided"
    exit 1
fi

while getopts s:d:t:f:l: flag; do
    case "${flag}" in
        s) species=${OPTARG};;
        d) species_dir=${OPTARG};;
        t) threads=${OPTARG};;
        f) fastq_dir=${OPTARG};;
        l) sample_file=${OPTARG};;
        *) echo "Usage: $0 -s species -d species_dir -t threads -f fastq_dir -l sample_file"
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

# Read the single sample from the temp file
sample_name=$(head -1 "$sample_file")
sample_name=$(echo "$sample_name" | tr -d '[:space:]')

echo "Processing FASTQ sample: $sample_name"
echo "Species directory: $species_dir"
echo "FASTQ directory: $fastq_dir"
echo "Threads: $threads"

# Find FASTQ files (try multiple naming conventions)
R1_FILE=""
R2_FILE=""

for pattern in "${sample_name}_R1" "${sample_name}_1" "${sample_name}.R1" "${sample_name}-R1"; do
    if [ -f "${fastq_dir}/${pattern}.fastq.gz" ]; then
        R1_FILE="${fastq_dir}/${pattern}.fastq.gz"
        break
    elif [ -f "${fastq_dir}/${pattern}.fastq" ]; then
        R1_FILE="${fastq_dir}/${pattern}.fastq"
        break
    fi
done

for pattern in "${sample_name}_R2" "${sample_name}_2" "${sample_name}.R2" "${sample_name}-R2"; do
    if [ -f "${fastq_dir}/${pattern}.fastq.gz" ]; then
        R2_FILE="${fastq_dir}/${pattern}.fastq.gz"
        break
    elif [ -f "${fastq_dir}/${pattern}.fastq" ]; then
        R2_FILE="${fastq_dir}/${pattern}.fastq"
        break
    fi
done

if [ -z "$R1_FILE" ] || [ -z "$R2_FILE" ]; then
    echo "Error: Could not find FASTQ pair for sample $sample_name"
    echo "Looked for: ${sample_name}_R[12], ${sample_name}_[12], ${sample_name}.R[12], ${sample_name}-R[12]"
    echo "In directory: $fastq_dir"
    echo "Available files:"
    ls -la "$fastq_dir" | head -10
    exit 1
fi

echo "Found FASTQ files:"
echo "  R1: $R1_FILE"
echo "  R2: $R2_FILE"

# Create temporary directory for this sample
temp_dir="/tmp/tepeak_${sample_name}_$$"
mkdir -p "$temp_dir"

# Align with BWA
echo "Running BWA alignment..."
bwa mem -M -t "$threads" -R "@RG\tID:${sample_name}\tSM:${sample_name}\tLB:${sample_name}\tPL:illumina" \
    "${species_dir}/${species}.fa" \
    "$R1_FILE" "$R2_FILE" \
    2> "${temp_dir}/bwa_${sample_name}.err" | \
    samtools sort -@ "$threads" -T "${temp_dir}/sort_${sample_name}" -o "${temp_dir}/${sample_name}.sorted.bam" -

# Index sorted BAM
echo "Indexing BAM file..."
samtools index "${temp_dir}/${sample_name}.sorted.bam"

# Fix mate information using local Picard build
echo "Fixing mate information..."
mkdir -p "${temp_dir}/picard_tmp"
java -jar "$PICARD_JAR" FixMateInformation \
    I="${temp_dir}/${sample_name}.sorted.bam" \
    ADD_MATE_CIGAR=true \
    O="${temp_dir}/${sample_name}.fixed.bam" \
    TMP_DIR="${temp_dir}/picard_tmp"

# Index fixed BAM
samtools index "${temp_dir}/${sample_name}.fixed.bam"

# Move final files to output location
mv "${temp_dir}/${sample_name}.fixed.bam" "${species_dir}/${sample_name}.bam"
mv "${temp_dir}/${sample_name}.fixed.bam.bai" "${species_dir}/${sample_name}.bam.bai"

echo "Successfully processed sample: $sample_name"

# Cleanup temporary files
rm -rf "$temp_dir"

echo "BWA alignment log saved to: ${temp_dir}/bwa_${sample_name}.err"
