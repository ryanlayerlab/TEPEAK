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

# Robust Picard detection - check multiple common locations
PICARD_JAR=""
POSSIBLE_PICARD_PATHS=(
    # Relative to TEPEAK directory (README instructions)
    "$TEPEAK_DIR/picard/build/libs/picard.jar"
    # Relative to current working directory
    "$(pwd)/picard/build/libs/picard.jar" 
    "./picard/build/libs/picard.jar"
    # Check if user hardcoded a path in environment
    "${PICARD_JAR_PATH:-}"
    # Common system locations
    "/usr/local/bin/picard.jar"
    "/opt/picard/picard.jar"
    "$HOME/picard/build/libs/picard.jar"
    # Check for picard command (conda or system install)
    "$(which picard 2>/dev/null || echo '')"
)

echo "Searching for Picard installation..."
for path in "${POSSIBLE_PICARD_PATHS[@]}"; do
    if [[ -n "$path" && -f "$path" ]]; then
        PICARD_JAR="$path"
        echo "Found Picard JAR at: $PICARD_JAR"
        break
    fi
done

# If no JAR found, try picard command
if [[ -z "$PICARD_JAR" ]]; then
    if command -v picard &> /dev/null; then
        PICARD_CMD="picard"
        echo "Using system picard command: $(which picard)"
    else
        echo "Error: Picard not found at any expected location:"
        printf '  %s\n' "${POSSIBLE_PICARD_PATHS[@]}"
        echo ""
        echo "Please install Picard following the README instructions:"
        echo "  git clone https://github.com/broadinstitute/picard.git"
        echo "  cd picard/"
        echo "  ./gradlew shadowJar"
        echo ""
        echo "Or set PICARD_JAR_PATH environment variable to the full path"
        exit 1
    fi
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
    
    # Use species directory for temporary files instead of /tmp
    sample_temp_dir="$data_path/${sra_example}_tmp"
    mkdir -p "$sample_temp_dir"
    cd "$sample_temp_dir"

    # Set cleanup trap
    trap "cd '$data_path' && rm -rf '$sample_temp_dir'" EXIT

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

    # Fix mate information using detected Picard
    echo "Fixing mate information..."
    mkdir -p picard_tmp
    
    # Set JAVA_OPTS to use our temp directory and increase memory
    export JAVA_OPTS="-Xmx8g -Djava.io.tmpdir=$(pwd)/picard_tmp"
    
    if [[ -n "$PICARD_JAR" ]]; then
        java -Djava.io.tmpdir="$(pwd)/picard_tmp" -Xmx8g -jar "$PICARD_JAR" FixMateInformation \
            I="${sra_example}".sorted.bam \
            ADD_MATE_CIGAR=true \
            O="${sra_example}".fixed.bam  \
            TMP_DIR="$(pwd)/picard_tmp"
    else
        picard -Djava.io.tmpdir="$(pwd)/picard_tmp" -Xmx8g FixMateInformation \
            I="${sra_example}".sorted.bam \
            ADD_MATE_CIGAR=true \
            O="${sra_example}".fixed.bam  \
            TMP_DIR="$(pwd)/picard_tmp"
    fi

    # Move final files to species directory
    mv "${sra_example}".fixed.bam ../"${sra_example}".bam
    mv "${sra_example}".fixed.bam.bai ../"${sra_example}".bam.bai
    
    echo "Completed $sra_example"
    
    # Cleanup
    cd "$data_path"
    rm -rf "$sample_temp_dir"
  
done < "$sra_file"
