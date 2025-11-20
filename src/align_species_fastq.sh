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
else
    PICARD_CMD="java -jar \"$PICARD_JAR\""
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

# Create temporary directory for this sample - USE SPECIES DIR INSTEAD OF /tmp
temp_dir="${species_dir}/tmp_${sample_name}_$$"
mkdir -p "$temp_dir"

# Ensure cleanup happens even if script fails
trap "rm -rf '$temp_dir'" EXIT

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

# Fix mate information using detected Picard
echo "Fixing mate information..."
mkdir -p "${temp_dir}/picard_tmp"

# Set JAVA_OPTS to use our temp directory and increase memory
export JAVA_OPTS="-Xmx8g -Djava.io.tmpdir=${temp_dir}/picard_tmp"

if [[ -n "$PICARD_JAR" ]]; then
    java -Djava.io.tmpdir="${temp_dir}/picard_tmp" -Xmx8g -jar "$PICARD_JAR" FixMateInformation \
        I="${temp_dir}/${sample_name}.sorted.bam" \
        ADD_MATE_CIGAR=true \
        O="${temp_dir}/${sample_name}.fixed.bam" \
        TMP_DIR="${temp_dir}/picard_tmp"
else
    picard -Djava.io.tmpdir="${temp_dir}/picard_tmp" -Xmx8g FixMateInformation \
        I="${temp_dir}/${sample_name}.sorted.bam" \
        ADD_MATE_CIGAR=true \
        O="${temp_dir}/${sample_name}.fixed.bam" \
        TMP_DIR="${temp_dir}/picard_tmp"
fi

# Index fixed BAM
samtools index "${temp_dir}/${sample_name}.fixed.bam"

# Move final files to output location
mv "${temp_dir}/${sample_name}.fixed.bam" "${species_dir}/${sample_name}.bam"
mv "${temp_dir}/${sample_name}.fixed.bam.bai" "${species_dir}/${sample_name}.bam.bai"

echo "Successfully processed sample: $sample_name"

# Cleanup is handled by trap
echo "Temporary files cleaned up from: $temp_dir"

echo "BWA alignment log saved to: ${temp_dir}/bwa_${sample_name}.err"
