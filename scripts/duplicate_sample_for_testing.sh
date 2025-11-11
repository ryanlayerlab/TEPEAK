#!/bin/bash
# Simple script to duplicate eFus sample for testing

# Get the original sample name from the sample file
ORIGINAL_SAMPLE=$(head -1 eFus_samples.txt)
NEW_SAMPLE="${ORIGINAL_SAMPLE}_test"

echo "Duplicating sample for parallelization testing..."
echo "Original sample: $ORIGINAL_SAMPLE"
echo "New sample: $NEW_SAMPLE"
echo

# Duplicate files in the species directory
echo "Copying BAM files..."
if [ -f "data/eFus/${ORIGINAL_SAMPLE}.bam" ]; then
    cp "data/eFus/${ORIGINAL_SAMPLE}.bam" "data/eFus/${NEW_SAMPLE}.bam"
    echo "  Copied ${ORIGINAL_SAMPLE}.bam -> ${NEW_SAMPLE}.bam"
fi

if [ -f "data/eFus/${ORIGINAL_SAMPLE}.bam.bai" ]; then
    cp "data/eFus/${ORIGINAL_SAMPLE}.bam.bai" "data/eFus/${NEW_SAMPLE}.bam.bai"
    echo "  Copied ${ORIGINAL_SAMPLE}.bam.bai -> ${NEW_SAMPLE}.bam.bai"
fi

# Add new sample to sample list
echo "$NEW_SAMPLE" >> eFus_samples.txt
echo "  Added $NEW_SAMPLE to eFus_samples.txt"

echo
echo "Sample duplication complete!"
echo "You now have 2 samples for testing parallelization:"
echo "  1. $ORIGINAL_SAMPLE (original)"  
echo "  2. $NEW_SAMPLE (duplicate)"
echo
echo "To test parallelization:"
echo "  1. Delete VCF outputs: rm -rf output/eFus/*/out.pass.vcf.gz"
echo "  2. Run pipeline: snakemake --configfile config_eFus.yaml --cores 64"
echo "  3. Watch for parallel execution in the logs"
