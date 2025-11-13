#!/bin/bash

# Horse-specific duplication script
ORIGINAL_SAMPLE="SRR1564422"
NEW_SAMPLE="SRR1564422_test"

echo "Duplicating horse sample for parallelization testing..."
echo "Original sample: $ORIGINAL_SAMPLE"
echo "New sample: $NEW_SAMPLE"
echo

# Check if BAM files exist
if [ ! -f "data/horse/${ORIGINAL_SAMPLE}.bam" ]; then
    echo "Error: BAM file data/horse/${ORIGINAL_SAMPLE}.bam not found"
    echo "Available files in data/horse/:"
    ls -la data/horse/ | head -10
    exit 1
fi

echo "Duplicating BAM files..."
cp "data/horse/${ORIGINAL_SAMPLE}.bam" "data/horse/${NEW_SAMPLE}.bam"
cp "data/horse/${ORIGINAL_SAMPLE}.bam.bai" "data/horse/${NEW_SAMPLE}.bam.bai"
echo "  ✓ Copied ${ORIGINAL_SAMPLE}.bam -> ${NEW_SAMPLE}.bam"
echo "  ✓ Copied ${ORIGINAL_SAMPLE}.bam.bai -> ${NEW_SAMPLE}.bam.bai"

# Update sample list
echo "$NEW_SAMPLE" >> horse_samples.txt
echo "  ✓ Added $NEW_SAMPLE to horse_samples.txt"

echo
echo "Sample duplication complete!"
echo "You now have 2 samples for testing parallelization:"
echo "  1. $ORIGINAL_SAMPLE (original)"  
echo "  2. $NEW_SAMPLE (duplicate)"
echo
echo "Current samples in horse_samples.txt:"
cat horse_samples.txt
echo
echo "To test parallelization:"
echo "  1. Delete VCF outputs: rm -rf output/horse/*/out.pass.vcf.gz"
echo "  2. Run pipeline: snakemake --configfile config_horse_no_smoove.yaml --cores 8 --forcerun call_insertions_one_sample"
echo "  3. Watch for parallel execution in the logs"
