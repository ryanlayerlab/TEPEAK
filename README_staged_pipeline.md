# TEPEAK Staged Pipeline

For large datasets or systems with limited `/tmp` space, TEPEAK can be run in two stages:

## Stage 1: Per-sample Processing
- Alignment (BWA)  
- Insertion calling (InsurVeyor)
- Can be run on sample subsets

## Stage 2: Merging and Analysis
- Global VCF creation
- TEPEAK clustering analysis
- Gene annotation
- Final output generation

## Usage

### Option 1: Split samples manually

1. **Split your sample file:**
```bash
# Split into chunks of 10 samples each
python scripts/split_samples.py horse_samples.txt chunks/horse --chunk-size 10

# Or split into 5 chunks total
python scripts/split_samples.py horse_samples.txt chunks/horse --num-chunks 5
```

2. **Run Stage 1 on each chunk:**
```bash
# For each chunk file created:
cp chunks/horse_chunk01.txt horse_samples.txt
snakemake -s Snakefile_stage1.py --configfile config_horse.yaml --cores 32

cp chunks/horse_chunk02.txt horse_samples.txt  
snakemake -s Snakefile_stage1.py --configfile config_horse.yaml --cores 32

# ... repeat for all chunks
```

3. **Run Stage 2 (merging):**
```bash
# This automatically detects all completed VCF files
snakemake -s Snakefile_stage2.py --configfile config_horse.yaml --cores 64
```

### Option 2: Manual subset selection

1. **Create subset sample files:**
```bash
# Create smaller sample files manually
head -5 horse_samples.txt > horse_samples_batch1.txt
tail -n +6 horse_samples.txt | head -5 > horse_samples_batch2.txt
```

2. **Run Stage 1 batches:**
```bash
# Batch 1
cp horse_samples_batch1.txt horse_samples.txt
snakemake -s Snakefile_stage1.py --configfile config_horse.yaml --cores 64

# Batch 2  
cp horse_samples_batch2.txt horse_samples.txt
snakemake -s Snakefile_stage1.py --configfile config_horse.yaml --cores 64
```

3. **Run Stage 2:**
```bash
snakemake -s Snakefile_stage2.py --configfile config_horse.yaml --cores 64
```

## Benefits

- **Disk space management**: Smaller batches use less `/tmp` space
- **Resource flexibility**: Can run batches on different systems
- **Fault tolerance**: If one batch fails, others continue
- **Scheduling friendly**: Easier to fit into HPC queue systems

## Stage 1 Outputs

Each sample produces:
- `data/{species}/{sample}.bam` - Aligned reads
- `data/{species}/{sample}.bam.bai` - BAM index  
- `output/{species}/{sample}/out.pass.vcf.gz` - Insertion calls

## Stage 2 Requirements

- All desired sample VCF files must exist
- Original sample file should contain all samples for proper merging
- Reference files from Stage 1 must be present

## Troubleshooting

**"No samples detected in Stage 2":**
- Ensure VCF files exist in `output/{species}/*/out.pass.vcf.gz`
- Check that Stage 1 completed successfully for all batches

**"Missing VCF files":**
- Some Stage 1 jobs may have failed
- Rerun failed batches or check logs

**"Reference files missing":**
- Stage 2 expects reference files from Stage 1
- If needed, run `snakemake process_reference` first
