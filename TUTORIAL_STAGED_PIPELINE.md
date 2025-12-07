# TEPEAK Staged Pipeline Tutorial
## Running Multiple Samples with Chunked Processing

This tutorial demonstrates how to run TEPEAK with multiple samples using the new staged pipeline approach that prevents memory and disk space issues.

---

## 📋 **Tutorial Overview**

We'll process **4 horse samples** split into **2 jobs** to demonstrate:
- Sample splitting and chunking
- Concurrent Stage 1 execution 
- Automatic Stage 2 sample detection and merging
- Resource efficiency and collision prevention

---

## 🔧 **Prerequisites**

```bash
# 1. Activate the conda environment
conda activate insurveyor-env

# 2. Ensure you're in the TEPEAK directory
cd /path/to/TEPEAK

# 3. Verify required files exist
ls config_horse_no_smoove.yaml  # Configuration file
ls horse_samples.txt            # Sample list (we'll expand this)
ls data/horse_reference.zip     # Reference genome
ls data/horse_gtf.zip          # GTF annotations
```

---

## 📝 **Step 1: Create Sample List with 4 Samples**

Let's create a sample file with 4 horse samples for this tutorial:

```bash
# Create a 4-sample list for the tutorial
cat > tutorial_horse_samples.txt << 'EOF'
SRR1564422
SRR1564423
SRR1564424
SRR1564425
EOF

echo "Created sample list:"
cat tutorial_horse_samples.txt
```

**Expected Output:**
```
SRR1564422
SRR1564423
SRR1564424
SRR1564425
```

---

## 🔀 **Step 2: Split Samples into Chunks**

Use the built-in splitting script to create 2 chunks of 2 samples each:

```bash
# Split 4 samples into chunks of 2
python scripts/split_samples.py tutorial_horse_samples.txt chunks/horse --chunk-size 2

echo "Created chunk files:"
ls -la chunks/horse_chunk*.txt
```

**Expected Output:**
```
Created 2 chunks with ~2 samples each
  Chunk 1: 2 samples → chunks/horse_chunk01.txt  
  Chunk 2: 2 samples → chunks/horse_chunk02.txt

-rw-r--r-- 1 user staff 22 Dec  7 horse_chunk01.txt
-rw-r--r-- 1 user staff 22 Dec  7 horse_chunk02.txt
```

**Verify chunk contents:**
```bash
echo "Chunk 1 samples:"
cat chunks/horse_chunk01.txt

echo "Chunk 2 samples:"
cat chunks/horse_chunk02.txt
```

---

## 🔧 **Step 3: Update Configuration**

Create a tutorial-specific config file:

```bash
# Copy and modify the horse config for tutorial
cp config_horse_no_smoove.yaml tutorial_config.yaml

# Update the sample list reference
sed -i '' 's/horse_samples.txt/tutorial_horse_samples.txt/g' tutorial_config.yaml

echo "Configuration updated:"
grep "sample_list" tutorial_config.yaml
```

---

## 🚀 **Step 4: Run Stage 1 - Concurrent Chunk Processing**

Now we'll demonstrate running both chunks simultaneously to show the collision prevention in action.

### **Terminal 1 - Process Chunk 1:**
```bash
# Copy chunk 1 to the expected sample file location
cp chunks/horse_chunk01.txt tutorial_horse_samples.txt

echo "Processing Chunk 1 (Samples: $(cat tutorial_horse_samples.txt | tr '\n' ' '))"

# Run Stage 1 for chunk 1
snakemake -s Snakefile_stage1.py \
  --configfile tutorial_config.yaml \
  --cores 4 \
  --printshellcmds
```

### **Terminal 2 - Process Chunk 2 (run simultaneously):**
```bash
# Wait a few seconds, then in a new terminal:
sleep 5

# Copy chunk 2 to a different sample file name to avoid conflicts
cp chunks/horse_chunk02.txt tutorial_horse_samples_chunk2.txt

# Temporarily modify config for chunk 2
sed 's/tutorial_horse_samples.txt/tutorial_horse_samples_chunk2.txt/g' tutorial_config.yaml > tutorial_config_chunk2.yaml

echo "Processing Chunk 2 (Samples: $(cat tutorial_horse_samples_chunk2.txt | tr '\n' ' '))"

# Run Stage 1 for chunk 2  
snakemake -s Snakefile_stage1.py \
  --configfile tutorial_config_chunk2.yaml \
  --cores 4 \
  --printshellcmds
```

---

## 👀 **What You'll Observe During Concurrent Execution**

### **Reference Building Coordination:**
```
# Terminal 1 output:
[timestamp] rule process_reference:
    ...
Created lock file: data/horse/.reference_building.lock
[bwa_index] Pack FASTA... 10.35 sec
[bwa_index] Construct BWT for the packed sequence...

# Terminal 2 output (slightly later):
[timestamp] rule process_reference:
    ...
Another process is building the reference. Waiting...
Reference completed by another process.
```

### **Serial Processing Within Each Chunk:**
```
# Each terminal will show:
Running serial alignment for 2 samples...
Using FASTQ input mode...
Processing FASTQ sample: SRR1564422
Processing FASTQ sample: SRR1564423  # (or SRR1564424, SRR1564425 for chunk 2)

Calling insertions for 2 samples serially...
Processing sample SRR1564422...
Processing sample SRR1564423...
```

---

## ✅ **Step 5: Verify Stage 1 Completion**

After both jobs complete, verify all outputs:

```bash
echo "=== Stage 1 Outputs ==="

# Check BAM files
echo "BAM files created:"
ls -la data/horse/*.bam

# Check VCF files  
echo "VCF files created:"
find output/horse/ -name "out.pass.vcf.gz" -ls

# Count total samples processed
echo "Total samples processed: $(find output/horse/ -name "out.pass.vcf.gz" | wc -l)"
```

**Expected Output:**
```
BAM files created:
-rw-r--r-- 1 user staff 4068933664 Dec  7 data/horse/SRR1564422.bam
-rw-r--r-- 1 user staff 4068933664 Dec  7 data/horse/SRR1564423.bam  
-rw-r--r-- 1 user staff 4068933664 Dec  7 data/horse/SRR1564424.bam
-rw-r--r-- 1 user staff 4068933664 Dec  7 data/horse/SRR1564425.bam

VCF files created:
output/horse/SRR1564422/out.pass.vcf.gz
output/horse/SRR1564423/out.pass.vcf.gz
output/horse/SRR1564424/out.pass.vcf.gz  
output/horse/SRR1564425/out.pass.vcf.gz

Total samples processed: 4
```

---

## 🔬 **Step 6: Run Stage 2 - Automatic Sample Detection**

Stage 2 automatically detects all completed samples from Stage 1:

```bash
echo "=== Running Stage 2 ==="

# Restore the full sample list for Stage 2
cp tutorial_horse_samples.txt tutorial_horse_samples.txt.backup

# Run Stage 2 (automatically detects all completed samples)
snakemake -s Snakefile_stage2.py \
  --configfile tutorial_config.yaml \
  --cores 8 \
  --dry-run

echo "Dry-run completed. Now running actual Stage 2:"

snakemake -s Snakefile_stage2.py \
  --configfile tutorial_config.yaml \
  --cores 8
```

**Stage 2 will show automatic detection:**
```
Stage 2: Auto-detected 4 completed samples from VCF files
Detected samples: SRR1564422, SRR1564423, SRR1564424, SRR1564425
Building DAG of jobs...
```

---

## 📊 **Step 7: Examine Final Results**

```bash
echo "=== Final Pipeline Outputs ==="

# Core analysis files
echo "Global VCF and analysis:"
ls -la output/horse/horse_global_vcf.txt
ls -la output/horse/tepeak_families_annotated.csv
ls -la output/horse/tepeak_clusters_plot.png
ls -la output/horse/horse_insertions_plot.png

# Gene annotation results  
echo "Gene annotation results:"
ls -la output/horse/peak_200-6400/

# Sample tracking files
echo "Internal pipeline files:"
ls -la output/horse/.stage1_complete
ls -la output/horse/.stage2_samples.txt

# Count insertions found
echo "Total insertions found:"
grep -v "^#" output/horse/horse_global_vcf.txt | wc -l
```

---

## 📚 **Key Learning Points**

### **1. Resource Efficiency**
- Only **one** reference gets built despite running concurrent jobs
- **Serial processing** within chunks prevents memory exhaustion  
- **Lock files** coordinate concurrent access safely

### **2. Fault Tolerance**
```bash
# If a chunk fails, you can rerun just that chunk:
cp chunks/horse_chunk01.txt tutorial_horse_samples.txt
snakemake -s Snakefile_stage1.py --configfile tutorial_config.yaml --cores 4
```

### **3. Scalability**
```bash
# For larger datasets, simply increase chunk size:
python scripts/split_samples.py large_sample_list.txt chunks/samples --chunk-size 20

# Or create more chunks:
python scripts/split_samples.py large_sample_list.txt chunks/samples --chunk-size 5  # Creates more, smaller chunks
```

### **4. HPC Compatibility**
```bash
# Each chunk can be submitted as a separate job:
for chunk in chunks/horse_chunk*.txt; do
    sbatch --job-name="tepeak_$(basename $chunk)" \
           --cpus-per-task=32 \
           --mem=64G \
           run_stage1_chunk.sh $chunk
done
```

---

## 🧹 **Cleanup**

```bash
echo "=== Cleaning Up Tutorial Files ==="

# Remove tutorial-specific files
rm -f tutorial_horse_samples.txt tutorial_horse_samples_chunk2.txt
rm -f tutorial_config.yaml tutorial_config_chunk2.yaml  
rm -rf chunks/horse_chunk*.txt

# Optionally remove outputs (keep if you want to examine them)
# rm -rf output/horse/
# rm -rf data/horse/SRR156442[3-5]*

echo "Tutorial cleanup complete!"
```

---

## 🎯 **Production Recommendations**

### **For Small Datasets (< 20 samples):**
```bash
# Chunk size: 5-10 samples
python scripts/split_samples.py samples.txt chunks/samples --chunk-size 8
```

### **For Medium Datasets (20-100 samples):**  
```bash
# Chunk size: 10-20 samples
python scripts/split_samples.py samples.txt chunks/samples --chunk-size 15
```

### **For Large Datasets (100+ samples):**
```bash  
# Chunk size: 20-50 samples
python scripts/split_samples.py samples.txt chunks/samples --chunk-size 25
```

### **Memory Considerations:**
- **Per chunk**: ~8GB RAM per sample during alignment
- **Chunk size calculation**: `Available RAM / 8GB = Max samples per chunk`
- **Disk space**: ~10GB temporary space per sample during processing

---

## 🆘 **Troubleshooting**

### **Issue: "Reference building lock timeout"**
```bash
# Remove stale lock file
rm -f data/*/reference_building.lock
```

### **Issue: "Sample not found in Stage 2"**
```bash
# Verify VCF files exist
find output/ -name "out.pass.vcf.gz" -ls

# Check Stage 2 sample detection
cat output/*/stage2_samples.txt
```

### **Issue: "Disk space full during alignment"**
```bash
# Reduce chunk size and ensure adequate /tmp space
python scripts/split_samples.py samples.txt chunks/samples --chunk-size 5
df -h /tmp  # Check available space
```

---

This tutorial demonstrates the complete staged pipeline workflow, showing how the new serial execution approach prevents resource conflicts while maintaining efficient parallel processing at the chunk level.