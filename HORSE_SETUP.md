# Horse Analysis Setup Instructions

## Current Setup Status

The directory structure has been created for horse analysis. You'll need to complete the following steps:

### 1. Move your horse FASTQ files to the correct location

```bash
# Move your horse FASTQ files to the correct directory with proper naming
mv SRR1564422.fastq.gz data/horse/fastq/horse1_R1.fastq.gz

# If you have the paired-end file (R2), also move it:
# mv SRR1564422_2.fastq.gz data/horse/fastq/horse1_R2.fastq.gz

# Or if you have a different naming convention, rename to match our standard:
# mv your_horse_file_1.fastq.gz data/horse/fastq/horse1_R1.fastq.gz
# mv your_horse_file_2.fastq.gz data/horse/fastq/horse1_R2.fastq.gz
```

### 2. Add the horse reference genome

You mentioned you have a horse reference genome. Place it in the data directory:

```bash
# Move your horse reference genome zip file to:
mv horse_reference_genome.zip data/horse_reference.zip

# Update the config file to point to the correct reference
# The config_horse_custom.yaml is already created and ready to use
```

### 3. Update the configuration

The file `config_horse_custom.yaml` has been created with these settings:

```yaml
species: horse
data_dir: data
output_dir: output
zipped_ref_genome_filepath: data/horse_reference.zip
zipped_gtf_filepath: data/horse_reference.zip  
input_type: fastq
fastq_input:
  sample_list: horse_samples.txt
  fastq_dir: data/horse/fastq
threads: 8
low: 0
high: 10000
gene: y
```

### 4. Verify the setup

Check that your directory structure looks like this:

```
TEPEAK/
├── config_horse_custom.yaml
├── horse_samples.txt
└── data/
    ├── horse_reference.zip          # Horse reference genome + GTF
    └── horse/
        └── fastq/
            ├── horse1_R1.fastq.gz    # Your horse FASTQ R1
            └── horse1_R2.fastq.gz    # Your horse FASTQ R2
```

### 5. Run the analysis

Once everything is in place:

```bash
conda activate insurveyor-env
snakemake --configfile config_horse_custom.yaml --cores 8 --latency-wait 60
```

## Notes

- The sample name "horse1" is used in `horse_samples.txt` - you can change this if needed
- Make sure your FASTQ files are paired-end and properly named
- The reference genome zip should contain both the genome FASTA and GTF annotation files
- Adjust the `threads` parameter in the config based on your system capabilities