"""
TEPEAK Stage 1: Per-sample alignment and insertion calling
This stage can be run on sample subsets to manage resource usage
"""

shell.prefix(
    'source ~/.bashrc 2>/dev/null || true; '
    'conda activate insurveyor-env 2>/dev/null || true; '
    'export PATH="$CONDA_PREFIX/bin:$PATH"; '
    'fastq-dump > /dev/null 2>&1 || export PATH=$PATH:$PWD/$(ls | grep "sratoolkit")/bin;'
)

configfile: 'config.yaml'
import pandas as pd
import os

REF_EXTENSIONS = ['fa', 'fa.fai', 'fa.amb', 'fa.ann', 'fa.bwt', 'fa.pac', 'fa.sa']

species = config['species']
data_dir = config['data_dir'].strip('/')
threads = config['threads']

species_dir = f'{data_dir}/{species}'
output_dir = f"{config['output_dir'].strip('/')}/{species}"

# Handle both SRA and custom FASTQ input modes
input_type = config.get('input_type', 'sra')

if input_type == 'fastq':
    sample_file = config['fastq_input']['sample_list']
else:
    sample_file = f'{species_dir}/{species}_samples.txt'

def get_samples(filename: str) -> list:
    if not os.path.exists(filename):
        print(f"Warning: Sample file {filename} not found")
        return []
    df = pd.read_csv(filename, names=['samples'])
    return list(df['samples'])

SAMPLES = get_samples(sample_file)
print(f"Stage 1: Processing {len(SAMPLES)} samples")

# Stage 1 outputs: BAM files and VCF files for each sample
STAGE1_ALL = []
if SAMPLES:
    STAGE1_ALL.extend([
        expand(f'{species_dir}/{{sample}}.bam', sample=SAMPLES),
        expand(f'{species_dir}/{{sample}}.bam.bai', sample=SAMPLES),
        expand(f'{output_dir}/{{sample}}/out.pass.vcf.gz', sample=SAMPLES)
    ])

rule all:
    input: STAGE1_ALL

rule process_reference:
    input:
        ref = config['zipped_ref_genome_filepath']
    params:
        species_dir = species_dir,
        species = species
    output:
        ref = expand(f'{species_dir}/{species}.{{ext}}', ext=REF_EXTENSIONS)
    script:
        "src/process_reference.py"

rule align_one_sample:
    input:
        ref_files = rules.process_reference.output.ref,
        sample_file = sample_file
    params:
        species = species,
        species_dir = species_dir,
        input_type = input_type,
        fastq_dir = config.get('fastq_input', {}).get('fastq_dir', '') if input_type == 'fastq' else ''
    threads: max(1, config['threads'] // max(1, len(SAMPLES)))
    output:
        bam = f'{species_dir}/{{sample}}.bam',
        bai = f'{species_dir}/{{sample}}.bam.bai'
    shell:
        """
        if [ "{params.input_type}" = "fastq" ]; then
            echo "Aligning sample {wildcards.sample} with FASTQ input..."
            echo "{wildcards.sample}" > {params.species_dir}/{wildcards.sample}_temp.txt
            bash src/align_species_fastq.sh \
                -s {params.species} \
                -d {params.species_dir} \
                -t {threads} \
                -f {params.fastq_dir} \
                -l {params.species_dir}/{wildcards.sample}_temp.txt
            rm -f {params.species_dir}/{wildcards.sample}_temp.txt
        else
            echo "Aligning sample {wildcards.sample} with SRA input..."
            echo "{wildcards.sample}" > {params.species_dir}/{wildcards.sample}_temp.txt
            bash src/align_species_sra.sh \
                -s {params.species} \
                -d {params.species_dir} \
                -t {threads} \
                -f {params.species_dir}/{wildcards.sample}_temp.txt
            rm -f {params.species_dir}/{wildcards.sample}_temp.txt
        fi
        """

rule align_serial:
    input:
        ref_files = rules.process_reference.output.ref,
        sample_file = sample_file
    params:
        species = species,
        species_dir = species_dir,
        input_type = input_type,
        fastq_dir = config.get('fastq_input', {}).get('fastq_dir', '') if input_type == 'fastq' else '',
        samples = SAMPLES
    threads: threads
    output:
        bam_files = expand(f'{species_dir}/{{sample}}.bam', sample=SAMPLES),
        bai_files = expand(f'{species_dir}/{{sample}}.bam.bai', sample=SAMPLES)
    shell:
        """
        echo "Running serial alignment for $(wc -l < {input.sample_file}) samples..."
        
        # Use the correct batch alignment script from scripts/ directory
        # This script is designed to process multiple samples from a file
        SCRIPT_PATH="$PWD/scripts/align_species.sh"
        
        if [ ! -f "$SCRIPT_PATH" ]; then
            echo "Error: align_species.sh script not found at $SCRIPT_PATH"
            echo "Available scripts:"
            ls -la "$PWD/scripts/" || echo "scripts/ directory not found"
            echo "Current working directory: $PWD"
            exit 1
        fi
        
        if [ "{params.input_type}" = "fastq" ]; then
            echo "Using FASTQ input mode..."
            bash "$SCRIPT_PATH" \
                -s {params.species} \
                -d {params.species_dir} \
                -t {threads} \
                -f {params.fastq_dir} \
                -l {input.sample_file}
        else
            echo "Using SRA input mode..."
            bash "$SCRIPT_PATH" \
                -s {params.species} \
                -d {params.species_dir} \
                -t {threads} \
                -l {input.sample_file}
        fi
        
        echo "Alignment script completed. Verifying outputs..."
        
        # List what files were actually created
        echo "Files created in {params.species_dir}:"
        ls -la {params.species_dir}/*.bam* 2>/dev/null || echo "No BAM files found"
        
        # Verify all expected output files exist
        missing_files=()
        for sample in {params.samples}; do
            bam_file="{params.species_dir}/$sample.bam"
            bai_file="{params.species_dir}/$sample.bam.bai"
            
            if [ ! -f "$bam_file" ]; then
                missing_files+=("$bam_file")
            fi
            
            if [ ! -f "$bai_file" ]; then
                missing_files+=("$bai_file")
            fi
        done
        
        if [ ${{#missing_files[@]}} -gt 0 ]; then
            echo "Error: Missing expected output files:"
            printf '  %s\n' "${{missing_files[@]}}"
            echo ""
            echo "This might be due to:"
            echo "1. Sample names in the sample file don't match expected format"
            echo "2. Alignment script failed for some samples"
            echo "3. File system latency (try increasing --latency-wait)"
            exit 1
        fi
        
        echo "All alignment outputs verified successfully"
        """

rule call_insertions_one_sample:
    input:
        ref = f'{species_dir}/{species}.fa',
        bam = f'{species_dir}/{{sample}}.bam',
        bai = f'{species_dir}/{{sample}}.bam.bai'
    params:
        output_dir = output_dir
    threads: 1
    output:
        vcf = f'{output_dir}/{{sample}}/out.pass.vcf.gz'
    shell:
        """
        # Use the correct insurveyor.py command
        INSURVEYOR_PATH="$CONDA_PREFIX/bin/insurveyor.py"
        
        if [ ! -f "$INSURVEYOR_PATH" ]; then
            echo "Error: insurveyor.py not found at $INSURVEYOR_PATH"
            echo "CONDA_PREFIX: $CONDA_PREFIX"
            exit 1
        fi
        
        mkdir -p {params.output_dir}/{wildcards.sample}
        echo "Calling insertions for sample {wildcards.sample}..."
        echo "Using insurveyor at: $INSURVEYOR_PATH"
        
        python "$INSURVEYOR_PATH" {input.bam} {params.output_dir}/{wildcards.sample} {input.ref}
        """
