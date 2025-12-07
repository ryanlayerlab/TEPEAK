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
        
        if [ "{params.input_type}" = "fastq" ]; then
            echo "Using FASTQ input mode..."
            bash src/align_species_fastq.sh \
                -s {params.species} \
                -d {params.species_dir} \
                -t {threads} \
                -f {params.fastq_dir} \
                -l {input.sample_file}
        else
            echo "Using SRA input mode..."
            bash src/align_species_sra.sh \
                -s {params.species} \
                -d {params.species_dir} \
                -t {threads} \
                -f {input.sample_file}
        fi
        """

rule call_insertions_serial:
    input:
        ref = f'{species_dir}/{species}.fa',
        bam_files = expand(f'{species_dir}/{{sample}}.bam', sample=SAMPLES),
        bai_files = expand(f'{species_dir}/{{sample}}.bam.bai', sample=SAMPLES),
        sample_file = sample_file
    params:
        output_dir = output_dir,
        species_dir = species_dir
    threads: 1
    output:
        vcf_files = expand(f'{output_dir}/{{sample}}/out.pass.vcf.gz', sample=SAMPLES)
    shell:
        """
        # Use the correct insurveyor.py command
        INSURVEYOR_PATH="$CONDA_PREFIX/bin/insurveyor.py"
        
        if [ ! -f "$INSURVEYOR_PATH" ]; then
            echo "Error: insurveyor.py not found at $INSURVEYOR_PATH"
            echo "CONDA_PREFIX: $CONDA_PREFIX"
            exit 1
        fi
        
        echo "Calling insertions for $(wc -l < {input.sample_file}) samples serially..."
        echo "Using insurveyor at: $INSURVEYOR_PATH"
        
        # Process each sample serially to avoid memory issues
        while IFS= read -r sample; do
            echo "Processing sample $sample..."
            mkdir -p {params.output_dir}/$sample
            
            if [ -f "{params.species_dir}/$sample.bam" ]; then
                python "$INSURVEYOR_PATH" "{params.species_dir}/$sample.bam" "{params.output_dir}/$sample" {input.ref}
            else
                echo "Warning: BAM file for sample $sample not found"
            fi
        done < {input.sample_file}
        """
