"""
TEPEAK Stage 2: Merging, analysis, and annotation
This stage requires all VCF files from Stage 1 to be complete
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
import hashlib
import json
import glob

ANN_EXTENSIONS = ['gene_annotate.txt', 'gtf_loci.txt', 'gtf.txt', 'pop_vcf_sorted.txt']

species = config['species']
data_dir = config['data_dir'].strip('/')
threads = config['threads']
low, high = config['low'], config['high']
run_smoove = config.get('run_smoove', False)

species_dir = f'{data_dir}/{species}'
output_dir = f"{config['output_dir'].strip('/')}/{species}"

final = f'{output_dir}/peak_{low}-{high}/{species}_{low}-{high}_merged'

# Auto-detect completed samples from Stage 1 VCF files
def get_completed_samples():
    vcf_pattern = f'{output_dir}/*/out.pass.vcf.gz'
    vcf_files = glob.glob(vcf_pattern)
    samples = []
    for vcf_file in vcf_files:
        # Extract sample name from path: output_dir/SAMPLE/out.pass.vcf.gz
        sample = os.path.basename(os.path.dirname(vcf_file))
        samples.append(sample)
    return sorted(samples)  # Sort for reproducible behavior

SAMPLES = get_completed_samples()
print(f"Stage 2: Auto-detected {len(SAMPLES)} completed samples from VCF files")

if SAMPLES:
    print("Detected samples:", ', '.join(SAMPLES[:10]) + ('...' if len(SAMPLES) > 10 else ''))

# Handle both SRA and custom FASTQ input modes for sample file
input_type = config.get('input_type', 'sra')
if input_type == 'fastq':
    sample_file = config['fastq_input']['sample_list']
else:
    sample_file = f'{species_dir}/{species}_samples.txt'

# Create a stage 2 sample file with detected samples
stage2_sample_file = f'{output_dir}/.stage2_samples.txt'

# Stage 2 outputs
STAGE2_ALL = [
    f'{output_dir}/{species}_global_vcf.txt', 
    f'{output_dir}/tepeak_families_annotated.csv',
    f'{output_dir}/tepeak_clusters_plot.png',
    f'{output_dir}/{species}_insertions_plot.png', 
    f'{final}_genes.txt' if config.get('gene', 'n').lower() in ('y', 'yes') else f'{final}.txt'
]

# Add optional outputs
if config.get('gene', '').lower() in ('y', 'yes') and config.get('run_enrichment', False):
    STAGE2_ALL.extend([
        f'{output_dir}/{species}_enrichment_results.csv',
        f'{output_dir}/{species}_enrichment_summary.txt'
    ])

if config.get('run_phylogeny', False):
    STAGE2_ALL.extend([
        f'{output_dir}/{species}_te_alignment.fasta',
        f'{output_dir}/{species}_te_phylogeny.newick',
        f'{output_dir}/{species}_phylogeny_summary.txt'
    ])

if run_smoove:
    STAGE2_ALL.extend([
        f'{output_dir}/{species}_smoove_plot.png'
    ])

rule all:
    input: STAGE2_ALL

# Create sample list file for Stage 2 processing
rule create_stage2_sample_file:
    output:
        sample_file = stage2_sample_file
    run:
        if not SAMPLES:
            raise ValueError("No completed samples detected from Stage 1. Ensure VCF files exist in output/{species}/*/out.pass.vcf.gz")
        
        with open(output.sample_file, 'w') as f:
            for sample in SAMPLES:
                f.write(f"{sample}\n")
        
        print(f"Created Stage 2 sample file with {len(SAMPLES)} samples")

# Verify all Stage 1 VCF files are present
rule verify_stage1_complete:
    input:
        vcf_files = expand(f'{output_dir}/{{sample}}/out.pass.vcf.gz', sample=SAMPLES) if SAMPLES else [],
        sample_file = rules.create_stage2_sample_file.output.sample_file
    output:
        marker = f'{output_dir}/.stage1_complete'
    params:
        num_samples = len(SAMPLES) if SAMPLES else 0,
        output_dir = output_dir
    shell:
        """
        echo "Verifying Stage 1 completion for {params.num_samples} samples..."
        
        if [ {params.num_samples} -eq 0 ]; then
            echo "Error: No samples detected from Stage 1"
            echo "Expected VCF files in: {params.output_dir}/*/out.pass.vcf.gz"
            exit 1
        fi
        
        echo "Found VCF files for samples:"
        for vcf in {input.vcf_files}; do
            echo "  $vcf"
            if [ ! -f "$vcf" ]; then
                echo "Error: Missing VCF file: $vcf"
                exit 1
            fi
        done
        echo "All Stage 1 VCF files verified successfully"
        touch {output.marker}
        """

rule get_global_vcf:
    input:
        marker = rules.verify_stage1_complete.output.marker,
        sample_file = stage2_sample_file
    params:
        species = species,
        output_dir = output_dir
    output:
        global_vcf = f'{output_dir}/{species}_global_vcf.txt'
    shell:
        "bash src/get_global_vcf.sh -s {params.species} -o {params.output_dir} -f {input.sample_file}"

rule tepeak_advanced: 
    input: 
        global_vcf = rules.get_global_vcf.output.global_vcf
    params: 
        output_dir = output_dir
    output: 
        families_csv = f'{output_dir}/tepeak_families_annotated.csv',
        clusters_plot = f'{output_dir}/tepeak_clusters_plot.png'
    script: 
        "src/tepeak_advanced.py"

rule extract_range:
    input:
        marker = rules.verify_stage1_complete.output.marker,
        sample_file = stage2_sample_file
    params:
        species = species,
        output_dir = output_dir,
        low = config['low'],
        high = config['high']
    output:
        range_file = f'{output_dir}/peak_{low}-{high}/{species}_{low}-{high}_pop_vcf.txt'
    shell:
        "bash src/extract_range.sh -s {params.species} -f {input.sample_file} -o {params.output_dir} -l {params.low} -h {params.high}"

rule get_species_gtf: 
    input: 
        zipped_gtf_file = config['zipped_gtf_filepath']
    params: 
        species = species, 
        species_dir = species_dir
    output: 
        gtf_file = f'{species_dir}/{species}.gtf'
    shell: 
        "bash src/get_species_gtf.sh -s {params.species} -d {params.species_dir} -f {input.zipped_gtf_file}"

rule annotate_genes: 
    input: 
        gtf_file = rules.get_species_gtf.output.gtf_file, 
        range_file = rules.extract_range.output.range_file
    params: 
        species = species, 
        output_dir = output_dir, 
        low = low, 
        high = high
    output:
        expand(f'{output_dir}/peak_{low}-{high}/{species}_{low}-{high}_{{ext}}', ext = ANN_EXTENSIONS)
    script: 
        "src/gene_helper.py"

rule write_output: 
    input:
        range_file = rules.extract_range.output.range_file
    params: 
        species = species,
        output_dir = output_dir,
        low = low, 
        high = high
    output: 
        f'{output_dir}/peak_{low}-{high}/{species}_{low}-{high}_pop_vcf_sorted.txt', 
        f'{output_dir}/peak_{low}-{high}/{species}_{low}-{high}_merged.txt'
    shell: 
        "bash src/write_output.sh -s {params.species} -o {params.output_dir} -r {input.range_file} -l {params.low} -h {params.high}"

rule write_output_genes:
    input: 
        annotated_file = f'{output_dir}/peak_{low}-{high}/{species}_{low}-{high}_gene_annotate.txt'
    params: 
        species = species, 
        output_dir = output_dir, 
        low = low, 
        high = high
    output: 
        f'{output_dir}/peak_{low}-{high}/{species}_{low}-{high}_merged_genes.txt'
    script: 
        "src/output_helper.py"

rule build_histogram_insertions: 
    input: 
        global_vcf_file = f'{output_dir}/{species}_global_vcf.txt'
    params: 
        low = low, 
        high = high
    output: 
        plot = f'{output_dir}/{species}_insertions_plot.png'
    script: 
        "src/build_histogram.py"

# Add a config checkpoint to trigger reruns when config changes
def get_config_hash():
    """Generate hash of relevant config parameters to detect changes."""
    relevant_config = {
        'run_enrichment': config.get('run_enrichment', False),
        'gene': config.get('gene', 'n'),
        'enrichment_organism': config.get('enrichment_organism', ''),
        'enrichment_sources': config.get('enrichment_sources', []),
        'enrichment_max_pvalue': config.get('enrichment_max_pvalue', 0.05)
    }
    config_str = json.dumps(relevant_config, sort_keys=True)
    return hashlib.md5(config_str.encode()).hexdigest()

CONFIG_HASH = get_config_hash()

# Only define enrichment rule if both gene annotation and enrichment are enabled
if config.get('gene', '').lower() in ('y', 'yes') and config.get('run_enrichment', False):
    rule functional_enrichment:
        input:
            annotation_file = f'{output_dir}/peak_{low}-{high}/{species}_{low}-{high}_gene_annotate.txt'
        params:
            species = species,
            output_dir = output_dir,
            config_hash = CONFIG_HASH
        output:
            results_csv = f'{output_dir}/{species}_enrichment_results.csv',
            summary_txt = f'{output_dir}/{species}_enrichment_summary.txt'
        script:
            "src/functional_enrichment.py"

# Only define phylogeny rule if enabled
if config.get('run_phylogeny', False):
    rule phylogenetic_analysis:
        input:
            families_csv = f'{output_dir}/tepeak_families_annotated.csv'
        params:
            species = species,
            output_dir = output_dir
        threads: config.get('threads', 4)
        output:
            alignment = f'{output_dir}/{species}_te_alignment.fasta',
            tree = f'{output_dir}/{species}_te_phylogeny.newick',
            summary = f'{output_dir}/{species}_phylogeny_summary.txt'
        script:
            "src/phylogenetic_analysis.py"
