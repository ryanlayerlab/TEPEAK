shell.prefix(
    'fastq-dump > /dev/null 2>&1 || export PATH=$PATH:$PWD/$(ls | grep "sratoolkit")/bin;'
)

configfile: 'config.yaml'
import pandas as pd
import os
import hashlib
import json

REF_EXTENSIONS = ['fa', 'fa.fai', 'fa.amb', 'fa.ann', 'fa.bwt', 'fa.pac', 'fa.sa']
ANN_EXTENSIONS = ['gene_annotate.txt', 'gtf_loci.txt', 'gtf.txt', 'pop_vcf_sorted.txt']

species = config['species']
data_dir = config['data_dir'].strip('/')
threads = config['threads']
low, high = config['low'], config['high']
run_smoove = config.get('run_smoove', False)

species_dir = f'{data_dir}/{species}'
output_dir = f"{config['output_dir'].strip('/')}/{species}"

final = f'{output_dir}/peak_{low}-{high}/{species}_{low}-{high}_merged'
ALL = [
    f'{output_dir}/{species}_global_vcf.txt', 
    f'{output_dir}/tepeak_families_annotated.csv',
    f'{output_dir}/tepeak_clusters_plot.png',
    f'{output_dir}/{species}_insertions_plot.png', 
    f'{final}_genes.txt' if config['gene'].lower() in ('y', 'yes') else f'{final}.txt'
]

# Add enrichment analysis outputs if enabled
if config.get('gene', '').lower() in ('y', 'yes') and config.get('run_enrichment', False):
    ALL.extend([
        f'{output_dir}/{species}_enrichment_results.csv',
        f'{output_dir}/{species}_enrichment_summary.txt'
    ])

# Add phylogenetic analysis outputs if enabled
if config.get('run_phylogeny', False):
    ALL.extend([
        f'{output_dir}/{species}_te_alignment.fasta',
        f'{output_dir}/{species}_te_phylogeny.newick',
        f'{output_dir}/{species}_phylogeny_summary.txt'
    ])

# Add smoove outputs only if run_smoove is True
if run_smoove:
    ALL.extend([
        f'{output_dir}/{species}_smoove_plot.png'
    ])


def get_samples(filename: str) -> list:
    df = pd.read_csv(filename, names = ['samples'])
    return list(df['samples'])

def get_sra_numbers(sra_file: str, sample_file: str, max_n: int) -> None:
    df = pd.read_csv(sra_file, low_memory = False)
    sorted_df = df.sort_values(by = 'Bases', ascending = False)
    max_n = len(sorted_df) if max_n is None else max_n
    top_N_runs = sorted_df.head(max_n)['Run']

    with open(sample_file, 'w') as f:
        for run in top_N_runs:
            f.write(f'{run}\n')

os.makedirs(species_dir, exist_ok = True)

# Handle both SRA and custom FASTQ input modes
input_type = config.get('input_type', 'sra')  # Default to 'sra' for backward compatibility

if input_type == 'fastq':
    # Custom FASTQ mode
    sample_file = config['fastq_input']['sample_list']
    if not os.path.exists(sample_file):
        raise FileNotFoundError(f"Sample list file not found: {sample_file}")
else:
    # Original SRA mode (default)
    sample_file = f'{species_dir}/{species}_samples.txt'
    if not os.path.exists(sample_file): 
        get_sra_numbers(
            config['sra_run']['filepath'], 
            sample_file, 
            config['sra_run']['number_of_runs']
        )

SAMPLES = get_samples(sample_file)
# number of samples available to smoove
N_SAMPLES = len(SAMPLES)

rule all:  
    input: 
        ALL

rule process_reference:
    input:
        ref = config['zipped_ref_genome_filepath']
    params:
        species_dir = species_dir,
        species = species
    output:
        ref = expand(f'{species_dir}/{species}.{{ext}}', ext = REF_EXTENSIONS)
    script:
        "src/process_reference.py"

# Per-sample alignment rule for parallelization
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

# Per-sample insertion calling
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
        mkdir -p {params.output_dir}/{wildcards.sample}
        echo "Calling insertions for sample {wildcards.sample}..."
        insurveyor {input.bam} {input.ref} {params.output_dir}/{wildcards.sample}
        """

# Aggregation rule - ensures all per-sample jobs complete before downstream rules
rule call_insertions_serial: 
    input: 
        sample_file = sample_file, 
        ref = f'{species_dir}/{species}.fa', 
        vcf_files = expand(f'{output_dir}/{{sample}}/out.pass.vcf.gz', sample=SAMPLES)
    params: 
        species_dir = species_dir, 
        output_dir = output_dir,
        samples = SAMPLES,
        num_samples = len(SAMPLES)  # Add this parameter
    output: 
        # Create a single aggregation marker file instead of the same files as input
        aggregation_complete = f'{output_dir}/.call_insertions_complete'
    shell:
        """
        echo "All per-sample insertion calling completed for {params.num_samples} samples"
        # Verify all files exist
        for sample in {params.samples}; do
            if [ ! -f {params.output_dir}/${{sample}}/out.pass.vcf.gz ]; then
                echo "Error: Missing VCF for sample ${{sample}}"
                exit 1
            fi
        done
        echo "All VCF files verified successfully"
        # Create completion marker
        touch {output.aggregation_complete}
        """

# Update downstream rules to depend on both the VCF files AND aggregation marker
rule check_insertions:
    input:
        sample_file = sample_file,
        vcf_files = expand(f'{output_dir}/{{sample}}/out.pass.vcf.gz', sample=SAMPLES),
        aggregation_complete = rules.call_insertions_serial.output.aggregation_complete
    params:
        species = species,
        output_dir = output_dir
    output:
        count_file = f'{output_dir}/{species}_count.txt'
    shell:
        "bash src/check_insertions.sh -s {params.species} -o {params.output_dir} -f {input.sample_file}"

rule get_global_vcf:
    input:
        sample_file = sample_file,
        vcf_files = expand(f'{output_dir}/{{sample}}/out.pass.vcf.gz', sample=SAMPLES),
        aggregation_complete = rules.call_insertions_serial.output.aggregation_complete
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

# Keep the original dfam_annotate rule for backward compatibility
rule dfam_annotate: 
    input: 
        global_vcf = rules.get_global_vcf.output.global_vcf, 
        ref = f'{species_dir}/{species}.fa'
    params: 
        output_dir = output_dir
    output: 
        f'{output_dir}/dfam_annotate.csv'
    script: 
        "src/dfam_annotate.py"

rule extract_range:
    input:
        sample_file = sample_file,
        vcf_files = expand(f'{output_dir}/{{sample}}/out.pass.vcf.gz', sample=SAMPLES),
        aggregation_complete = rules.call_insertions_serial.output.aggregation_complete
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
            config_hash = CONFIG_HASH  # This will trigger rerun when config changes
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

# Only define smoove rules if run_smoove is True
if run_smoove:
    rule smoove_1: 
        input: 
            bam_files = f'{species_dir}/{{sample}}.bam',
            bai_files = f'{species_dir}/{{sample}}.bam.bai',
            reference_file = f'{species_dir}/{species}.fa'
        params: 
            output_dir = f'{output_dir}/results-smoove/',
            species_dir = species_dir
        output:
            results_smoove = f'{output_dir}/results-smoove/{{sample}}-smoove.genotyped.vcf.gz'
        shell: 
            """
            mkdir -p {params.output_dir}
            smoove call --outdir {params.output_dir} --name {wildcards.sample} \
                --fasta {input.reference_file} -p 1 --genotype {input.bam_files}
            """

    rule smoove_2: 
        input: 
            reference_file = f'{species_dir}/{species}.fa',
            results = expand(rules.smoove_1.output.results_smoove, sample = SAMPLES)
        params:
            prev_output_dir = rules.smoove_1.params.output_dir,
            output_dir = output_dir
        output: 
            merged_sites_file = f'{output_dir}/merged.sites.vcf.gz'
        shell:
            """
            mkdir -p {params.output_dir}
            smoove merge --name merged -f {input.reference_file} \
                --outdir {params.output_dir} {params.prev_output_dir}/*.genotyped.vcf.gz
            """

    rule smoove_3: 
        input: 
            reference_file = f'{species_dir}/{species}.fa',
            merged_file = rules.smoove_2.output.merged_sites_file, 
            bam_files = f'{species_dir}/{{sample}}.bam'
        params: 
            output_dir = f'{output_dir}/results-genotyped/',
            species_dir = species_dir
        output: 
            results_genotyped = f'{output_dir}/results-genotyped/{{sample}}-joint-smoove.genotyped.vcf.gz'
        shell: 
            """
            mkdir -p {params.output_dir}
            smoove genotype -d -x -p 1 --name {wildcards.sample}-joint \
                --outdir {params.output_dir} --fasta {input.reference_file} \
                --vcf {input.merged_file} {input.bam_files}
            """

    rule smoove_4:
        input:
            vcf_files = expand(f'{output_dir}/results-genotyped/{{sample}}-joint-smoove.genotyped.vcf.gz',
                               sample=SAMPLES)
        params:
            species = species,
            n_samples = N_SAMPLES,
            genotyped_dir = f'{output_dir}/results-genotyped/',
            output_dir = output_dir
        output:
            f'{output_dir}/{species}.smoove.square.vcf.gz'
        shell:
            """
            set -euo pipefail
            
            if [ "{params.n_samples}" -eq 1 ]; then
                # Single-sample: copy and index
                cp -f {input.vcf_files} {output}
                if [ ! -s {output}.tbi ]; then
                    tabix -f -p vcf {output}
                fi
            else
                # Multi-sample: paste
                smoove paste --name {params.species} --outdir {params.output_dir} \
                    {params.genotyped_dir}/*.genotyped.vcf.gz
                # Ensure index exists
                if [ ! -s {output}.tbi ]; then
                    tabix -f -p vcf {output}
                fi
            fi
            """

    rule smoove_global_vcf:
        input:
            f'{output_dir}/{species}.smoove.square.vcf.gz'
        params:
            input_file = f'{output_dir}/{species}.smoove.square.vcf'
        output:
            global_vcf = f'{output_dir}/{species}_smoove_global.vcf'
        shell:
            "gzip -f -d {input} && "
            """bcftools query -f "%CHROM\t%POS\t%INFO/END\t%SVLEN\n" {params.input_file} > {output} """

    rule build_histogram_deletions: 
        input: 
            global_vcf_file = rules.smoove_global_vcf.output.global_vcf
        params: 
            low = low, 
            high = high
        output: 
            plot = f'{output_dir}/{species}_smoove_plot.png'
        script: 
            "src/build_histogram_deletions.py"
