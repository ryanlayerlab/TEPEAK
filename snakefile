shell.prefix(
    'fastq-dump > /dev/null 2>&1 || export PATH=$PATH:$PWD/$(ls | grep "sratoolkit")/bin;'
)

configfile: 'config.yaml'
import pandas as pd

REF_EXTENSIONS = ['fa', 'fa.fai', 'amb', 'ann', 'bwt', 'pac', 'sa']
ANN_EXTENSIONS = ['gene_annotate.txt', 'gtf_loci.txt', 'gtf.txt', 'pop_vcf_sorted.txt']

species = config['species']
data_dir = config['data_dir'].strip('/')
threads = config['threads']
low, high = config['low'], config['high']

species_dir = f'{data_dir}/{species}'
output_dir = f"{config['output_dir'].strip('/')}/{species}"

final = f'{output_dir}/peak_{low}-{high}/{species}_{low}-{high}_merged'
ALL = [
    f'{output_dir}/{species}_global_vcf.txt', 
    f'{output_dir}/dfam_annotate.csv', 
    f'{output_dir}/{species}_insertions_plot.svg', 
    f'{output_dir}/{species}_smoove_plot.svg', 
    f'{final}_genes.txt'if config['gene'].lower() in ('y', 'yes') else f'{final}.txt'
]

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
sample_file = f'{species_dir}/{species}_samples.txt'
if not os.path.exists(sample_file): 
    get_sra_numbers(
        config['sra_run']['filepath'], 
        sample_file, 
        config['sra_run']['number_of_runs']
    )
SAMPLES = get_samples(sample_file)

rule all:  
    input: 
        ALL

rule process_reference: 
    input: 
        sample_file = sample_file, 
        ref = config['zipped_ref_genome_filepath']
    params: 
        species_dir = species_dir, 
        species = species, 
    output: 
        ref = expand(f'{species_dir}/{species}.{{ext}}', ext = REF_EXTENSIONS)
    script: 
        "src/process_reference.py"

rule align_species: 
    input: 
        ref_files = rules.process_reference.output.ref, 
        sample_file = sample_file
    params: 
        species = species, 
        species_dir = species_dir, 
        threads = config['threads']
    threads: config['threads']
    output: 
        expand(
            [
                f'{species_dir}/{{sample}}.bam', 
                f'{species_dir}/{{sample}}.bam.bai', 
            ], 
            sample = SAMPLES
        )
    shell: 
        "bash src/align_species.sh -s {params.species} -d {params.species_dir} -t {params.threads} -f {input.sample_file}"

rule call_insertions_serial: 
    input: 
        sample_file = sample_file, 
        ref = f'{species_dir}/{species}.fa', 
        bam_file = expand(f'{species_dir}/{{sample}}.bam', sample = SAMPLES)
    params: 
        threads = threads, 
        species_dir = species_dir, 
        output_dir = output_dir
    output: 
        vcf_file = expand(f'{output_dir}/{{sample}}/out.pass.vcf.gz', sample = SAMPLES)
    script: 
        "src/call_insertions_serial.py"

rule check_insertions: 
    input: 
        sample_file = sample_file, 
        vcf_file = rules.call_insertions_serial.output.vcf_file
    params: 
        species = species, 
        output_dir = output_dir
    output: 
        count_file = f'{output_dir}/{species}_count.txt', 
        vcf_file = expand(f'{output_dir}/{{sample}}/out.pass.vcf', sample = SAMPLES)
    shell:
        "bash src/check_insertions.sh -s {params.species} -o {params.output_dir} -f {input.sample_file}"

rule get_global_vcf: 
    input: 
        sample_file = sample_file, 
        decompressed_file = rules.check_insertions.output.vcf_file
    params: 
        species = species, 
        output_dir = output_dir
    output: 
        global_vcf = f'{output_dir}/{species}_global_vcf.txt'
    shell: 
        "bash src/get_global_vcf.sh -s {params.species} -o {params.output_dir} -f {input.sample_file}"

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
        vcf_file = rules.check_insertions.output.vcf_file
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
        plot = f'{output_dir}/{species}_insertions_plot.svg'
    script: 
        "src/build_histogram.py"

rule smoove_1: 
    input: 
        bam_files = f'{species_dir}/{{sample}}.bam',
        smoove_file = 'smoove_latest.sif',
        reference_file = f'{species_dir}/{species}.fa', 
        input_dir = f'{species_dir}'
    params: 
        output_dir = f'{output_dir}/results-smoove/'
    output:
        results_smoove = f'{output_dir}/results-smoove/{{sample}}-smoove.genotyped.vcf.gz'
    shell: 
        "singularity exec {input.smoove_file} smoove call --outdir {params.output_dir} --name {wildcards.sample} "
        "--fasta {input.reference_file} -p 1 --genotype {input.input_dir}/{wildcards.sample}.bam"

rule smoove_2: 
    input: 
        smoove_file = 'smoove_latest.sif', 
        reference_file = f'{species_dir}/{species}.fa',
	    results = expand(rules.smoove_1.output.results_smoove, sample = SAMPLES)
    params:
        prev_output_dir = rules.smoove_1.params.output_dir,
        output_dir = output_dir
    output: 
        merged_sites_file = f'{output_dir}/merged.sites.vcf.gz'
    shell: 
        "singularity exec {input.smoove_file} smoove merge --name merged -f {input.reference_file} --outdir {params.output_dir} " 
        "{params.prev_output_dir}*.genotyped.vcf.gz"

rule smoove_3: 
    input: 
        smoove_file = 'smoove_latest.sif', 
        reference_file = f'{species_dir}/{species}.fa',
        merged_file = rules.smoove_2.output.merged_sites_file, 
        bam_files = f'{species_dir}/{{sample}}.bam'
    params: 
        output_dir = f'{output_dir}/results-genotyped/',
        species_dir = species_dir,
	    samples = SAMPLES
    output: 
        results_genotyped = f'{output_dir}/results-genotyped/{{sample}}-joint-smoove.genotyped.vcf.gz'
    shell: 
        "singularity exec {input.smoove_file} smoove genotype -d -x -p 1 --name {wildcards.sample}-joint "
        "--outdir {params.output_dir} --fasta {input.reference_file} --vcf {input.merged_file} "
        "{params.species_dir}/{wildcards.sample}.bam"

rule smoove_4: 
    input: 
        smoove_file = 'smoove_latest.sif',
        vcf_files = expand(rules.smoove_3.output.results_genotyped, sample = SAMPLES)
    params: 
        species = species
    output: 
        f'{species}.smoove.square.vcf.gz'
    shell: 
        "singularity exec {input.smoove_file} smoove paste --name {params.species} {input.vcf_files}"

rule smoove_global_vcf: 
    input: 
        f'{species}.smoove.square.vcf.gz'
    params: 
        input_file = f'{species}.smoove.square.vcf'
    output:
        global_vcf = f'{species}_smoove_global.vcf'
    shell: 
        "gzip -d {input};"
        """bcftools query -f "%CHROM\t%POS\t%INFO/END\t%SVLEN\n" {params.input_file} >> {output} """


rule build_histogram_deletions: 
    input: 
        global_vcf_file = rules.smoove_global_vcf.output.global_vcf
    params: 
        low = low, 
        high = high
    output: 
        plot = f'{output_dir}/{species}_smoove_plot.svg'
    script: 
        "src/build_histogram_deletions.py"