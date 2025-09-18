# TEPEAK
A pipeline for identifying and characterizing polymorphic transposable elements in non-model species populations.

---

## Setup

We recommend running everything inside a dedicated `conda` environment.

### Create the environment
```bash
conda env create -f environment.yaml
conda activate insurveyor-env
```

### Additional requirements
TEPEAK requires a few tools that are not bundled into the conda environment:

#### Picard
```bash
git clone https://github.com/broadinstitute/picard.git
cd picard/
./gradlew shadowJar
java -jar build/libs/picard.jar
```

If you run into Java errors, make sure you have a modern Java JDK installed.

#### Smoove

We run Smoove entirely via Docker. You must have Docker installed. It is not currently possible to run via the binaries. We will implement this feature in the future. 

```bash
docker pull brentp/smoove
docker run -it brentp/smoove smoove -h
```

## Data requirements

TEPEAK can run in two modes:

1. NCBI SRA mode (default)

Inputs:

A zipped genome reference (e.g. `ncbi_dataset.zip`)

A zipped genome GTF (for gene annotations, e.g. `ncbi_dataset_gtf.zip`)

An SRA run table (CSV/TSV from NCBI)

Example project structure:

```
TEPEAK/
  Snakefile
  config.yaml
  SraRunTable.txt
  data/
    ncbi_dataset.zip
    ncbi_dataset_gtf.zip

```


2. Custom FASTQ mode

Inputs:

- A zipped genome reference

- A zipped genome GTF

- A directory of FASTQs

- A plain text sample list (one sample ID per line)

FASTQs must be paired-end, named with any of the following conventions:

- <sample>_R1.fastq(.gz) / <sample>_R2.fastq(.gz)

- <sample>_1.fastq(.gz) / <sample>_2.fastq(.gz)

- <sample>.R1.fastq(.gz) / <sample>.R2.fastq(.gz)

- <sample>-R1.fastq(.gz) / <sample>-R2.fastq(.gz)

```
TEPEAK/
  Snakefile
  config_custom.yaml
  custom_samples.txt
  data/
    my_species_reference.zip
    my_species_gtf.zip
    my_species/
      fastq/
        sample1_R1.fastq.gz
        sample1_R2.fastq.gz
        sample2_1.fastq.gz
        sample2_2.fastq.gz

```

### Configuration

All runs are controlled by a config.yaml file.
You must specify whether you are using SRA or FASTQ input via the input_type field.

Example: SRA input
```
species: ecoli
data_dir: data
output_dir: output
zipped_ref_genome_filepath: data/ncbi_dataset.zip
zipped_gtf_filepath: data/ncbi_dataset_gtf.zip
sra_run: 
  filepath: SraRunTable.txt
  number_of_runs: 4
input_type: sra
threads: 10
low: 0
high: 10000
gene: y
```

Example: custom FASTQ Input

```
species: horse
data_dir: data
output_dir: output
zipped_ref_genome_filepath: data/horse_reference.zip
zipped_gtf_filepath: data/horse_gtf.zip
fastq_input:
  sample_list: horse_samples.txt
  fastq_dir: data/horse/fastq
input_type: fastq
threads: 8
low: 0
high: 10000
gene: n

```

## Running TEPEAK

From the TEPEAK directory with the environment activated:
```
snakemake --configfile config.yaml --cores 8
```

## Output

All results are written to <output_dir>/<species>/.

Key outputs include:

- <species>_insertions_plot.svg — insertion histogram

- <species>_smoove_plot.svg — deletion histogram

- <species>_global_vcf.txt — global VCF of insertions

- <species>.smoove.square.vcf.gz — squared Smoove VCF across samples

Gene annotation tables (if gene: y)

## Notes

- Always activate the `insurveyor-env` environment before running.

- For large references, indexing (`bwa index`) may take hours.

- Delete `prefetch_tmp` when finished to save space.

- If running on Windows, use WSL
