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

**Important**: The Picard jar file should be located at `picard/build/libs/picard.jar` relative to the TEPEAK directory. The alignment scripts will automatically detect and use this path.

**Troubleshooting Picard Issues:**

If you encounter Picard path errors, the scripts will search for Picard in this order:
1. `TEPEAK/picard/build/libs/picard.jar` (recommended location)
2. `$(pwd)/picard/build/libs/picard.jar` (current directory)
3. `$PICARD_JAR_PATH` environment variable
4. Common system locations (`/usr/local/bin/picard.jar`, etc.)
5. `picard` command in PATH (conda or system install)

**Solutions for Picard issues:**

1. **Standard installation** (recommended):
   ```bash
   cd TEPEAK
   git clone https://github.com/broadinstitute/picard.git
   cd picard/
   ./gradlew shadowJar
   ```

2. **Custom location** - set environment variable:
   ```bash
   export PICARD_JAR_PATH="/path/to/your/picard.jar"
   ```

3. **System installation**:
   ```bash
   # Using conda (alternative)
   conda install -c bioconda picard
   
   # Or download pre-built JAR
   wget https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar
   ```

**Java Requirements:**
- Java 8 or higher is required for Picard
- Check version: `java -version`
- If missing: Install OpenJDK (`sudo apt install openjdk-11-jdk` on Ubuntu, or use conda)

If you run into Java errors:
- Ensure you have a modern Java JDK installed (Java 8+)
- Try: `./gradlew clean && ./gradlew shadowJar` to rebuild Picard
- Check Java version: `java -version`

## Quick Start
To run TEPEAK:
```bash
snakemake --configfile config.yaml --cores 8
```
## Directory Structure

**The directory structure is flexible** and controlled by your config file parameters. You can organize your files however you prefer. The examples below are just suggestions - adjust the paths in your config.yaml to match your preferred structure.

### Minimal Required Structure
```
your_project/
  Snakefile              # (from TEPEAK)
  config.yaml            # Your configuration
  your_reference.zip     # Reference genome 
  your_gtf.zip          # GTF file
  your_samples.txt      # Sample list (for FASTQ mode)
```

### Example Organized Structure
```
TEPEAK/
  Snakefile
  config.yaml
  references/
    ncbi_dataset.zip
    ncbi_dataset_gtf.zip
  samples/
    SraRunTable.txt
    custom_samples.txt
  fastq_data/
    sample1_R1.fastq.gz
    sample1_R2.fastq.gz
  results/              # Output directory (configurable)
```

## Data requirements

TEPEAK can run in two modes:

### 1. NCBI SRA mode (default)

**Required inputs:**
- Zipped genome reference (e.g. `ncbi_dataset.zip`)
- Zipped genome GTF (for gene annotations, e.g. `ncbi_dataset_gtf.zip`)  
- SRA run table. This can be a list of SRA Accession numbers or the SRARunInfo file downloaded from NCBI SRA

**Need help getting the SRA Run Table?** See our [detailed guide on downloading SRA Run Tables from NCBI](https://github.com/ryanlayerlab/TEPEAK/wiki/Getting-the-SRARunInfo-file).


### 2. Custom FASTQ mode

**Required inputs:**
- Zipped genome reference
- Zipped genome GTF
- Directory containing paired-end FASTQs
- Plain text sample list (one sample ID per line)

**FASTQ naming conventions** (any of these work):
- `<sample>_R1.fastq(.gz)` / `<sample>_R2.fastq(.gz)`
- `<sample>_1.fastq(.gz)` / `<sample>_2.fastq(.gz)`
- `<sample>.R1.fastq(.gz)` / `<sample>.R2.fastq(.gz)`
- `<sample>-R1.fastq(.gz)` / `<sample>-R2.fastq(.gz)`

## Configuration

All runs are controlled by a `config.yaml` file. You must specify your input type and configure paths to match your directory structure.

### Core Configuration Parameters

| Parameter | Required | Description | Default |
|-----------|----------|-------------|---------|
| `species` | Yes | Species name for output files | - |
| `data_dir` | Yes | Directory for reference files and working data | - |
| `output_dir` | Yes | Directory for all results | - |
| `zipped_ref_genome_filepath` | Yes | Path to zipped reference genome | - |
| `zipped_gtf_filepath` | Yes | Path to zipped GTF file | - |
| `input_type` | Yes | `"sra"` or `"fastq"` | `"sra"` |
| `threads` | Yes | Number of CPU threads to use | - |
| `low` | Yes | Minimum insertion size (bp) | - |
| `high` | Yes | Maximum insertion size (bp) | - |

### Analysis Control Parameters

| Parameter | Required | Description | Default |
|-----------|----------|-------------|---------|
| `gene` | No | Include gene annotations (`"y"` or `"n"`) | `"n"` |
| `run_smoove` | No | Run smoove deletion analysis (`true` or `false`) | `false` |
| `run_enrichment` | No | Run functional enrichment analysis (`true` or `false`) | `false` |
| `run_phylogeny` | No | Run phylogenetic analysis (`true` or `false`) | `false` |

### Phylogenetic Analysis Parameters

| Parameter | Required | Description | Default |
|-----------|----------|-------------|---------|
| `phylo_min_sequences` | No | Minimum sequences required for phylogeny | `3` |
| `phylo_max_sequences` | No | Maximum sequences to include (for performance) | `50` |

### TEPEAK Advanced Parameters

| Parameter | Required | Description | Default |
|-----------|----------|-------------|---------|
| `tepeak_min_cluster_size` | No | Minimum insertions per size cluster | `10` |
| `tepeak_percentile_threshold` | No | Density threshold percentile (0-100) | `75` |
| `tepeak_window_size` | No | Sliding window size for density calculation | `50` |
| `tepeak_merge_distance` | No | Max bp distance to merge nearby clusters | `100` |
| `tepeak_pid_threshold` | No | Sequence similarity threshold (0-1) | `0.85` |
| `tepeak_max_clusters` | No | Maximum number of size clusters to process | `50` |
| `tepeak_max_sample_seqs` | No | Max sequences to sample per cluster | `200` |
| `tepeak_dfam_evalue` | No | Dfam e-value cutoff | `10.0` |
| `tepeak_dfam_batch_size` | No | Sequences per batch for rate limiting | `5` |
| `tepeak_dfam_delay` | No | Delay between batches (seconds) | `3.0` |
| `tepeak_dfam_organism` | No | Organism for Dfam searches (scientific name) |

**Dfam Organism Examples**: Use scientific names like "Homo sapiens", "Mus musculus", "Equus caballus", "Drosophila melanogaster". Check the [Dfam organism list](https://www.dfam.org/browse) or use NCBI taxonomy IDs.

### Functional Enrichment Parameters

| Parameter | Required | Description | Default |
|-----------|----------|-------------|---------|
| `enrichment_organism` | No | g:Profiler organism code |
| `enrichment_sources` | No | Data sources to query | `['GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC']` |
| `enrichment_max_pvalue` | No | P-value threshold for significance | `0.05` |
| `enrichment_max_terms` | No | Maximum number of terms to report | `50` |

**Supported organisms**: See species supported by g:Profiler.



#### Complete Configuration Examples

#### SRA Input Example
```yaml
species: ecoli
data_dir: data
output_dir: output
zipped_ref_genome_filepath: references/ncbi_dataset.zip
zipped_gtf_filepath: references/ncbi_dataset_gtf.zip
sra_run: 
  filepath: samples/SraRunTable.txt
  number_of_runs: 4
input_type: sra
threads: 10
low: 200
high: 6400
gene: y
run_smoove: false
run_enrichment: true  # Enable functional enrichment
enrichment_organism: hsapiens  
enrichment_max_pvalue: 0.01
# TEPEAK Advanced (optional - these are defaults)
tepeak_min_cluster_size: 10
tepeak_percentile_threshold: 75
tepeak_dfam_organism: "Escherichia coli"  # Scientific name for E. coli
```

#### Custom FASTQ Input Example
```yaml
species: horse
data_dir: data
output_dir: results
zipped_ref_genome_filepath: /path/to/horse_reference.zip
zipped_gtf_filepath: /path/to/horse_gtf.zip
fastq_input:
  sample_list: samples/horse_samples.txt
  fastq_dir: /path/to/horse_fastqs
input_type: fastq
threads: 8
low: 200
high: 6400
gene: y
run_smoove: false
run_enrichment: true
# Custom TEPEAK parameters
tepeak_percentile_threshold: 80
tepeak_max_clusters: 100
tepeak_dfam_evalue: 5.0
tepeak_dfam_organism: "Equus caballus"  # Scientific name for horse
# Custom enrichment parameters
enrichment_organism: Equus caballus  
enrichment_max_pvalue: 0.01
```

## Running TEPEAK

From the TEPEAK directory with the environment activated:
```bash
snakemake --configfile config.yaml --cores 8
```

For dry run (see what will be executed):
```bash
snakemake --configfile config.yaml --cores 8 --dry-run
```

## Output

All results are written to `<output_dir>/<species>/`.

### Core Outputs

- `<species>_global_vcf.txt` — Global VCF of all insertions
- `<species>_insertions_plot.png` — Insertion size histogram
- `<species>_count.txt` — Per-sample insertion counts


- `tepeak_families_annotated.csv` —  Detailed TE family analysis with Dfam annotations
- `tepeak_clusters_plot.png` — Visualization of size-based clusters with detected peaks


### Gene Annotation Outputs (if `gene: y`)

- `peak_<low>-<high>/<species>_<low>-<high>_merged_genes.txt` — Insertions with gene context
- Additional annotation files in `peak_<low>-<high>/` directory

### Functional Enrichment Outputs (if `gene: y` and `run_enrichment: true`)

- `<species>_enrichment_results.csv` — Detailed GO/pathway enrichment statistics
- `<species>_enrichment_summary.txt` — Summary of top enriched functional categories

### Smoove Outputs (if `run_smoove: true`)

- `<species>_smoove_plot.svg` — Deletion size histogram  
- `<species>.smoove.square.vcf.gz` — Multi-sample structural variant calls

### Phylogenetic Outputs (if `run_phylogeny: true`)

- `<species>_te_alignment.fasta` — Multiple sequence alignment of representative TE sequences
- `<species>_te_phylogeny.newick` — Phylogenetic tree showing evolutionary relationships
- `<species>_phylogeny_summary.txt` — Analysis summary with interpretation guide


