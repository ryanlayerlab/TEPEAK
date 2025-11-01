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

Install Smoove binaries inside your conda bin directory (i.e. `/opt/anaconda3/envs/insurveyor-env/bin`)

```bash
conda activate insurveyor-env
cd /opt/anaconda3/envs/insurveyor-env/bin
wget https://github.com/brentp/smoove/releases/download/v0.2.8/smoove
chmod +x smoove
# Verify
smoove --help
```

## Directory Structure

**The directory structure is completely flexible** and controlled by your config file parameters. You can organize your files however you prefer. The examples below are just suggestions - adjust the paths in your config.yaml to match your preferred structure.

### Minimal Required Structure
```
your_project/
  Snakefile              # (from TEPEAK)
  config.yaml            # Your configuration
  your_reference.zip     # Reference genome (can be anywhere)
  your_gtf.zip          # GTF file (can be anywhere)
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

The key point: **Set `data_dir`, `output_dir`, and file paths in config.yaml to match your preferred organization.**

## Data requirements

TEPEAK can run in two modes:

### 1. NCBI SRA mode (default)

**Required inputs:**
- Zipped genome reference (e.g. `ncbi_dataset.zip`)
- Zipped genome GTF (for gene annotations, e.g. `ncbi_dataset_gtf.zip`)  
- SRA run table (CSV/TSV from NCBI)

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
| `tepeak_dfam_organism` | No | Organism for Dfam searches (scientific name) | Auto-detect |

**Dfam Organism Examples**: Use scientific names like "Homo sapiens", "Mus musculus", "Equus caballus", "Drosophila melanogaster". If auto-detection fails, check the [Dfam organism list](https://www.dfam.org/browse) or use NCBI taxonomy IDs.

### Functional Enrichment Parameters

| Parameter | Required | Description | Default |
|-----------|----------|-------------|---------|
| `enrichment_organism` | No | g:Profiler organism code (auto-detected if not set) | Auto-detect |
| `enrichment_sources` | No | Data sources to query | `['GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC']` |
| `enrichment_max_pvalue` | No | P-value threshold for significance | `0.05` |
| `enrichment_max_terms` | No | Maximum number of terms to report | `50` |

**Supported organisms**: Human (hsapiens), Mouse (mmusculus), Rat (rnorvegicus), Fly (dmelanogaster), Worm (celegans), Yeast (scerevisiae), Zebrafish (drerio), Arabidopsis (athaliana), Rice (osativa), and others supported by g:Profiler.

### Input Type Specific Parameters

**For SRA mode (`input_type: sra`):**
```yaml
sra_run: 
  filepath: path/to/SraRunTable.txt
  number_of_runs: 4
```

**For FASTQ mode (`input_type: fastq`):**
```yaml
fastq_input:
  sample_list: path/to/samples.txt
  fastq_dir: path/to/fastq_directory
```

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
enrichment_organism: hsapiens  # Override if needed
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
enrichment_organism: hsapiens  # Use human annotations for horse (if desired)
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

### TEPEAK Advanced Outputs (New!)

- `tepeak_families_annotated.csv` — **Main result**: Detailed TE family analysis with Dfam annotations
- `tepeak_clusters_plot.png` — Visualization of size-based clusters with detected peaks
- `dfam_annotate.csv` — Legacy Dfam annotation results (for compatibility)

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

## TEPEAK Advanced Analysis

The new TEPEAK advanced analysis provides sophisticated TE family identification through:

1. **Adaptive Size Clustering**: Identifies insertion size peaks using local density analysis
2. **Sequence Similarity Clustering**: Groups insertions by sequence similarity within each size cluster  
3. **Representative Selection**: Chooses consensus sequences for each family
4. **Dfam Annotation**: Annotates representatives against the Dfam database with relaxed thresholds

### Key Features

- **Adaptive thresholds** that adjust to local insertion density
- **Hierarchical clustering** that first groups by size, then by sequence
- **Smart sampling** that oversamples near size peaks
- **Relaxed Dfam matching** with configurable e-value cutoffs
- **Rate limiting** for Dfam API compliance
- **Comprehensive visualization** of detected clusters

### Interpreting Results

The main output `tepeak_families_annotated.csv` contains:

- **Size clusters**: Groups of insertions with similar sizes
- **Sequence families**: Subgroups within each size cluster with similar sequences  
- **Abundance estimates**: Projected total insertions for each family
- **Dfam annotations**: TE class, family, and confidence metrics
- **Representative sequences**: Consensus sequence for each family

## Functional Enrichment Analysis (Optional)

When `gene: y` and `run_enrichment: true` are set, TEPEAK performs functional enrichment analysis on genes intersecting with TE insertions using the g:Profiler web service.

### Features

- **Automatic organism detection** from species name or manual override
- **Multiple data sources**: Gene Ontology (BP, MF, CC), KEGG pathways, Reactome pathways
- **Statistical correction**: FDR-adjusted p-values
- **Comprehensive output**: Detailed results CSV + human-readable summary

### Interpretation

The enrichment analysis identifies biological processes, molecular functions, cellular components, and pathways that are overrepresented among genes with nearby TE insertions. This can reveal:

- **Functional bias**: Whether TEs preferentially insert near genes of specific functions
- **Pathway disruption**: Metabolic or signaling pathways potentially affected by TE activity
- **Regulatory impact**: Enrichment in transcription factors, chromatin modifiers, etc.

### Output Files

- `<species>_enrichment_results.csv` — Detailed enrichment statistics for all significant terms
- `<species>_enrichment_summary.txt` — Human-readable summary of top enriched terms by category

### Requirements

- Requires `gene: y` (gene annotation must be enabled)
- Requires internet connection for g:Profiler API access
- Minimum 5 genes recommended for meaningful results

## Phylogenetic Analysis (Optional)

When `run_phylogeny: true` is set, TEPEAK performs phylogenetic reconstruction of representative TE sequences using MAFFT for alignment and FastTree for tree building.

### Features

- **Automatic sequence selection** from annotated TE families
- **Multiple sequence alignment** with MAFFT
- **Fast phylogenetic reconstruction** with FastTree (approximate maximum-likelihood)
- **Informative sequence labels** including family info and TE annotations
- **Performance optimization** with configurable sequence limits

### Output Files

- `<species>_te_alignment.fasta` — Multiple sequence alignment of representative TEs
- `<species>_te_phylogeny.newick` — Phylogenetic tree in Newick format
- `<species>_phylogeny_summary.txt` — Analysis summary and interpretation guide

### Requirements

- Requires MAFFT and FastTree (included in environment.yaml)
- Minimum 3 representative sequences (configurable)
- Works best with 10-50 sequences for meaningful phylogenies

### Interpretation

The phylogenetic tree reveals:
- **Evolutionary relationships** among TE families within your dataset
- **Sequence divergence** through branch lengths
- **Family clustering** of related transposable elements
- **Novel or divergent TEs** appearing as long branches

Use tree visualization software like FigTree, iTOL, or R packages (ggtree, ape) to explore results.

## Notes

- Always activate the `insurveyor-env` environment before running
- For large references, indexing (`bwa index`) may take hours
- TEPEAK advanced analysis may take significant time due to Dfam API rate limits
- Delete working directories when finished to save space
- If running on Windows, use WSL
- Directory structure is completely flexible - organize files however you prefer and adjust config paths accordingly

## Troubleshooting

### Directory Structure Issues
- **Problem**: "File not found" errors
- **Solution**: Check that all paths in config.yaml are correct relative to where you run snakemake

### TEPEAK Advanced Issues  
- **Problem**: Dfam annotation timeouts or organism errors
- **Solution**: Set `tepeak_dfam_organism` to the correct scientific name (e.g., "Equus caballus" for horse)
- **Problem**: No clusters detected
- **Solution**: Lower `tepeak_percentile_threshold` or `tepeak_min_cluster_size`

### Performance Issues
- **Problem**: Analysis is very slow
- **Solution**: Reduce `tepeak_max_clusters`, increase `threads`, or skip smoove with `run_smoove: false`

### Dfam API Issues
- **Problem**: All Dfam queries fail
- **Solution**: Check internet connection and verify organism name at https://www.dfam.org/browse
