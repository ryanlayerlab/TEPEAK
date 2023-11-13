# TEPEAK
A novel method for identifying and characterizing polymorphic transposable elements in  non-model species populations.
## Setup
We recommend running everything inside a `conda` virtual environment using the latest version of `conda` (`23.9.0` at the time of writing). This is because some packages installed by older versions of `conda` may not work properly. 
### Windows
A lot of tools used in this project are only available on UNIX based systems such as macOS and Linux. For proper functionality, we recommend running everything through a Windows Subsytem for Linux (WSL). The instructions to install and set up a WSL on your system are available at https://learn.microsoft.com/en-us/windows/wsl/install. 

### macOS/Linux
Continue following the rest of the setup documentation. 

### Creating INSurVeyor `conda` environment and installing dependecies. 
Run the following command to create the `insurveyor-env` `conda` environment and install INSurVeyor along with other dependencies inside the environment. 
```
mamba env create -f environment.yaml
```
It is highly recommended to run the command above with `mamba` as it is much faster than `conda`. If you do not have `mamba` installed, you can install it by running  the following command inside your `conda` `base` environment. 
```
conda install -c conda-forge mamba
```
After the environment has been created, activate it by running 
```
conda activate insurveyor-env
```

### Other environment requirements
Please view the [wiki](https://github.com/ryanlayer/TEPEAK/wiki/Species-Name-and-SRA-List-Startup) for instructions on how to install and set up the NCBI SDK, Picard, and verifying that you have an appropriate version of Java installed. 

## Running TEPEAK
Create a config file for the species you'd like to analyse.
```
bash src/create_config.sh -s <species> -d <data_dir> -n <number of threads>
```
This creates the file `configs/config_<species>.yaml` and the directory `data_dir/species/`.
<br>`data_dir` is the data directory where all reference and BAM files are going to be produced/stored. 

### Data structure requirements

For SRA and BAM list options, TEPEAK requires a txt/csv file input where each line is a unique sample identifier. This unique sample needs to be attached to a BAM and BAI in the working directory if using the BAM list option. TEPEAK also requires an indexed reference fasta file and an optional GTF file in the same directory. 
```
DATA_DIR/
  SPECIES_DIR/
    species_samples.txt

    SAMPLE1.BAM
    SAMPLE1.BAM.BAI
    SAMPLE2.BAM
    SAMPLE2.BAM.BAI

    species.fa
    species.fa.fai

    species.gtf (optional)
``` 
where `species_samples.csv` looks like 
```
species_samples.txt
  SAMPLE1
  SAMPLE2
```
Note: The samples, reference, and GTF files need to be named after the species. 

### Start Options
#### 1. SRA List

Start with a list of SRA run numbers for a species. 

**Requirements**: 
1. Find reference genome in NCBI's [reference genome database](https://www.ncbi.nlm.nih.gov/datasets/genome/) (See [wiki](https://github.com/ryanlayer/TEPEAK/wiki/Choosing-a-reference-genome) for guide)
2. Setup NCBI SDK and Picard (See [wiki](https://github.com/ryanlayer/TEPEAK/wiki/Species-Name-and-SRA-List-Startup) for instructions)

Once the requirements are setup and a reference and species name is chosen navigate to SRA List Start below.

#### 2. Aligned BAMs

Insurveyor requires coordinate sorted, indexed with MC and MQ tags, BAM files.

Ensure your BAM files meet the data structure requirements and then proceed to Aligned Bams Start below. 

### OPTION 1: SRA Run List Start

- Required data: zipped reference genome downloaded and SRA run list. 
- Required environment setup: NCBI SDK, picard (See [wiki](https://github.com/ryanlayer/TEPEAK/wiki/Species-Name-and-SRA-List-Startup))

The SRA run list must be a txt/csv file with each line being one SRR accession. Name this file  `<species>_samples.txt` and move it to `data_dir/species/`. 

---
If you'd like to generate an SRA run list instead of creating your own, follow the steps below. 
1. Using NCBI's [SRA run selector tool](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/Traces/study/), download the metadata for the desired accession for the chosen species. The downloaded file should be called `SraRunTable.txt`. 

2. Run `python src/get_sra_numbers.py -f <SraRunTable.txt file path and name> -n <max no. of samples> -s <species>`. 
   - This script produces a file called `<species>_samples.txt` inside `data_dir/species` with `<max no. of samples>` samples
   - If you'd like to include all samples, run the script without the `-n` flag. 
---
Now that you have an SRA run list, process the reference genome. 
```
python src/process_reference.py -s <species> -f <zipped genome file path and name> 
```
This produces the `<species>.fa` and `<species>.fa.fai` files inside `data_dir/species/`. 

Ater acquiring the reference files, download the BAM files associated with each sample in `<species>_samples.csv` and align them.
```
bash src/align_species.sh -s <species>
```
This produces `SAMPLE.bam` and `SAMPLE.bam.bai` files inside `data_dir/species/` for each `SAMPLE` in `<species>_samples.txt`. 

Now that you have the aligned BAM files, follow the steps in OPTION 2: Aligned BAMs start. 

### OPTION 2: Aligned BAMs Start
- Required data: reference genome named as `<species>.fa` and Aligned BAM files. 
- Ensure your data matches the data structure requirements.

There are two different options for calling insertions, serial and parallel. If you have a large sample size its highly recommended to use the parallel method instead of the serial. 
- For the serial method, run 
```
python src/call_insertions_serial.py -s <species>
```
- For the parallel method (<ins>please note that this has not been tested yet and may not work</ins>), run 
```
bash src/spawn_parallel.sh -s <species> -n <number of threads> 
```

Insertion call quality depends highly on the sample quality. The following will check the number of insertions per sample.
```
bash src/check_insertions.sh -s <species>
```

This creates a tab deliminated file called `output/count_{species}.txt` where each line is a sample name and the respective number of insertions. Remove unsatisfactory samples from the samplename file before continuing. 

Run the following to generate the global vcf information file. This will also result in the file `output/dfam_annotate.csv` containing the DFAM annotations for any significant peak found in the histogram.
```
bash src/get_global_vcf.sh -s <species>
```

You can get the size frequency histogram for specific ranges by running the following. Omitting the ranges will set the default as 0-10,000bp.

```
python src/build_histogram.py -f <global VCF filepath and name> -l <lower range> -u <upper range>
```

Note: INSurVeyor generates a number of files not directly used in TEPEAK and TEPEAK does not have any garbage collection features. 

Now that you have a range of interest in the histogram extract all sequences with sizes that match.
```
bash src/extract_range.sh -s <species> -l <lower bp range> -u <upper bp range>
```

Annotate loci for genes (requires gtf named as `<species>.gtf`. You can get the gtf file from the NCBI [reference genome database](https://www.ncbi.nlm.nih.gov/datasets/genome/))
```
bash src/annotate_genes.sh -s <species> -l <lower bp range> -u <upper bp range> 
```

Write final output files for a range. Use the ` -g ` flag to include gene annotations (y or n). This will output a final with merged loci information in ` output/species/peak_l-h/` as `_merged.txt`, `_merged_genes.txt` , and `_pop_vcf_.txt`  

```
bash src/write_output.sh -s species -l <lower bp range> -u <upper bp range> -g <(y/n) include genes>
```
---
Note: delete the contents of `prefetch_tmp` when finished
