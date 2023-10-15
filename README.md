# TEPEAK
novel method for identifying and characterizing polymorphic transposable elements in  non-model species populations
---
TODO: bcftools, pandas, bwa
 conda install -c conda-forge selenium
 picard

### Setup

##### Installing INSurVeyor 
Conda is preferred method 

```
conda activate insurveyor-env
conda config --append channels bioconda
conda install insurveyor
``` 
Singularity and source options also available https://github.com/kensung-lab/INSurVeyor

### Start Options

1. Species Name

Simply supply a species name and the location of a desired reference genome and TEPEAK will do the rest. 

Reference genome database https://www.ncbi.nlm.nih.gov/datasets/genome/ (please see Wiki for details)

See Wiki for details on setting up the requirements for this option.

Once the requirements are setup and a reference and species name is chosen navigate to Species Name Start below

2. SRA List

Start with a list of SRA number for the same species

Find the reference genome Reference genome database https://www.ncbi.nlm.nih.gov/datasets/genome/ (please see Wiki for details

See Wiki for details on setting up the requirements for this option

Once the requirements are setup and a reference and species name is chosen navigate to SRA List Start below


3. Aligned BAMs

Insurveyor requires coordinate sorted, indexed with MC and MQ tags, BAM files.

Ensure your BAM meet the data requirements below then proceed to Aligned Bams Start below

---

### Data Requirements

For SRA and BAM list options, TEPEAK requires a txt file input where each line is a unique sample identifier. This unique sample needs to be attached to a BAM and BAI in the working 
directory if using the BAM list option. TEPEAK also requires an indexed reference fasta file and optional GTF file in the same directory. 
```
horse_samples.csv
  SAMPLE1
  SAMPLE2


DATA_DIR/
  SAMPLE1.BAM
  SAMPLE1.BAM.BAI
  SAMPLE2.BAM
  SAMPLE2.BAM.BAI

  horse.fa
  horse.fa.fai

  horse.gtf (optional)
``` 
Note: The reference and GTF file need to be named after the species. 

---
## Running TEPEAK
---
## OPTION 1: Species Name Start

Required data: zipped reference genome downloaded and species name 

Required environment setup: NCBI Scraper, NCBI SDK, picard (SEE WIKI) 

Begin by creating a config file

1. ``` bash species_start_config.sh -s <species> -d <data_dir> -n <number of threads> ```
 
2. ``` python ncbi_scrape.py -s <species name> ```

This will download a file ```SraRunInfo.csv``` to your downloads directory. Once finished downloading either move to TEPEAK directory or copy its path as input to the next script.

2. ``` python get_sra_numbers.py -f <SraRunInfo.csv file path and name> -n <max no. of samples> -s <species name> ```

Prepare reference genome 

3. ``` bash process_reference.sh -s <species name> -f <zipped genome file path and namemvc > ```

Download SRA data and align to reference 

4.  ``` bash align_species.sh -s <species>```

Call insertions (Note if you would like to run parallel jobs see Parallel Run section below)

5. ```bash call_insertions_serial.sh -s <species name> ```

Insertion call quality depends highly on sample quality. The following will check the number of insertions per samples

6. ```bash checkInsertions.sh -s <species>```

Output will be a tab deliminated file ```count_{species}.txt``` where each line is sample name and respective number of insertions. Remove unsatisfactory samples from samplename file before continuiing. 

Run the following to generate the global vcf information file and overall size-frequency histogram.

7. ```bash getGlobalVCF.sh -s <species>```

You can get the histogram for specific ranges by running the following. Omitting the ranges will set the default as 0-10,000bp.

```python buildHistogram.py -f <sample_filename> -s <species> -l <lower range> -u <upper range>```

Note: delete the contents of ```prefetch_tmp``` when finished

---

## OPTION 2: SRA List Start

If you do not already have a reference genome see the wiki about selecting one from the NCBI database 

Prepare reference genome 

SRA list must be txt file with each line being one SRA accessions. Name this file  ```<species>_samples.txt``` and move it to TEPEAK directory

Begin by creating a config file

1. ``` bash sra_start_config.sh -s <species> -d <data_dir> -n <number of threads> ```

Process reference genome

2. ``` bash process_reference.sh -s <species name> -f <zipped genome file path and namemvc > ```

Download and align 

3. ``` bash align_species.sh -s <species>```

Call insertions (Note if you would like to run parallel jobs see Parallel Run section below)

5. ```bash call_insertions_serial.sh -s <species name> ```

Insertion call quality depends highly on sample quality. The following will check the number of insertions per samples

6. ```bash checkInsertions.sh -s <species>```

Output will be a tab deliminated file ```count_{species}.txt``` where each line is sample name and respective number of insertions. Remove unsatisfactory samples from samplename file before continuiing. 

Run the following to generate the global vcf information file and overall size-frequency histogram.

7. ```bash getGlobalVCF.sh -s <species>```

You can get the histogram for specific ranges by running the following. Omitting the ranges will set the default as 0-10,000bp.

```python buildHistogram.py -f <sample_filename> -s <species> -l <lower range> -u <upper range>```

---

## OPTION 3: Aligned BAMs Start

Ensure your data matches the data structure requirements. Name your list of BAM ids as  ```<species>_samples.txt``` and move it to TEPEAK directory

As of now the BAMs need to be inside the TEPEAK directory, stored in ```<data directory>```

There are two different options for calling insertions, serial and parallel. If you have a large sample size its highly reccomenmded the parallel method 
is used. See Parallel Insertion Calling section below

##### Serial Run

Begin by creating a config file

1. ``` bash bam_start_config.sh -s <species> -d <data_dir> -n <number of threads> ```

2. ```bash call_insertions_serial.sh -s <species name> ```

Insertion call quality depends highly on sample quality. The following will check the number of insertions per samples

3. ```bash checkInsertions.sh -s <species>```

Output will be a tab deliminated file ```count_{species}.txt``` where each line is sample name and respective number of insertions. Remove unsatisfactory samples from samplename file before continuiing. 

Run the following to generate the global vcf information file and overall size-frequency histogram.

4. ```bash getGlobalVCF.sh -s <species>```

You can get the histogram for specific ranges by running the following. Omitting the ranges will set the default as 0-10,000bp.

5. ```python buildHistogram.py -f <sample_filename> -s <species> -l <lower range> -u <upper range>```

Note: INSurVeyor generates a number of files not directly used in TEPEAK. TEPEAK also does not have any garbage collection feature. 

### Analysis

#### Checking quality of insertions calls
Insertion call quality depends highly on sample quality. Run the following to check the number of insertions per samples

```bash checkInsertions.sh -f <sample_filename> -s <species>```

Output will be a tab deliminated file ```count_{species}.txt``` where each line is sample name and respective number of insertions. Remove unsatisfactory samples from samplename file before continuiing. 

#### Generating Size-Frequency Histogram and getting preliminary DFAM Annotations
Run the following to generate the global vcf information file and overall size-frequency histogram.

```bash getGlobalVCF.sh -f <sample_filename> -s <species>```

You can get the histogram for specific ranges by running the following. Omitting the ranges will set the default as 0-10,000bp.

```python buildHistogram.py -f <sample_filename> -s <species> -l <lower range> -u <upper range>```

---

### Parallel Insertion Calling

TODO CONFIG FILE

Requires xargs
Determine how many separate jobs you want to start as -p flag (this parameter will divide number of lines in your input sample number file), also easily 
extendible to sbatch script

```bash spawn_parallel.sh -f <sample filename> -d <data directory> -n <number of threads per process> -s <species name> -p <number of jobs>```

