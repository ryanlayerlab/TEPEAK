# TEPEAK
novel method for identifying and characterizing polymorphic transposable elements in  non-model species populations

TODO: bcftools, pandas, bwa
 conda install -c conda-forge selenium

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

2. SRA Numbers

TODO

3. Aligned BAMs

Insurveyor requires coordinate sorted, indexed with MC and MQ tags, BAM files.

Ensure your BAM meet the data requirements below then proceed to Aligned Bams Start below

### Data Requirements


TEPEAK requires a txt file input where each line is a unique sample identifier. This unique sample needs to be attached to a BAM and BAI in the working 
directory. TEPEAK also requires an indexed reference fasta file and optional GTF file in the same directory. 
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


## Running TEPEAK

### Species Name Start

``` python ncbi_scrape.py -s <species name> ```

This will download a file ```SraRunInfo.csv``` to your downloads directory. Once finished downloading either move to TEPEAK directory or copy its path as input to the next script

``` python get_sra_numbers.py -f <SRA file path> -n <max no. of samples> -o <output file name> ```

To prepare and index the reference genome run 

``` bash process_reference.sh -s <species name> -f <zipped genome file path> ```

You can now proceed to the SRA List Start

### SRA List Start

If you do not already have a reference genome see the wiki about selecting one from the NCBI database and run 


``` bash process_reference.sh -s <species name> -f <zipped genome file path> ```


### Aligned BAMs Start

#### Calling Insertions
There are two different options for calling insertions, serial and parallel. If you have a large sample size its highly reccomenmded the parallel method 
is used.

##### Serial Run

```bash call_insertions_serial.sh -f <sample filename> -d <data directory> -n <number of threads> -s <species name> ```

##### Parallel Run
Requires xargs
Determine how many separate jobs you want to start as -p flag (this parameter will divide number of lines in your input sample number file), also easily 
extendible to sbatch script

```bash spawn_parallel.sh -f <sample filename> -d <data directory> -n <number of threads per process> -s <species name> -p <number of jobs>```

Note: INSurVeyor generates a number of files not directly used in TEPEAK. TEPEAK also does not have any garbage collection feature. 

#### Checking quality of insertions calls
Insertion call quality depends highly on sample quality. Run the following to check the number of insertions per samples

```bash checkInsertions.sh -f <sample_filename> -s <species>```

Output will be a tab deliminated file ```count_{species}.txt``` where each line is sample name and respective number of insertions. Remove unsatisfactory samples from samplename file before continuiing. 

#### Generating Size-Frequency Histogram and getting preliminary DFAM Annotations
Run the following to generate the global vcf information file and overall size-frequency histogram.

```bash getGlobalVCF.sh -f <sample_filename> -s <species>```

You can get the histogram for specific ranges by running the following. Omitting the ranges will set the default as 0-10,000bp.

```python buildHistogram.py -f <sample_filename> -s <species> -l <lower range> -u <upper range>```

####

