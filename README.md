# TEPEAK
novel method for identifying and characterizing polymorphic transposable elements in  non-model species populations

TODO: bcftools

### Setup

##### Installing INSurVeyor 
Conda is preferred method 

```
conda activate insurveyor-env
conda config --append channels bioconda
conda install insurveyor
``` 
Singularity and source options also available https://github.com/kensung-lab/INSurVeyor


### Data Requirements
Singularity requires coordinate sorted, indexed with MC and MQ tags, BAM files. 

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


### Running TEPEAK


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



####

