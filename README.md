# TEPEAK
novel method for identifying and characterizing polymorphic transposable elements in  non-model species populations
---
### Setup

#### Installing INSurVeyor and other requirements

Create the insurveyor conda virtual environment and install insurveyor inside the environment. 
```
conda create -n insurveyor-env -c bioconda -c conda-forge insurveyor
```
Activate the insurveyor environment. 
```
conda activate insurveyor-env
```
Conda is the preferred method for installing insurveyor. Singularity and source options are also available at https://github.com/kensung-lab/INSurVeyor

Finish setup by running
```
bash dep/setup.sh
```
### Start Options

#### 1. Species Name

Simply supply a species name and the location of a desired reference genome and TEPEAK will do the rest.


**Requirements**: Find reference genome in NCBI database [Reference genome database](https://www.ncbi.nlm.nih.gov/datasets/genome/) (See [wiki](https://github.com/mrburke00/TEPEAK/wiki/Choosing-a-reference-genome) for guide)

Setup webscraper, NCBI SDK, and picard (See [wiki](https://github.com/mrburke00/TEPEAK/wiki/Species-Name-and-SRA-List-Startup) for instructions)

Once the requirements are setup and a reference and species name is chosen navigate to Species Name Start below

#### 2. SRA List

Start with a list of SRA number for the same species

**Requirements**: Find reference genome in NCBI database [Reference genome database](https://www.ncbi.nlm.nih.gov/datasets/genome/) (See [wiki](https://github.com/mrburke00/TEPEAK/wiki/Choosing-a-reference-genome) for guide)

Setup NCBI SDK, and picard (See [wiki](https://github.com/mrburke00/TEPEAK/wiki/Species-Name-and-SRA-List-Startup) for instructions)

Once the requirements are setup and a reference and species name is chosen navigate to SRA List Start below


#### 3. Aligned BAMs

Insurveyor requires coordinate sorted, indexed with MC and MQ tags, BAM files.

Ensure your BAM meet the data requirements below then proceed to Aligned Bams Start below

---

### Data Requirements

For SRA and BAM list options, TEPEAK requires a txt file input where each line is a unique sample identifier. This unique sample needs to be attached to a BAM and BAI in the working 
directory if using the BAM list option. TEPEAK also requires an indexed reference fasta file and an optional GTF file in the same directory. 
```
horse_samples.csv
  SAMPLE1
  SAMPLE2


DATA_DIR/
  SAMPLE1.BAM
  SAMPLE1.BAM.BAI
  SAMPLE2.BAM
  SAMPLE2.BAM.BAI

  horse.fa  (ref must be named <species>.fa)
  horse.fa.fai

  horse.gtf (optional)
``` 
Note: The reference and GTF file need to be named after the species. 

---
## Running TEPEAK
---
## OPTION 1: Species Name Start

Required data: zipped reference genome downloaded and species name 

Required environment setup: webscraper, NCBI SDK, picard (See [wiki](https://github.com/mrburke00/TEPEAK/wiki/Species-Name-and-SRA-List-Startup))

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

Output will be a tab deliminated file ```count_{species}.txt``` where each line is sample name and respective number of insertions. Remove unsatisfactory samples from samplename file before continuing. 

Run the following to generate the global vcf information file and overall size-frequency histogram. This will also result in the file ```output/dfam_annotate.csv``` containing the DFAM annotations for any significant peak found in the histogram.

7. ```bash get_global_vcf.sh -s <species>```

You can get the histogram for specific ranges by running the following. Omitting the ranges will set the default as 0-10,000bp.

```python build_histogram.py -f <global VCF filename> -s <species> -l <lower range> -u <upper range>```

Now that you have a range of interest in the histogram extract all sequences with sizes that match

8. ``` bash extract_range.sh -s <species> -l <lower bp range> -u <upper bp range> ```

Annotate loci for genes (requires gtf named as ```<species>.gtf```)

9. ```bash annotate_genes -s <species> -l <lower bp range> -u <upper bp range> ``` 

Write final output files for a range. Use the ``` -g ``` flag to include gene annotations (y or n). This will output a final with merged loci information in ``` output/species/peak_l-h/``` as ```_merged.txt```, ```_merged_genes.txt``` , and ```_pop_vcf_.txt```  

10. ```bash write_output -s species -l <lower bp range> -u <upper bp range> -g <(y/n) include genes> ``` 

Note: delete the contents of ```prefetch_tmp``` when finished

---

## OPTION 2: SRA List Start

Required data: zipped reference genome downloaded and SRA list

Required environment setup: NCBI SDK, picard (See [wiki](https://github.com/mrburke00/TEPEAK/wiki/Species-Name-and-SRA-List-Startup))

Prepare reference genome 

SRA list must be txt file with each line being one SRA accession. Name this file  ```<species>_samples.txt``` and move it to TEPEAK directory

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

Run the following to generate the global vcf information file and overall size-frequency histogram. This will also result in the file ```output/dfam_annotate.csv``` containing the DFAM annotations for any significant peak found in the histogram.

7. ```bash get_global_vcf.sh -s <species>```

You can get the histogram for specific ranges by running the following. Omitting the ranges will set the default as 0-10,000bp.

```python build_histogram.py -f <global VCF filename> -s <species> -l <lower range> -u <upper range>```
Now that you have a range of interest in the histogram extract all sequences with sizes that match

8. ``` bash extract_range.sh -s <species> -l <lower bp range> -u <upper bp range> ```

Annotate loci for genes (requires gtf named as ```<species>.gtf```)

9. ```bash annotate_genes -s <species> -l <lower bp range> -u <upper bp range> ``` 

Write final output files for a range. Use the ``` -g ``` flag to include gene annotations (y or n). This will output a final with merged loci information in ``` output/species/peak_l-h/``` as ```_merged.txt```, ```_merged_genes.txt``` , and ```_pop_vcf_.txt```  

10. ```bash write_output -s species -l <lower bp range> -u <upper bp range> -g <(y/n) include genes> ``` 

Note: delete the contents of ```prefetch_tmp``` when finished


---

## OPTION 3: Aligned BAMs Start

Required data: reference genome named as ```<species>.fa``` and Aligned BAMs list

Ensure your data matches the data structure requirements. Name your list of BAM ids as  ```<species>_samples.txt``` and move it to TEPEAK directory

As of now the BAMs need to be inside the TEPEAK directory, stored in ```<data directory>```

There are two different options for calling insertions, serial and parallel. If you have a large sample size its highly recommended to use the parallel method. See Parallel Insertion Calling section below

##### Serial Run

Begin by creating a config file

1. ``` bash bam_start_config.sh -s <species> -d <data_dir> -n <number of threads> ```

2. ```bash call_insertions_serial.sh -s <species name> ```

Insertion call quality depends highly on sample quality. The following will check the number of insertions per sample

3. ```bash checkInsertions.sh -s <species>```

Output will be a tab deliminated file ```count_{species}.txt``` where each line is a sample name and the respective number of insertions. Remove unsatisfactory samples from the samplename file before continuing. 

Run the following to generate the global vcf information file and overall size-frequency histogram. This will also result in the file ```output/dfam_annotate.csv``` containing the DFAM annotations for any significant peak found in the histogram.

4. ```bash get_global_vcf.sh -s <species>```

You can get the histogram for specific ranges by running the following. Omitting the ranges will set the default as 0-10,000bp.

5. ```python build_histogram.py -f <global VCF filename> -s <species> -l <lower range> -u <upper range>```

Note: INSurVeyor generates a number of files not directly used in TEPEAK. TEPEAK also does not have any garbage collection features. 

Now that you have a range of interest in the histogram extract all sequences with sizes that match

8. ``` bash extract_range.sh -s <species> -l <lower bp range> -u <upper bp range> ```

Annotate loci for genes (requires gtf named as ```<species>.gtf```)

9. ```bash annotate_genes -s <species> -l <lower bp range> -u <upper bp range> ``` 

Write final output files for a range. Use the ``` -g ``` flag to include gene annotations (y or n). This will output a final with merged loci information in ``` output/species/peak_l-h/``` as ```_merged.txt```, ```_merged_genes.txt``` , and ```_pop_vcf_.txt```  

10. ```bash write_output -s species -l <lower bp range> -u <upper bp range> -g <(y/n) include genes> ``` 


---
### Parallel Insertion Calling

Requires xargs
Determine how many separate jobs you want to start as -p flag (this parameter will divide number of lines in your input sample number file), also easily 
extendible to sbatch script. Requires aligned BAMs

```bash spawn_parallel.sh -f <sample filename> -d <bam data directory> -n <number of threads per process> -s <species name> -p <number of jobs>```

