# TEPEAK
A novel method for identifying and characterizing polymorphic transposable elements in  non-model species populations.
## Setup
We recommend running everything inside a `conda` virtual environment using the latest version of `conda` (`23.9.0` at the time of writing). This is because some packages installed by older versions of `conda` may not work properly. 
### Windows
Most of the tools used in this project are only available on UNIX based systems such as macOS and Linux. For proper functionality, we recommend running everything through a Windows Subsytem for Linux (WSL). The instructions to install and set up a WSL on your system are available at https://learn.microsoft.com/en-us/windows/wsl/install. 

### macOS/Linux
Continue following the rest of the setup documentation. 

### Creating INSurVeyor `conda` environment and installing dependecies. 
Create the environment with `mamba` as it's much faster than `conda`. You can check if you have mamba installed by running
```
mamba --version
```
This should output the version of `mamba` installed. If `mamba` is not installed, run 
```
conda install -n base -c conda-forge mamba
```
to install it into the `base` environment. 

Run the following command to create the `insurveyor-env` environment and install INSurVeyor along with other dependencies inside the environment. 
```
mamba env create -f environment.yaml
```
After the environment has been created, activate it by running 
```
conda activate insurveyor-env
```

### Other environment requirements
Please view the [wiki](https://github.com/ryanlayer/TEPEAK/wiki/Species-Name-and-SRA-List-Startup) for instructions on how to install and set up the NCBI SDK, Picard, and verifying that you have an appropriate version of Java installed. 

## Running TEPEAK
Create a `config.yaml` file in the TEPEAK directory for the species you'd like to analyse with the format
```
species: # name of species
data_dir: # name of data directory
output_dir: # name of output directory
zipped_ref_genome_filepath: # path to zipped genome reference file
zipped_gtf_filepath: # path to zipped gtf reference file (optional)
sra_run: 
  filepath: # path to SRA run file 
  number_of_runs: # number of SRA runs
threads: # number of threads to download .bam files
low: 0 
high: 10000 
gene:  # ('y' or 'n')
```
For example, to analyse ecoli, a sample `config.yaml` file would look like 
```
species: ecoli
data_dir: data
output_dir: output
zipped_ref_genome_filepath: data/ncbi_dataset.zip 
zipped_gtf_filepath: data/ncbi_dataset_gtf.zip 
sra_run: 
  filepath: SraRunTable.txt 
  number_of_runs: 4
threads: 1
low: 0 
high: 10000 
gene: 'y'
```
and the TEPEAK directory would look like 
```
TEPEAK/
  snakefile
  config.yaml
  SraRunTable.txt
  data/
    ncbi_dataset.zip
    ncbi_dataset_gtf.zip
  ...
```
where `ncbi_dataset.zip` is the zipped reference genome for ecoli, `ncbi_dataset_gtf.zip` is the zipped gtf for ecoli, and `SraRunTable.txt` is the SRA run list for ecoli. 

After the config file and the data directory directory structure has been set up, run the TEPEAK pipeline with 
```
snakemake --cores
```
If your config file is named something other than `config.yaml`, run the pipeline with 
```
snakemake --configfile <path to config file> --cores 
```
This allows you to define multiple config files for different species and run the pipeline for a specific species. 

---
Note: delete the contents of `prefetch_tmp` when finished