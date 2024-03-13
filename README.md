# TEPEAK
A novel method for identifying and characterizing polymorphic transposable elements in  non-model species populations.
## Setup
We recommend running everything inside a `conda` virtual environment using the latest version of `conda` (`23.9.0` at the time of writing). This is because some packages installed by older versions of `conda` may not work properly. 

## OS requirements
We recommend running everything on a Linux machine. If you're on Windows, you can run everything through a Windows Subsystem for Linux (WSL). The instructions to install and set up a WSL on a Windows machine are available at https://learn.microsoft.com/en-us/windows/wsl/install.  

### Creating INSurVeyor `conda` environment and installing dependecies. 
Create the environment with `mamba` as it's much faster than `conda`. You can check if you have mamba installed by running
```console
mamba --version
```
This should output the version of `mamba` installed. If `mamba` is not installed, run 
```bash
conda install -n base -c conda-forge mamba
```
to install it into the `base` environment. 

Run the following command to create the `insurveyor-env` environment and install INSurVeyor along with other dependencies inside the environment. 
```console
mamba env create -f environment.yaml
```
After the environment has been created, activate it by running 
```bash
conda activate insurveyor-env
```
### Other environment requirements
Please view the [wiki](https://github.com/ryanlayerlab/TEPEAK/wiki/Installing-NCBI-SDK%2C-Picard%2C-and-Smoove) for instructions on how to install and set up the NCBI SDK, Picard, Smoove, and verifying that you have an appropriate version of Java installed. 

## Data requirements
TEPEAK requires the following files to run properly: 
- Zipped genome reference file 
- Zipped genome gtf file (for gene annotations)
- Sra run list file

The Sra run list file needs to be placed in the TEPEAK directory while the zipped files need to be placed inside the data directory. For example, 
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
where `data` is the data directory, `ncbi_dataset.zip` is the zipped reference genome, `ncbi_dataset_gtf.zip` is the zipped gtf, and `SraRunTable.txt` is the SRA run list. 
## Running TEPEAK
After the required files have been placed in the appropriate locations, edit `config.yaml` to configure it for the species you'd like to analyse. 

For example, to analyse ecoli, the `config.yaml` file would look like 
```yaml
species: ecoli
data_dir: data #name and path of the data directory
output_dir: output #name and path of the output directory
zipped_ref_genome_filepath: data/ncbi_dataset.zip #zipped reference genome file for ecoli
zipped_gtf_filepath: data/ncbi_dataset_gtf.zip  #zipped gtf file for ecoli
sra_run: 
  filepath: SraRunTable.txt #sra run list for ecoli
  number_of_runs: 4 #number of runs to analyse
threads: 10
low: 0 #lower bp range
high: 10000 #upper bp range
gene: y #(y/n) -- includes gene annotations if 'y'
```
Note: If you want to include gene annotations, you must have a zipped gtf file. 

After the config file and the data directory directory structure has been set up, run the TEPEAK pipeline with 
```
snakemake --cores
```
If your config file is named something other than `config.yaml`, run the pipeline with 
```
snakemake --configfile <name and path to config file> --cores 
```
This allows you to define multiple config files for different species and run the pipeline for a specific species. Make sure that the new config file follows the exact same format as above. 

Make sure to always run the pipeline with the `insurveyor-env` activated and from the TEPEAK directory. The pipeline produces two histograms, both of which are located in `<output_dir>/<species>/`. `<species>_

If you're running the pipeline from a remote machine, the pipeline creates the file `<output_dir>/<species>/<species>_plot.svg` which stores the plot. 

---
Note: Delete the contents of `prefetch_tmp` when finished
