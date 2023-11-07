conda config --append channels bioconda
conda config --append channels conda-forge
mamba install -c bioconda bedtools=2.31.0
mamba install -c conda-forge selenium=4.9.0
mamba install -c bioconda -c conda-forge bcftools=1.17
mamba install -c bioconda bwa=0.7.17
mamba install -c bioconda samtools=1.18
mamba install jq=1.6
pip install -r dep/requirements.txt
fastq-dump > /dev/null 2>&1 || export PATH=$PATH:$PWD/$(ls | grep "sratoolkit")/bin