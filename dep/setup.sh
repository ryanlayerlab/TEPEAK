conda config --append channels bioconda
conda config --append channels conda-forge
conda install -c bioconda bedtools=2.31.0
conda install -c conda-forge selenium=4.9.0
conda install -c bioconda -c conda-forge bcftools=1.17
conda install -c bioconda bwa=0.7.17
conda install -c bioconda samtools=1.18
pip install -r dep/requirements.txt
fastq-dump > /dev/null 2>&1 || export PATH=$PATH:$PWD/$(ls | grep "sratoolkit")/bin