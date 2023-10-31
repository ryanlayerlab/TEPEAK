conda config --append channels bioconda
conda config --append channels conda-forge
conda install -c bioconda bedtools
conda install -c conda-forge selenium
conda install -c bioconda -c conda-forge bcftools
conda install -c bioconda bwa
conda install -c bioconda samtools
conda install jq
pip install -r dep/requirements.txt