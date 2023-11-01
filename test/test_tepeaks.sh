test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

# check if $PATH points to sratoolkit/bin, if not point it
fastq-dump > /dev/null 2>&1 || export PATH=$PATH:$PWD/$(ls | grep "sratoolkit")/bin

# fresh start by removing all the directories being produced
# rm -rf configs/config_ecoli.yaml
# rm -rf data/ecoli output count_ecoli.txt
# echo "Files have been removed. Starting afresh."; echo 

run test_species_start_config bash src/species_start_config.sh -s ecoli -d data -n 1
#test that data_dir is made 
assert_equal "true" $(test -d data && echo true)
#test that a data_dir/species_dir is made
assert_equal "true" $(test -d data/ecoli && echo true)
#test that a config file is made 
assert_equal configs/config_ecoli.yaml $(ls configs/config_ecoli.yaml)
#test that "species: $species" is in $config_file
assert_equal "ecoli" $(cat configs/config_ecoli.yaml | grep species | cut -d " " -f 2)
#test that "data_directory: $data_dir" is in $config_file
assert_equal "data" $(cat configs/config_ecoli.yaml | grep data_directory | cut -d " " -f 2)
#test that "threads: $threads" is in $config_file
assert_equal "1" $(cat configs/config_ecoli.yaml | grep thread | cut  -d " " -f 2)
# ===============

run test_get_sra_numbers python src/get_sra_numbers.py -f SraRunTable.txt -n 4 -s ecoli
assert_equal "ERR10355883" $(cat data/ecoli/ecoli_samples.txt | grep ERR10355883)
assert_equal "ERR10355891" $(cat data/ecoli/ecoli_samples.txt | grep ERR10355891)
assert_equal "ERR10355911" $(cat data/ecoli/ecoli_samples.txt | grep ERR10355911)
assert_equal "ERR10355944" $(cat data/ecoli/ecoli_samples.txt | grep ERR10355944)
# ===============

# run test_process_reference bash src/process_reference.sh -s ecoli -f data/ncbi_dataset.zip 
run test_process_reference python src/process_reference.py -s ecoli -f data/ncbi_dataset.zip
# testing the creation of files inside data_dir/species_dir
assert_equal "data/ecoli/ecoli.fa" $(ls data/ecoli/ecoli.fa)
assert_equal "data/ecoli/ecoli.amb" $(ls data/ecoli/ecoli.amb)
assert_equal "data/ecoli/ecoli.ann" $(ls data/ecoli/ecoli.ann)
assert_equal "data/ecoli/ecoli.bwt" $(ls data/ecoli/ecoli.bwt)
assert_equal "data/ecoli/ecoli.fa.fai" $(ls data/ecoli/ecoli.fa.fai)
assert_equal "data/ecoli/ecoli.pac" $(ls data/ecoli/ecoli.pac)
assert_equal "data/ecoli/ecoli.sa" $(ls data/ecoli/ecoli.sa)

# testing the contents of files created inside data_dir/species_dir
assert_equal "Escherichia coli" $(cat data/ecoli/ecoli.fa | head -n 1 | cut -d " " -f 2,3)
assert_equal "4641652 1 0" $(cat data/ecoli/ecoli.amb)
assert_equal "Escherichia coli" $(cat data/ecoli//ecoli.ann | grep NC | cut -d " " -f 3,4)
assert_equal "NC_000913.3	4641652	72	80	81" $(car data/ecoli/ecoli.fa.fai)
# ===============

run test_align_species bash src/align_species.sh -s ecoli # this script takes a while. 
# testing the creation of .bam and .bam.bai files
assert_equal "data/ecoli/ERR10355883.bam" $(ls data/ecoli/ERR10355883.bam)
assert_equal "data/ecoli/ERR10355883.bam.bai" $(ls data/ecoli/ERR10355883.bam.bai)

assert_equal "data/ecoli/ERR10355891.bam" $(ls data/ecoli/ERR10355891.bam)
assert_equal "data/ecoli/ERR10355891.bam.bai" $(ls data/ecoli/ERR10355891.bam.bai)

assert_equal "data/ecoli/ERR10355911.bam" $(ls data/ecoli/ERR10355911.bam)
assert_equal "data/ecoli/ERR10355911.bam.bai" $(ls data/ecoli/ERR10355911.bam.bai)

assert_equal "data/ecoli/ERR10355944.bam" $(ls data/ecoli/ERR10355944.bam)
assert_equal "data/ecoli/ERR10355944.bam.bai" $(ls data/ecoli/ERR10355944.bam.bai)
# ===============

## for some reason the .py script produces the desired output while the .sh script does not. 
run test_call_insertions_serial python src/call_insertions_serial.py -s ecoli
# run test_call_insertions_serial bash src/call_insertions_serial.sh -s ecoli
# testing that output/ is made
assert_equal "true" $(test -d output && echo true)
# testing that output/species is made
assert_equal "true" $(test -d output/ecoli && echo true)
# testing creation of other folders 
assert_equal "true" $(test -d output/ecoli/ERR10355891 && echo true)
assert_equal "true" $(test -d output/ecoli/ERR10355944 && echo true)
assert_equal "true" $(test -d output/ecoli/ERR10355883 && echo true)
assert_equal "true" $(test -d output/ecoli/ERR10355911 && echo true)

# test the existence of out.pass.vcf.gz files inside the above directories -- file is used later on. 
assert_equal "output/ecoli/ERR10355891/out.pass.vcf.gz" $(ls output/ecoli/ERR10355891/out.pass.vcf.gz)
assert_equal "output/ecoli/ERR10355944/out.pass.vcf.gz" $(ls output/ecoli/ERR10355944/out.pass.vcf.gz)
assert_equal "output/ecoli/ERR10355883/out.pass.vcf.gz" $(ls output/ecoli/ERR10355883/out.pass.vcf.gz)
assert_equal "output/ecoli/ERR10355911/out.pass.vcf.gz" $(ls output/ecoli/ERR10355911/out.pass.vcf.gz)
# ===============

run test_check_insertions bash src/check_insertions.sh -s ecoli 
# test to see if count_{species}.txt has been created
assert_equal "count_ecoli.txt" $(ls count_ecoli.txt)
# test the contents of counts_{species}.txt
assert_equal "2" $(cat count_ecoli.txt | grep ERR10355883 | cut -f 2)
assert_equal "2" $(cat count_ecoli.txt | grep ERR10355891 | cut -f 2)
assert_equal "2" $(cat count_ecoli.txt | grep ERR10355911 | cut -f 2)
assert_equal "2" $(cat count_ecoli.txt | grep ERR10355944 | cut -f 2)
# ===============

run test_get_global_vcf bash src/get_global_vcf.sh -s ecoli 
# test the creation of output/{species}_global_vcf.txt
assert_equal "output/ecoli_global_vcf.txt" $(ls output/ecoli_global_vcf.txt)
# test the creation of output/dfam_annotate.csv
assert_equal "output/dfam_annotate.csv" $(ls output/dfam_annotate.csv) # this file only contains the column headers for this particular dataset
# ===============

run test_extract_range bash src/extract_range.sh -s ecoli -l 0 -u 10000
# test the creation of output/peak_low-high folder
assert_equal "true" $(test -d output/peak_0-10000 && echo true)
# test the creation of output/peak_low-high/species_low-high_pop_vcf.txt file
assert_equal "output/ecoli/peak_0-10000/ecoli_0-10000_pop_vcf.txt" $(ls output/ecoli/peak_0-10000/ecoli_0-10000_pop_vcf.txt)