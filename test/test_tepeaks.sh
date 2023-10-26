test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

rm -rf configs/config_ecoli.yaml

run test_species_start_config bash src/species_start_config.sh -s ecoli -d data -n 1
#test that data_dir is made 
data_exists=false
if [ -d "data" ]; then 
    data_exists=true
fi
assert_equal "true" "$data_exists"
# assert_equal true $(-d "data")

#test that a data_dir/species_dir is made
species_exists=false
if [ -d "data/ecoli" ]; then 
    species_exists=true
fi
assert_equal "true" "$species_exists"

#test that a config file is made 
assert_equal configs/config_ecoli.yaml $(ls configs/config_ecoli.yaml)
#test that "species: $species" is in  $config_file
assert_equal "ecoli" $(cat configs/config_ecoli.yaml | grep species | cut -d " " -f 2)
#test that "data_directory: $data_dir" is in $config_file
assert_equal "data" $(cat configs/config_ecoli.yaml | grep data_directory | cut -d " " -f 2)
#test that "threads: $threads" is in $config_file
assert_equal "1" $(cat configs/config_ecoli.yaml | grep thread | cut  -d " " -f 2)


run test_get_sra_numbers python src/get_sra_numbers.py -f SraRunTable.txt -n 4 -s ecoli
assert_equal "ERR10355883" $(cat data/ecoli/ecoli_samples.txt | grep ERR10355883)
assert_equal "ERR10355891" $(cat data/ecoli/ecoli_samples.txt | grep ERR10355891)
assert_equal "ERR10355911" $(cat data/ecoli/ecoli_samples.txt | grep ERR10355911)
assert_equal "ERR10355944" $(cat data/ecoli/ecoli_samples.txt | grep ERR10355944)

# run test_process_reference.sh bash src/process_reference.sh -s ecoli -f data/ncbi_dataset.zip 
run test_process_reference.sh python src/process_reference.py -s ecoli -f data/ncbi_dataset.zip
# testing the creation of files innside data_dir/species_dir
assert_equal data/ecoli/ecoli.fa $(ls data/ecoli/ecoli.fa)
assert_equal data/ecoli/ecoli.amb $(ls data/ecoli/ecoli.amb)
assert_equal data/ecoli/ecoli.ann $(ls data/ecoli/ecoli.ann)
assert_equal data/ecoli/ecoli.bwt $(ls data/ecoli/ecoli.bwt)
assert_equal data/ecoli/ecoli.fa.fai $(ls data/ecoli/ecoli.fa.fai)
assert_equal data/ecoli/ecoli.pac $(ls data/ecoli/ecoli.pac)
assert_equal data/ecoli/ecoli.sa $(ls data/ecoli/ecoli.sa)

# testing the contents of files created inside data_dir/species_dir
