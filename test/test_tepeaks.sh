test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

rm -rf config_ecoli.yaml

run test_species_start_config bash species_start_config.sh -s ecoli -d data -n 1
#test that data_dir is made 
data_exists=false
if [ -d "data" ]; then 
    data_exists=true
fi
assert_equal "true" "$data_exists"
# assert_equal true $(-d "data")

#test that a specices sub dir is made
species_exists=false
if [ -d "data/ecoli" ]; then 
    species_exists=true
fi
assert_equal "true" "$species_exists"

#test that a config file is made 
assert_equal config_ecoli.yaml $(ls config_ecoli.yaml)
#test that "species: $species" is in  $config_file
assert_equal "ecoli" $(cat config_ecoli.yaml | grep species | cut -d " " -f 2)
#test that "data_directory: $data_dir" is in $config_file
assert_equal "data" $(cat config_ecoli.yaml | grep data_directory | cut -d " " -f 2)
#test that "threads: $threads" is in $config_file
assert_equal "1" $(cat config_ecoli.yaml | grep thread | cut  -d " " -f 2)


run test_get_sra_numbers python get_sra_numbers.py -f SraRunTable.txt -n 4 -s ecoli
assert_equal "ERR10355883" $(cat data/ecoli/ecoli_samples.txt | grep ERR10355883)
assert_equal "ERR10355891" $(cat data/ecoli/ecoli_samples.txt | grep ERR10355891)
assert_equal "ERR10355911" $(cat data/ecoli/ecoli_samples.txt | grep ERR10355911)
assert_equal "ERR10355944" $(cat data/ecoli/ecoli_samples.txt | grep ERR10355944)