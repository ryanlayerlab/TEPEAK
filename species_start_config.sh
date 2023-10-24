#!/bin/bash

# Initialize variables
species=""
data_dir=""
threads=""

# Parse command-line arguments
while getopts "s:d:n:" flag; do
    case "${flag}" in
        s) species=${OPTARG};;
        d) data_dir=${OPTARG};;
        n) threads=${OPTARG};;
        *) echo "Usage: $0 -s species_name -d data_directory -n number_of_threads"
           exit 1
           ;;
    esac
done

# Check if all flags are provided
if [[ -z "$species" || -z "$data_dir" || -z "$threads" ]]; then
    echo "Missing arguments. Usage: $0 -s species_name -d data_directory -n number_of_threads"
    exit 1
fi

mkdir -p $data_dir
mkdir -p $data_dir/$species

# Create the config_<species>.yaml file
config_file="configs/config_${species}.yaml"

echo "species: $species" > $config_file
echo "data_directory: $data_dir" >> $config_file
echo "threads: $threads" >> $config_file

echo "Config file $config_file created successfully!"

