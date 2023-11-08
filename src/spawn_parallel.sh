#!/bin/bash
set -eu 
set -o pipefail

if [ $# -eq 0 ]; then
    >&2 echo "No arguments provided"
    exit 1
fi

while getopts s:p:n: flag
do
    case "${flag}" in
        s) species=${OPTARG};;
        p) procs=${OPTARG};;
        n) threads=${OPTARG};;
        *) echo "Usage: $0 -s species_name -p number_of_processes -n number_of_threads_per_process"
           exit 1;;
    esac
done

data_dir=$(grep 'data_directory:' configs/config_${species}.yaml | awk '{print $2}')
sra_file="$data_dir/${species}_samples.txt"
cat $sra_file | xargs -I {} -L 2 ./call_insertions_parallel.sh -s $species -n $threads -d $data_dir -l {}