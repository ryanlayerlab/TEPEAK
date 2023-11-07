#!/bin/bash
set -eu 
set -o pipefail
#name=$(cat /scratch/Shares/layer/workspace/devin_sra/sv_step/config.yaml | shyaml get-va$
#procs=$(cat /scratch/Shares/layer/workspace/devin_sra/sv_step/config.yaml | shyaml get-v$

if [ $# -eq 0 ]; then
    >&2 echo "No arguments provided"
    exit 1
fi

while getopts s:p:n:d:f: flag
do
    case "${flag}" in
        s) name=${OPTARG};;
        p) procs=${OPTARG};;
        n) threads=${OPTARG};;
        d) data_dir=${OPTARG};;
        f) sra_file=${OPTARG};;
    esac
done


sra_path="$data_dir"/"$sra_file"
lines=$(wc -l < $sra_path)
x=$((lines / procs))
cat $sra_path | xargs -I {} -L 2 ./call_insertions_parallel.sh -s $name -n $threads -d $data_dir -l {}

#-n $x -P $procs ./call_insertions_parallel.sh -l {}


#cat $sra_file | gargs --n $x -p $procs "sbatch combine.sbatch {}"
#cat $sra_file | gargs --n $x -p $procs "sbatch getTESeqs.sbatch {}"
