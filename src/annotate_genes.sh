#!/bin/bash
if [ $# -eq 0 ]; then
    >&2 echo "No arguments provided"
    exit 1
fi

while getopts s:l:u: flag
do
    case "${flag}" in
        s) species=${OPTARG};;
        l) low=${OPTARG};;
        u) high=${OPTARG};;
    esac
done



data_dir=$(grep 'data_directory:' config_${species}.yaml | awk '{print $2}')
echo $data_dir
data_path="$(pwd)/$data_dir/${species}"
threads=$(grep 'threads:' config_${species}.yaml | awk '{print $2}')

sra_file=$data_dir/$species/${species}_samples.txt

range_file=output/"${species}"/peak_"$low"-"$high"/"${species}"_"$low"-"$high"_pop_vcf.txt

bedtools sort -i $range_file > output/"${species}"/peak_"$low"-"$high"/"${species}"_"$low"-"$high"_pop_vcf_sorted.txt

sorted_range_file=output/"${species}"/peak_"$low"-"$high"/"${species}"_"$low"-"$high"_pop_vcf_sorted.txt
gtf_file=output/"${species}"/peak_"$low"-"$high"/"${species}"_"$low"-"$high"_gtf.txt

bedtools intersect -a $data_dir/${species}.gtf -b $sorted_range_file -wb > $gtf_file

loci_file=output/"${species}"/peak_"$low"-"$high"/"${species}"_"$low"-"$high"_gtf_loci.txt

python3 gene_helper.py -s $species -l $low -u $high
