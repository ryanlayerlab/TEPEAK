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
# threads=$(grep 'threads:' config_${species}.yaml | awk '{print $2}')

sra_file=$data_dir/$species/${species}_samples.txt


mkdir output/"${species}"/
mkdir output/"${species}"/peak_"$low"-"$high"/

rm -f output/"${species}"/peak_"$low"-"$high"/output_seqs.vcf

while read -r line; do
        sra_example=$line

        bcftools view -i  'SVLEN>='${low}' && SVLEN<='${high}'' output/"${species}"/"${sra_example}"/out.pass.vcf -o  \
		output/"${species}"/peak_"$low"-"$high"/output_seqs.vcf

        bcftools query -f '%CHROM\t%POS\t%INFO/END\t%SVINSSEQ\t'${sra_example}'\n' output/"${species}"/peak_"$low"-"$high"/output_seqs.vcf >> \
                output/"${species}"/peak_"$low"-"$high"/"${species}"_"$low"-"$high"_pop_vcf.txt

       rm -f output/"${species}"/peak_"$low"-"$high"/output_seqs.vcf

done < "$sra_file"



