#!/bin/bash
if [ $# -eq 0 ]; then
    >&2 echo "No arguments provided"
    exit 1
fi

while getopts s:o:r:l:h: flag; do
    case "${flag}" in
        s) species=${OPTARG};;
        o) output_dir=${OPTARG};;
        r) range_file=${OPTARG};;
        l) low=${OPTARG};;
        h) high=${OPTARG};;
        *) echo echo "Usage: $0 -s species_name -o output_dir -r range_file -l low -h high"
        exit 1;;
    esac
done

sorted_vcf_file="$output_dir/peak_$low-$high/${species}_$low-${high}_pop_vcf_sorted.txt"
bedtools sort -i $range_file > $sorted_vcf_file

merged_range_file="$output_dir/peak_$low-$high/${species}_$low-${high}_tmp.txt"
final_file="$output_dir/peak_$low-$high/${species}_$low-${high}_merged.txt"

bedtools merge -i $sorted_vcf_file -c 4,5 -o collapse,collapse > $merged_range_file

awk 'BEGIN {FS=OFS="\t"}
{
    split($4, sequences, ",");
    split($5, sampleIDs, ",");
    
    max_seq_len = 0;
    for (seq in sequences) {
        len = length(sequences[seq]);
        if (len > max_seq_len) {
            max_seq_len = len;
        }
    }

    if (($3 - $2) > max_seq_len) {
        $3 = $2 + max_seq_len;
    }

    print $1, $2, $3, join(sampleIDs, ",");
}
function join(arr, sep) {
    result = arr[1];
    for (i=2; i<=length(arr); i++) {
        result = result sep arr[i];
    }
    return result;
}
' $merged_range_file > $final_file
rm $merged_range_file
