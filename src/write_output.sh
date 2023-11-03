if [ $# -eq 0 ]; then
    >&2 echo "No arguments provided"
    exit 1
fi

while getopts s:l:u:g: flag
do
    case "${flag}" in
        s) species=${OPTARG};;
        l) low=${OPTARG};;
        u) high=${OPTARG};;
        g) gene=${OPTARG};;
    esac
done

range_file=output/"${species}"/peak_"$low"-"$high"/"${species}"_"$low"-"$high"_pop_vcf.txt

if [ "$g_arg" == "n" ]; then
    bedtools sort -i $range_file > output/"${species}"/peak_"$low"-"$high"/"${species}"_"$low"-"$high"_pop_vcf_sorted.txt

    sorted_range_file=output/"${species}"/peak_"$low"-"$high"/"${species}"_"$low"-"$high"_pop_vcf_sorted.txt 
    merged_range_file=output/"${species}"/peak_"$low"-"$high"/"${species}"_"$low"-"$high"_tmp.txt
    final_file=output/"${species}"/peak_"$low"-"$high"/"${species}"_"$low"-"$high"_merged.txt

    bedtools merge -i $sorted_range_file -c 4,5 -o collapse,collapse > $merged_range_file

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
else
    python3 src/output_helper.py -s $species -l $low -u $high
fi