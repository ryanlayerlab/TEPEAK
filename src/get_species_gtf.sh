set -eu
set -o pipefail 

species=""
zipped_gtf_dataset=""

while getopts "s:f:" flag; do
    case "${flag}" in 
        s) species=${OPTARG};;
        f) zipped_gtf_dataset=${OPTARG};;
        *) echo "Usage: $0 -s species_name -f zipped_gtf_dataset"
            exit 1;;
    esac
done 

if [[ -z "$species" || -z "$zipped_gtf_dataset" ]]; then 
    echo "Missing arguments. Usage: $0 -s species_name -f zipped_gtf_dataset"
    exit 1
fi

data_dir=$(grep 'data_directory:' configs/config_${species}.yaml | awk '{print $2}')
echo "n" | unzip $zipped_gtf_dataset
# extracting the genomic.gtf refseq file
gcf_folder=$(ls ncbi_dataset/data | grep GCF)
gtf_file=$(ls "ncbi_dataset/data/$gcf_folder/")

# moving the gtf file to data_dir/species and renaming as species.gtf
mv "ncbi_dataset/data/$gcf_folder/$gtf_file" "$data_dir/$species/${species}.gtf"
rm -rf ncbi_dataset