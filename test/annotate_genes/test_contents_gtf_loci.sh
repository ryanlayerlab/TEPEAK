#!/bin/bash
set -eu
set -o pipefail 

gtf_file="output/ecoli/peak_0-10000/ecoli_0-10000_gtf.txt"
gtf_loci_file="output/ecoli/peak_0-10000/ecoli_0-10000_gtf_loci.txt"

while read -r -u 3 lineA && read -r -u 4 lineB; do 
    gene_nameA=$(echo $lineA | cut -d ";" -f 3 | cut -d " " -f 3 | tr -d '"')
    sample_nameA=$(echo $lineA | rev | cut -d " " -f 1 | rev)
    chromA=$(echo $lineA | cut -d " " -f 1)
    gene_idA=$(echo $lineA | cut -d " " -f 10 | tr -d '";')
    gene_typeA=$(echo $lineA | cut -d ";" -f 1 | cut -d " " -f 3)

    gene_nameB=$(echo $lineB | cut -d " " -f 6)
    sample_nameB=$(echo $lineB | cut -d " " -f 4)
    chromB=$(echo $lineB | cut -d " " -f 1)
    gene_idB=$(echo $lineB | cut -d " " -f 5)
    gene_typeB=$(echo $lineB | cut -d " " -f 7)

    [[ $gene_nameA =~ ^(ASAP:ABE-0005979|UniProtKB/Swiss-Prot:P64488)$ ]] || exit 1
    [[ $gene_nameB =~ ^(ASAP:ABE-0005979|UniProtKB/Swiss-Prot:P64488)$ ]] || exit 1
    [ $sample_nameA = $sample_nameB ] || exit 1
    [ "$chromA" = "$chromB" ] || exit 1
    [ "$chromA" = "NC_000913.3" ] || exit 1
    [ "$gene_idA" = "b1797" -a "$gene_idB" = "b1797" ] || exit 1
done 3<$gtf_file 4<$gtf_loci_file
echo true