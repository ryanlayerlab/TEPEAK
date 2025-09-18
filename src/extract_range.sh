#!/bin/bash
# Extract insertions within a size range from each sample's VCF and
# append to a combined pop_vcf.txt. Handles .vcf.gz or .vcf.
set -euo pipefail

if [ $# -eq 0 ]; then
  echo "Usage: $0 -s species -f sample_file -o output_dir -l low -h high" >&2
  exit 1
fi

species=""; sample_file=""; output_dir=""; low=""; high=""
while getopts "s:f:o:l:h:" flag; do
  case "${flag}" in
    s) species="${OPTARG}" ;;
    f) sample_file="${OPTARG}" ;;
    o) output_dir="${OPTARG}" ;;
    l) low="${OPTARG}" ;;
    h) high="${OPTARG}" ;;
    *) echo "Usage: $0 -s species -f sample_file -o output_dir -l low -h high" >&2; exit 1 ;;
  esac
done

# basic checks
: "${species:?missing -s species}"
: "${sample_file:?missing -f sample_file}"
: "${output_dir:?missing -o output_dir}"
: "${low:?missing -l low}"
: "${high:?missing -h high}"

if ! command -v bcftools >/dev/null 2>&1; then
  echo "Error: bcftools not found in PATH." >&2
  exit 1
fi

# output locations
mkdir -p "${output_dir}"
peak_dir="${output_dir}/peak_${low}-${high}"
mkdir -p "${peak_dir}"
out_txt="${peak_dir}/${species}_${low}-${high}_pop_vcf.txt"

# start fresh
: > "${out_txt}"

# read each sample name and process
while IFS= read -r line || [ -n "$line" ]; do
  # skip blanks and comments
  case "${line}" in
    ""|\#*) continue ;;
  esac
  sra="${line}"

  vcfgz="${output_dir}/${sra}/out.pass.vcf.gz"
  vcf="${output_dir}/${sra}/out.pass.vcf"
  if [ -s "${vcfgz}" ]; then
    in_vcf="${vcfgz}"
  elif [ -s "${vcf}" ]; then
    in_vcf="${vcf}"
  else
    echo "Warning: no out.pass.vcf(.gz) for ${sra} under ${output_dir}/${sra} — skipping" >&2
    continue
  fi

  # filter by SVLEN and append with a sample column
  bcftools view -i "SVLEN>=${low} && SVLEN<=${high}" "${in_vcf}" \
  | bcftools query -f "%CHROM\t%POS\t%INFO/END\t%SVINSSEQ\t${sra}\n" \
  >> "${out_txt}"

done < "${sample_file}"
