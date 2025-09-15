#!/bin/bash
set -euo pipefail

# usage: -s species -o output_dir -f sample_file
species=""; output_dir=""; sample_file=""
while getopts "s:o:f:" flag; do
  case "${flag}" in
    s) species="${OPTARG}" ;;
    o) output_dir="${OPTARG}" ;;
    f) sample_file="${OPTARG}" ;;
    *) echo "Usage: $0 -s species -o output_dir -f sample_file" >&2; exit 1 ;;
  esac
done

: "${species:?missing -s species}"
: "${output_dir:?missing -o output_dir}"
: "${sample_file:?missing -f sample_file}"

if ! command -v bcftools >/dev/null 2>&1; then
  echo "Error: bcftools not found in PATH." >&2
  exit 1
fi

mkdir -p "${output_dir}"
out_file="${output_dir}/${species}_global_vcf.txt"
: > "${out_file}"

while IFS= read -r line || [[ -n "$line" ]]; do
  # skip blanks and comments
  [[ -z "${line// }" || "${line#\#}" != "$line" ]] && continue
  sample="$line"

  vcf_gz="${output_dir}/${sample}/out.pass.vcf.gz"
  vcf_plain="${output_dir}/${sample}/out.pass.vcf"
  if   [[ -s "${vcf_gz}"   ]]; then src="${vcf_gz}"
  elif [[ -s "${vcf_plain}" ]]; then src="${vcf_plain}"
  else
    echo "Warning: no out.pass.vcf(.gz) for ${sample} under ${output_dir}/${sample} — skipping" >&2
    continue
  fi

  # CHROM  POS  END  SVLEN  SVINSSEQ  (no sample column here, by design)
  bcftools query -f "%CHROM\t%POS\t%INFO/END\t%SVLEN\t%SVINSSEQ\n" "${src}" >> "${out_file}"
done < "${sample_file}"
