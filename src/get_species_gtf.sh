#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

while getopts "s:d:f:" flag; do
  case "${flag}" in
    s) species="${OPTARG}" ;;
    d) species_dir="${OPTARG}" ;;
    f) zip_file="${OPTARG}" ;;
    *) exit 1 ;;
  esac
done

if [[ -z "${species:-}" || -z "${species_dir:-}" || -z "${zip_file:-}" ]]; then
  echo "Usage: $0 -s species -d species_dir -f zip" >&2
  exit 1
fi

mkdir -p "${species_dir}"

# Unzip directly into species_dir (creates species_dir/ncbi_dataset/…)
if [[ ! -f "${species_dir}/ncbi_dataset/data/dataset_catalog.json" ]]; then
  unzip -o "${zip_file}" -d "${species_dir}" > /dev/null
else
  echo "Dataset already extracted, proceeding..."
fi

gtf_path="$(find "${species_dir}/ncbi_dataset" -type f -name '*.gtf' | head -n 1 || true)"
if [[ -z "${gtf_path}" ]]; then
  echo "No .gtf found under ${species_dir}/ncbi_dataset" >&2
  exit 1
fi

cp "${gtf_path}" "${species_dir%/}/${species}.gtf"
