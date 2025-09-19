#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

while getopts "s:d:f:" flag; do
  case "${flag}" in
    s) species="${OPTARG}" ;;
    d) species_dir="${OPTARG}" ;;
    f) zip_file="${OPTARG}" ;;
    *) echo "Usage: $0 -s species -d species_dir -f zip" >&2; exit 1 ;;
  esac
done

if [[ -z "${species:-}" || -z "${species_dir:-}" || -z "${zip_file:-}" ]]; then
  echo "Usage: $0 -s species -d species_dir -f zip" >&2
  exit 1
fi

mkdir -p "${species_dir}"

# Extract to a temp dir so we don't care about the zip's internal layout
tmpd="$(mktemp -d "${species_dir%/}/.gtf_extract.XXXXXX")"
cleanup() { rm -rf "${tmpd}" || true; }
trap cleanup EXIT

unzip -oq "${zip_file}" -d "${tmpd}"

# Find a GTF anywhere inside (prefer uncompressed .gtf; fall back to .gtf.gz)
gtf_path="$(find "${tmpd}" -type f \( -name '*.gtf' -o -name '*.gtf.gz' \) | head -n 1 || true)"
if [[ -z "${gtf_path}" ]]; then
  echo "No .gtf or .gtf.gz found inside ${zip_file}" >&2
  exit 1
fi

out="${species_dir%/}/${species}.gtf"
if [[ "${gtf_path}" == *.gtf.gz ]]; then
  gunzip -c "${gtf_path}" > "${out}"
else
  cp -f "${gtf_path}" "${out}"
fi

echo "Wrote ${out}"
