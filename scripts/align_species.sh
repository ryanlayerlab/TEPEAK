#!/usr/bin/env bash


IFS=$'\n\t'
shopt -s nullglob

species=""
data_path=""
threads="1"
fastq_dir=""
sample_file=""

while getopts "s:d:t:f:l:" flag; do
  case "${flag}" in
    s) species="${OPTARG}" ;;
    d) data_path="${OPTARG}" ;;
    t) threads="${OPTARG}" ;;
    f) fastq_dir="${OPTARG}" ;;   # if set => FASTQ mode
    l) sample_file="${OPTARG}" ;;
    *) exit 1 ;;
  esac
done

if [[ -z "${species}" || -z "${data_path}" || -z "${sample_file}" ]]; then
  echo "Missing required args: -s species -d data_path -l sample_file" >&2
  exit 1
fi

PICARD_JAR="$(dirname "$0")/../picard/build/libs/picard.jar"
ref="${data_path}/${species}.fa"

declare -a R1_LIST=()
declare -a R2_LIST=()

pair_for_r1() {
  local r1="$1" cand=""
  [[ -z "$cand" && "$r1" == *"_R1."*  ]] && cand="${r1/_R1./_R2.}"
  [[ -z "$cand" && "$r1" == *"_R1_"*  ]] && cand="${r1/_R1_/_R2_}"
  [[ -z "$cand" && "$r1" == *"_1."*   ]] && cand="${r1/_1./_2.}"
  [[ -z "$cand" && "$r1" == *"_1_"*   ]] && cand="${r1/_1_/_2_}"
  [[ -z "$cand" && "$r1" == *".R1."*  ]] && cand="${r1/.R1./.R2.}"
  [[ -z "$cand" && "$r1" == *".1."*   ]] && cand="${r1/.1./.2.}"
  [[ -z "$cand" && "$r1" == *"-R1."*  ]] && cand="${r1/-R1./-R2.}"
  [[ -z "$cand" && "$r1" == *"-1."*   ]] && cand="${r1/-1./-2.}"
  [[ -n "$cand" && -f "$cand" ]] && { printf '%s\n' "$cand"; return 0; }
  printf '%s\n' ""
}

arr_contains() { local n="$1"; shift; local x; for x in "$@"; do [[ "$x" == "$n" ]] && return 0; done; return 1; }

safe_append_pair() {
  local r1="$1" r2="$2"
  arr_contains "$r1" "${R1_LIST[@]}" || { R1_LIST+=("$r1"); R2_LIST+=("$r2"); }
}

collect_pairs_for_sample() {
  local sample="$1" dir="$2"
  R1_LIST=(); R2_LIST=()
  local ext r1 r2
  for ext in fastq.gz fq.gz fastq fq; do
    for r1 in \
      ${dir}/${sample}*_R1_*.${ext} \
      ${dir}/${sample}*_R1.${ext} \
      ${dir}/${sample}*_1_*.${ext} \
      ${dir}/${sample}*_1.${ext} \
      ${dir}/${sample}.R1.${ext} \
      ${dir}/${sample}.1.${ext} \
      ${dir}/${sample}-R1.${ext} \
      ${dir}/${sample}-1.${ext}
    do
      [[ -e "$r1" ]] || continue
      r2="$(pair_for_r1 "$r1")"
      [[ -n "$r2" ]] || continue
      safe_append_pair "$r1" "$r2"
    done
  done
}

concat_to_fifo() {
  local fifo="$1"; shift
  (
    local f
    for f in "$@"; do
      if [[ "$f" == *.gz ]]; then gzip -cd "$f"; else cat "$f"; fi
    done
  ) > "$fifo"
}

align_from_fastq() {
  local sample="$1"
  collect_pairs_for_sample "$sample" "$fastq_dir"

  if (( ${#R1_LIST[@]} == 0 || ${#R2_LIST[@]} == 0 )); then
    echo "No FASTQ pairs found for ${sample} in ${fastq_dir}" >&2
    exit 1
  fi
  if (( ${#R1_LIST[@]} != ${#R2_LIST[@]} )); then
    echo "Unequal R1/R2 pair counts for ${sample}" >&2
    exit 1
  fi


  local tmpd r1_fifo r2_fifo pid1 pid2
  tmpd="$(mktemp -d)"; r1_fifo="${tmpd}/r1.fq"; r2_fifo="${tmpd}/r2.fq"
  mkfifo "$r1_fifo" "$r2_fifo"

  concat_to_fifo "$r1_fifo" "${R1_LIST[@]}" & pid1=$!
  concat_to_fifo "$r2_fifo" "${R2_LIST[@]}" & pid2=$!

  bwa mem -M -t "${threads}" -R "@RG\tID:${sample}\tSM:${sample}" "${ref}" "$r1_fifo" "$r2_fifo" \
    | samtools sort -@ "${threads}" -o "${data_path}/${sample}.bam" -

  wait "$pid1" "$pid2" || true
  rm -f "$r1_fifo" "$r2_fifo"; rmdir "$tmpd" || true

  samtools index "${data_path}/${sample}.bam"

  java -jar "${PICARD_JAR}" AddOrReplaceReadGroups \
    I="${data_path}/${sample}.bam" \
    O="${data_path}/${sample}.rg.bam" \
    RGID="${sample}" RGLB="lib1" RGPL="ILLUMINA" RGPU="unit1" RGSM="${sample}"

  samtools index "${data_path}/${sample}.rg.bam"
  mv "${data_path}/${sample}.rg.bam"     "${data_path}/${sample}.bam"
  mv "${data_path}/${sample}.rg.bam.bai" "${data_path}/${sample}.bam.bai"


  local mapped
  mapped=$(samtools view -c -F 4 "${data_path}/${sample}.bam")
  if [[ "${mapped}" -eq 0 ]]; then
    echo "No mapped reads in ${sample}.bam — check FASTQs (pairing/counts/species)." >&2
    exit 1
  fi
}

align_from_sra() {
  local sample="$1"
  local tmpdir r1 r2 mapped
  tmpdir="$(mktemp -d)"
  fasterq-dump --split-files -e "${threads}" -O "${tmpdir}" "${sample}" >/dev/null
  r1="${tmpdir}/${sample}_1.fastq"; r2="${tmpdir}/${sample}_2.fastq"
  if [[ ! -s "$r1" || ! -s "$r2" ]]; then
    echo "SRA download produced empty reads for ${sample}" >&2
    exit 1
  fi
  bwa mem -M -t "${threads}" -R "@RG\tID:${sample}\tSM:${sample}" "${ref}" "$r1" "$r2" \
    | samtools sort -@ "${threads}" -o "${data_path}/${sample}.bam" -
  rm -rf "${tmpdir}"
  samtools index "${data_path}/${sample}.bam"

  java -jar "${PICARD_JAR}" AddOrReplaceReadGroups \
    I="${data_path}/${sample}.bam" \
    O="${data_path}/${sample}.rg.bam" \
    RGID="${sample}" RGLB="lib1" RGPL="ILLUMINA" RGPU="unit1" RGSM="${sample}"

  samtools index "${data_path}/${sample}.rg.bam"
  mv "${data_path}/${sample}.rg.bam"     "${data_path}/${sample}.bam"
  mv "${data_path}/${sample}.rg.bam.bai" "${data_path}/${sample}.bam.bai"

  mapped=$(samtools view -c -F 4 "${data_path}/${sample}.bam")
  if [[ "${mapped}" -eq 0 ]]; then
    echo "No mapped reads in ${sample}.bam — check SRA run (pairing/species)." >&2
    exit 1
  fi
}

# --- main ---
if [[ -n "${fastq_dir}" ]]; then
  while IFS= read -r sample || [[ -n "${sample}" ]]; do
    [[ -z "${sample// }" ]] && continue
    align_from_fastq "${sample}"
  done < "${sample_file}"
else
  while IFS= read -r sample || [[ -n "${sample}" ]]; do
    [[ -z "${sample// }" ]] && continue
    align_from_sra "${sample}"
  done < "${sample_file}"
fi
