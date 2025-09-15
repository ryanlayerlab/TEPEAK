# 1) Current helper scripts (the ones those rules run)
sed -n '1,200p' src/get_global_vcf.sh   > /tmp/get_global_vcf.sh.txt
sed -n '1,200p' src/extract_range.sh    > /tmp/extract_range.sh.txt

# 2) Sample list used by both rules
cat horse_samples.txt                   > /tmp/horse_samples.txt

# 3) Output tree + sizes (portable on macOS)
mkdir -p /tmp/outdump
find output/horse -maxdepth 3 -type f 2>/dev/null \
  -exec stat -f "%N\t%z" {} \; | sort   > /tmp/output_horse_files.txt
ls -ld output/horse output/horse/peak_0-10000 2>/dev/null \
                                        > /tmp/output_dirs.txt

# 4) For the first sample, show what files exist and a header peek
SAMPLE="$(head -n1 horse_samples.txt | tr -d '[:space:]')"
ls -lh "output/horse/${SAMPLE}/" 2>&1  > "/tmp/${SAMPLE}_dir.txt"

# bcftools env
which bcftools                         > /tmp/bcftools_which.txt 2>&1
bcftools --version                    >> /tmp/bcftools_which.txt 2>&1

# VCF headers (gz & plain if present)
if [ -s "output/horse/${SAMPLE}/out.pass.vcf.gz" ]; then
  bcftools view -h "output/horse/${SAMPLE}/out.pass.vcf.gz" | head -n 40 \
    > "/tmp/${SAMPLE}_out.pass.vcfgz.header.txt"
fi
if [ -s "output/horse/${SAMPLE}/out.pass.vcf" ]; then
  head -n 40 "output/horse/${SAMPLE}/out.pass.vcf" \
    > "/tmp/${SAMPLE}_out.pass.vcf.head.txt"
fi

# 5) Dry-run (trace) each script manually to capture their stderr/stdout
#    (This helps if a path is missing or a bcftools call fails)
bash -x src/get_global_vcf.sh  -s horse -o output/horse -f horse_samples.txt \
  > /tmp/get_global_vcf.run.out 2>&1 || true

bash -x src/extract_range.sh   -s horse -f horse_samples.txt -o output/horse -l 0 -h 10000 \
  > /tmp/extract_range.run.out 2>&1 || true

# 6) Latest Snakemake log tail for context
LOG="$(ls -t .snakemake/log/*.log 2>/dev/null | head -n1)"
[ -n "$LOG" ] && tail -n 200 "$LOG"   > /tmp/snakemake.tail.txt || true

# 7) Bash version (useful for script semantics)
bash --version | head -n1            > /tmp/bash_version.txt

