import os, subprocess  # noqa: F821

def main(): 
    # noqa: F821  # snakemake provides the `snakemake` object
    threads = snakemake.params.threads
    ref_file = snakemake.input.ref
    species_dir = snakemake.params.species_dir
    sra_file = snakemake.input.sample_file
    output_dir = snakemake.params.output_dir
    os.makedirs(output_dir, exist_ok = True)    

    cmd = [
        'insurveyor.py',
        '--threads', str(threads),
        'BAM_FILE',
        'WORKDIR',
        ref_file
    ]

    # Run insurveyor per sample; ensure output exists even on failure
    with open(sra_file) as f:
        for line in f:
            sample = line.strip()
            if not sample:
                continue
            work_dir = os.path.join(output_dir, sample)
            bam_file = os.path.join(species_dir, f'{sample}.bam')
            os.makedirs(work_dir, exist_ok=True)
            # set BAM and working directory args
            cmd[3], cmd[4] = bam_file, work_dir
            try:
                print(f"Running InsurVeyor for sample {sample}...")
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError:
                print(f"Warning: insurveyor failed for sample {sample}, creating empty VCF")
            # ensure Snakemake sees an output file
            vcf_path = os.path.join(work_dir, 'out.pass.vcf.gz')
            if not os.path.exists(vcf_path):
                open(vcf_path, 'a').close()

if __name__ == '__main__':
    main()