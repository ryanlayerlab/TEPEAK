import os, subprocess

def main(): 
    threads = snakemake.params.threads
    ref_file = snakemake.input.ref
    species_dir = snakemake.params.species_dir
    sra_file = snakemake.input.sample_file
    output_dir = snakemake.params.output_dir
    os.makedirs(output_dir, exist_ok = True)    

    cmd = ['insurveyor.py', 
            '--threads', str(threads), 
            'BAM_FILE', 
            'WORKDIR', 
            ref_file
        ]
    
    with open(sra_file) as sra_file: 
        for line in sra_file: 
            sra_example = line.strip()
            work_dir = os.path.join(output_dir, sra_example)
            bam_file = os.path.join(species_dir, f'{sra_example}.bam')

            os.makedirs(work_dir, exist_ok = True)
            cmd[3], cmd[4] = bam_file, work_dir 
            subprocess.run(cmd)

if __name__ == '__main__':
    main()