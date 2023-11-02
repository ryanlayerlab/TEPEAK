from argparse import ArgumentParser
import os, subprocess, yaml

def parse_args():
    parser = ArgumentParser(description = "Process and run insurveyor.py on .bam files")
    parser.add_argument('-s', '--species', required = True, help = "Species name")
    return parser.parse_args()

def main(): 
    species = parse_args().species
    os.makedirs(f'output/{species}', exist_ok = True)    

    with open(f'configs/config_{species}.yaml') as config_file: 
        config_file = yaml.safe_load(config_file)
    
    threads = config_file.get('threads')
    data_dir = config_file.get('data_directory')
    reference_dir = os.path.join(data_dir, species, f'{species}.fa')
    cmd = ['insurveyor.py', 
            '--threads', str(threads), 
            'BAM_FILE', 
            'WORKDIR', 
            reference_dir
        ]
    
    sra_file = os.path.join(data_dir, species, f'{species}_samples.txt')
    with open(sra_file) as sra_file: 
        for line in sra_file: 
            sra_example = line.strip()
            output_dir = os.path.join('output', species, sra_example)
            bam_file = os.path.join(data_dir, species, f'{sra_example}.bam')

            os.makedirs(output_dir, exist_ok = True)
            cmd[3], cmd[4] = bam_file, output_dir 
            subprocess.run(cmd)

if __name__ == '__main__':
    main()