import os, subprocess, argparse, yaml

def parse_args():
    parser = argparse.ArgumentParser(description = "Process and run insurveyor.py on .bam files")
    parser.add_argument('-s', '--species', required = True, help = "Species name")
    return parser.parse_args()

def main(): 
    species = parse_args().species
    os.makedirs(f'output/{species}', exist_ok = True)    

    with open(f'configs/config_{species}.yaml') as config_file: 
        config_file = yaml.safe_load(config_file)
    
    threads = config_file.get('threads')
    bam_dir = config_file.get('bam_dir')
    data_dir = config_file.get('data_directory') if bam_dir is None else config_file.get('data_directory')

    sra_file = os.path.join(data_dir, species, species + '_samples.txt')
    with open(sra_file) as sra_file: 
        for line in sra_file: 
            sra_example = line.strip()
            output_dir = os.path.join('output', species, sra_example)
            os.makedirs(output_dir, exist_ok = True)
            cmd = [
                'insurveyor.py', 
                '--threads', str(threads), 
                os.path.join(data_dir, sra_example + '.bam'), 
                output_dir, 
                os.path.join(data_dir, species + '.fa')
            ]
            print(" ".join(cmd))
            subprocess.run(cmd)

if __name__ == '__main__':
    main()