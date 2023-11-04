from argparse import ArgumentParser
import pandas as pd
import yaml, os.path

def parse_args():
    parser = ArgumentParser(description = "Process some arguments")
    parser.add_argument('-f', '--sra_file', help = "Path to sra runinfo file")
    parser.add_argument('-n', '--max_n', help = "Max number of sra numbers returned")
    parser.add_argument('-s', '--species', help = "species")
    return parser.parse_args()

def main():
    args = parse_args()
    species = args.species
    config_file = f'configs/config_{species}.yaml'

    with open(config_file, 'r') as stream:
        data = yaml.safe_load(stream)
    data_dir = data['data_directory']

    df = pd.read_csv (args.sra_file, low_memory = False)
    sorted_df = df.sort_values(by = 'Bases', ascending = False)
    top_N_runs = sorted_df.head(int(args.max_n))['Run']

    with open(os.path.join(data_dir, species, f'{species}_samples.txt'), 'w') as f:
        for run in top_N_runs:
            f.write(str(run) + '\n')

if __name__ == '__main__':
    main()