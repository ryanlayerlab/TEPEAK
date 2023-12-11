import pandas as pd

def main():
    df = pd.read_csv(snakemake.input.sra_file, low_memory = False)
    sorted_df = df.sort_values(by = 'Bases', ascending = False)
    max_n = snakemake.params.max_n
    top_N_runs = sorted_df.head(int(max_n))['Run']

    with open(snakemake.output.sample_file, 'w') as f:
        for run in top_N_runs:
            f.write(f'{run}\n')

if __name__ == '__main__':
    main()