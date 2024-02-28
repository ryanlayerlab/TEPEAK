import pandas as pd
import matplotlib.pyplot as plt

def main(): 
    lower = 0 if snakemake.params.low is None else snakemake.params.low
    upper = 10000 if snakemake.params.high is None else snakemake.params.high
    min_size, max_size = str(lower), str(upper)

    sv_info_file = snakemake.input.global_vcf_file

    df = pd.read_csv(sv_info_file, sep='\t', lineterminator='\n', header=None)
    df.columns = ['chrom', 'start', 'end', 'length']

    df = df[df['length'] != '.']
    df.dropna(subset = ['length'], inplace = True)
    df['length'] = df['length'].astype(int)

    t_rows = df.query('length >= ' + min_size)
    t_rows = t_rows.query('length <= ' + max_size)
    t = df['length'].value_counts()

    plt.hist(t_rows['length'], density=False, bins=len(t))
    plt.ylabel('Frequency')
    plt.xlabel('Deletion size (bp)')
    plt.show()

    plt.savefig(snakemake.output.plot)

if __name__ == '__main__':
    main()