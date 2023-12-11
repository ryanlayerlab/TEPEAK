import pandas as pd
import re, os

def main():
    species = snakemake.params.species
    gtf_file = snakemake.input.gtf_file 
    range_file = snakemake.input.range_file
    output_dir = snakemake.params.output_dir
    low, high = snakemake.params.low, snakemake.params.high

    peak_path = os.path.join(output_dir, f'peak_{low}-{high}')
    peak_species_filename = f'{species}_{low}-{high}'

    gcf_annotated_file = os.path.join(peak_path, f'{peak_species_filename}_gtf.txt')
    loci_annotated_file = os.path.join(peak_path, f'{peak_species_filename}_gtf_loci.txt')
    sorted_range_file = os.path.join(peak_path, f'{peak_species_filename}_pop_vcf_sorted.txt')
   
    os.system(f'bedtools sort -i {range_file} > {sorted_range_file}')
    os.system(f'bedtools intersect -a {gtf_file} -b {sorted_range_file} -wb > {gcf_annotated_file}')

    with open(gcf_annotated_file) as f, open(loci_annotated_file, 'w') as w:
        lines = f.readlines()
        for line in lines:
            try:
                gene_id = re.search('gene_id (.*?)(?:;|$)', line).group(1).strip('"') # extracts the value after 'gene_id' 
                gene_name = re.search(r'db_xref "(.*?)"', line).group(1).strip('"') # extracts the value after the first 'db_xref'
                gene_type = re.search('(.*?);', line).group(1).split('\t', 3)[2] # extracts value in the third column of line 
                s = re.search(r'([^;]*$)', line).group(1).split() # extracts the last item in line after a ';'

                s = '\t'.join(s[:3]) + '\t' + s[4]
                w.write('\t'.join([s, gene_id, gene_name, gene_type]) + '\n')
            except Exception as e:
                print(e)
                              
    df = pd.read_csv(loci_annotated_file, sep='\t', header=None)
                    
    header = ['chrom', 'chromStart', 'chromEnd','sname', 'gene_id', 'gene_name', 'type']
    df.columns = header[:len(df.columns)]
    df_full = pd.read_csv(range_file, sep='\t', lineterminator='\n')
    df_full.columns = ['chrom','chromStart','chromEnd','seq','sname']

    merged_df = pd.merge(
        df_full, df[['chrom', 'chromStart', 'chromEnd', 'sname', 'gene_name', 'type']],
        on = ['chrom', 'chromStart', 'chromEnd', 'sname'], 
        how = 'left'
    )
    merged_df['gene_name'].fillna('None', inplace=True)
    merged_df['type'].fillna('None', inplace=True)
    out_file = os.path.join(peak_path, f'{peak_species_filename}_gene_annotate.txt')
    merged_df.to_csv(out_file,  sep='\t',index=False)

if __name__ == '__main__':
    main()