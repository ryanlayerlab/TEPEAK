import pandas as pd
import re, os
import sys

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
                # FIX: Add null checks for all regex matches
                gene_id_match = re.search('gene_id (.*?)(?:;|$)', line)
                gene_id = gene_id_match.group(1).strip('"') if gene_id_match else 'NA'
                
                gene_name_match = re.search(r'db_xref "(.*?)"', line)
                gene_name = gene_name_match.group(1).strip('"') if gene_name_match else 'NA'
                
                # Extract value in the third column
                parts = line.split('\t')
                gene_type = parts[2] if len(parts) > 2 else 'NA'
                
                # Extract the last item in line after a ';'
                last_match = re.search(r'([^;]*$)', line)
                if last_match:
                    s = last_match.group(1).split()
                    if len(s) >= 5:
                        s = '\t'.join(s[:3]) + '\t' + s[4]
                    else:
                        continue  # Skip malformed lines
                else:
                    continue
                    
                w.write('\t'.join([s, gene_id, gene_name, gene_type]) + '\n')
            except Exception as e:
                print(f"Warning: Error parsing line: {e}")
                continue
    
    # FIX: Check if loci_annotated_file is empty or has no valid data
    if os.path.getsize(loci_annotated_file) == 0:
        print(f"Warning: No gene annotations found for range {low}-{high}")
        # Create empty output files with proper structure
        out_file = os.path.join(peak_path, f'{peak_species_filename}_gene_annotate.txt')
        
        # Read the full range file
        df_full = pd.read_csv(range_file, sep='\t', lineterminator='\n')
        df_full.columns = ['chrom','chromStart','chromEnd','seq','sname']
        
        # Add empty gene annotation columns
        df_full['gene_name'] = 'None'
        df_full['type'] = 'None'
        
        df_full.to_csv(out_file, sep='\t', index=False)
        
        # Create empty other output files
        for suffix in ['gtf_loci.txt', 'gtf.txt', 'pop_vcf_sorted.txt']:
            empty_file = os.path.join(peak_path, f'{peak_species_filename}_{suffix}')
            if not os.path.exists(empty_file) or suffix == 'gtf_loci.txt':
                open(empty_file, 'a').close()
        
        sys.exit(0)
                              
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