import pandas as pd
import re, os
import sys

def main():
    # Get file paths from Snakemake
    gtf_file = snakemake.input.gtf_file
    range_file = snakemake.input.range_file
    species = snakemake.params.species
    output_dir = snakemake.params.output_dir
    low = snakemake.params.low
    high = snakemake.params.high
    
    # Create output directory
    os.makedirs(f'{output_dir}/peak_{low}-{high}', exist_ok=True)
    
    # Output files
    loci_annotated_file = f'{output_dir}/peak_{low}-{high}/{species}_{low}-{high}_gtf_loci.txt'
    gtf_output_file = f'{output_dir}/peak_{low}-{high}/{species}_{low}-{high}_gtf.txt'
    gene_annotate_file = f'{output_dir}/peak_{low}-{high}/{species}_{low}-{high}_gene_annotate.txt'
    
    print("Starting gene annotation...")
    print(f"GTF file: {gtf_file}")
    print(f"Range file: {range_file}")
    
    # Check if input files exist and are not empty
    if not os.path.exists(gtf_file) or os.path.getsize(gtf_file) == 0:
        print(f"Error: GTF file {gtf_file} is missing or empty")
        # Create empty output files
        for output_file in [loci_annotated_file, gtf_output_file, gene_annotate_file]:
            with open(output_file, 'w') as f:
                f.write("")
        return
    
    if not os.path.exists(range_file) or os.path.getsize(range_file) == 0:
        print(f"Error: Range file {range_file} is missing or empty")
        # Create empty output files
        for output_file in [loci_annotated_file, gtf_output_file, gene_annotate_file]:
            with open(output_file, 'w') as f:
                f.write("")
        return
    
    # Run bedtools intersect
    cmd = f'bedtools intersect -a {range_file} -b {gtf_file} -wa -wb > {loci_annotated_file}'
    print(f"Running: {cmd}")
    os.system(cmd)
    
    # Check if intersect produced results
    if not os.path.exists(loci_annotated_file) or os.path.getsize(loci_annotated_file) == 0:
        print("Warning: No intersections found between insertions and genes")
        print("This could mean:")
        print("  1. No insertions overlap with annotated genes")
        print("  2. GTF file format is incompatible")
        print("  3. Chromosome naming mismatch between VCF and GTF")
        
        # Create empty output files with headers
        with open(gene_annotate_file, 'w') as f:
            f.write("chrom\tstart\tend\tlength\tseq\tgene_name\tgene_type\n")
        
        with open(gtf_output_file, 'w') as f:
            f.write("")
            
        return
    
    try:
        # Read the intersect results
        df = pd.read_csv(loci_annotated_file, sep='\t', header=None)
        print(f"Found {len(df)} intersection records")
        
        # Continue with processing...
        # ...existing code for processing intersections...
        
    except pd.errors.EmptyDataError:
        print("Warning: Intersection file is empty - no genes overlap with insertions")
        # Create empty output files with headers
        with open(gene_annotate_file, 'w') as f:
            f.write("chrom\tstart\tend\tlength\tseq\tgene_name\tgene_type\n")
        
        with open(gtf_output_file, 'w') as f:
            f.write("")
        
        return

    with open(gtf_annotated_file) as f, open(loci_annotated_file, 'w') as w:
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