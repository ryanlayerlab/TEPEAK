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
    peak_path = f'{output_dir}/peak_{low}-{high}'
    os.makedirs(peak_path, exist_ok=True)
    
    # Output files
    loci_annotated_file = f'{peak_path}/{species}_{low}-{high}_gtf_loci.txt'
    gtf_output_file = f'{peak_path}/{species}_{low}-{high}_gtf.txt'
    gene_annotate_file = f'{peak_path}/{species}_{low}-{high}_gene_annotate.txt'
    pop_vcf_sorted_file = f'{peak_path}/{species}_{low}-{high}_pop_vcf_sorted.txt'
    
    print("Starting gene annotation...")
    print(f"GTF file: {gtf_file}")
    print(f"Range file: {range_file}")
    
    # Check if input files exist and are not empty
    if not os.path.exists(gtf_file) or os.path.getsize(gtf_file) == 0:
        print(f"Error: GTF file {gtf_file} is missing or empty")
        # Create empty output files
        for output_file in [loci_annotated_file, gtf_output_file, gene_annotate_file, pop_vcf_sorted_file]:
            with open(output_file, 'w') as f:
                f.write("")
        return
    
    if not os.path.exists(range_file) or os.path.getsize(range_file) == 0:
        print(f"Error: Range file {range_file} is missing or empty")
        # Create empty output files
        for output_file in [loci_annotated_file, gtf_output_file, gene_annotate_file, pop_vcf_sorted_file]:
            with open(output_file, 'w') as f:
                f.write("")
        return
    
    # Run bedtools intersect
    cmd = f'bedtools intersect -a {range_file} -b {gtf_file} -wa -wb > {loci_annotated_file}'
    print(f"Running: {cmd}")
    exit_code = os.system(cmd)
    
    # Check if intersect produced results
    if not os.path.exists(loci_annotated_file) or os.path.getsize(loci_annotated_file) == 0:
        print("Warning: No intersections found between insertions and genes")
        print("This could mean:")
        print("  1. No insertions overlap with annotated genes")
        print("  2. GTF file format is incompatible")
        print("  3. Chromosome naming mismatch between VCF and GTF")
        
        # Create empty output files with headers for gene_annotate
        with open(gene_annotate_file, 'w') as f:
            f.write("chrom\tchromStart\tchromEnd\tseq\tsname\tgene_name\ttype\n")
        
        # Create empty other files
        for output_file in [gtf_output_file, pop_vcf_sorted_file]:
            with open(output_file, 'w') as f:
                f.write("")
        
        return
    
    try:
        # Process the intersection results
        processed_lines = []
        
        with open(loci_annotated_file, 'r') as f:
            lines = f.readlines()
            
        print(f"Processing {len(lines)} intersection records...")
        
        for line_num, line in enumerate(lines):
            try:
                line = line.strip()
                if not line:
                    continue
                    
                # Split the line into parts
                parts = line.split('\t')
                if len(parts) < 13:  # Need at least VCF (5) + GTF (8) columns
                    print(f"Warning: Line {line_num+1} has insufficient columns ({len(parts)}), skipping")
                    continue
                
                # Extract VCF parts (first 5 columns)
                vcf_part = '\t'.join(parts[:5])
                
                # Extract GTF attributes (last column of GTF part)
                gtf_attributes = parts[-1]
                
                # Extract gene information from GTF attributes
                gene_id_match = re.search(r'gene_id\s+"([^"]+)"', gtf_attributes)
                gene_id = gene_id_match.group(1) if gene_id_match else 'NA'
                
                # Try different attribute patterns for gene name
                gene_name = 'NA'
                
                # Standard patterns first
                for pattern in [r'gene_name\s+"([^"]+)"', r'db_xref\s+"([^"]+)"', r'Name\s*=\s*([^;]+)']:
                    match = re.search(pattern, gtf_attributes)
                    if match:
                        gene_name = match.group(1)
                        break
                
                # FIX 1: Extract gene symbol from transcript_id for TOGA GTF format
                if gene_name == 'NA':
                    transcript_id_match = re.search(r'transcript_id\s+"([^"]+)"', gtf_attributes)
                    if transcript_id_match:
                        transcript_id = transcript_id_match.group(1)
                        # Extract symbol between hashes: "ENST00000711184.1#PLCXD1#447" -> "PLCXD1"
                        symbol_match = re.search(r'#([^#]+)#', transcript_id)
                        if symbol_match:
                            gene_name = symbol_match.group(1)

                # Get feature type (3rd column in GTF part)
                gtf_start_idx = 5  # Where GTF columns start
                if len(parts) > gtf_start_idx + 2:
                    gene_type = parts[gtf_start_idx + 2]  # 3rd GTF column (feature type)
                else:
                    gene_type = 'unknown'
                
                # Format the output line
                processed_line = f"{vcf_part}\t{gene_id}\t{gene_name}\t{gene_type}"
                processed_lines.append(processed_line)
                
            except Exception as e:
                print(f"Warning: Error processing line {line_num+1}: {e}")
                continue
        
        # Write processed results to a temporary file
        temp_processed_file = f"{loci_annotated_file}.processed"
        with open(temp_processed_file, 'w') as f:
            for line in processed_lines:
                f.write(line + '\n')
        
        # Check if we have any processed results
        if len(processed_lines) == 0:
            print("Warning: No valid annotations found after processing")
            # Create empty output files
            with open(gene_annotate_file, 'w') as f:
                f.write("chrom\tchromStart\tchromEnd\tseq\tsname\tgene_name\ttype\n")
            
            for output_file in [gtf_output_file, pop_vcf_sorted_file]:
                with open(output_file, 'w') as f:
                    f.write("")
            return
        
        # Read the processed annotations
        df = pd.read_csv(temp_processed_file, sep='\t', header=None)
        
        # Assign column names based on expected structure
        expected_cols = ['chrom', 'chromStart', 'chromEnd', 'seq', 'sname', 'gene_id', 'gene_name', 'type']
        df.columns = expected_cols[:len(df.columns)]
        
        # Read the original range file
        df_full = pd.read_csv(range_file, sep='\t', header=None)
        df_full.columns = ['chrom','chromStart','chromEnd','seq','sname']
        
        # FIX 2: Include 'seq' in merge keys to prevent mis-attachment
        merge_cols = ['chrom', 'chromStart', 'chromEnd', 'seq', 'sname']
        merged_df = pd.merge(
            df_full, 
            df[merge_cols + ['gene_name', 'type']],
            on=merge_cols, 
            how='left'
        )
        
        # Fill missing values
        merged_df['gene_name'].fillna('None', inplace=True)
        merged_df['type'].fillna('None', inplace=True)
        
        # Save the final gene annotation file
        merged_df.to_csv(gene_annotate_file, sep='\t', index=False)
        print(f"Saved {len(merged_df)} records to {gene_annotate_file}")
        
        # Create the other required output files
        # Copy the original annotation results to gtf_loci.txt
        if os.path.exists(temp_processed_file):
            os.rename(temp_processed_file, gtf_output_file)
        
        # Create sorted pop_vcf file
        df_full_sorted = df_full.sort_values(['chrom', 'chromStart'])
        df_full_sorted.to_csv(pop_vcf_sorted_file, sep='\t', index=False, header=False)
        
        print(f"Gene annotation completed successfully")
        print(f"  Total insertions: {len(df_full)}")
        print(f"  Annotated insertions: {len(merged_df[merged_df['gene_name'] != 'None'])}")
        
    except pd.errors.EmptyDataError:
        print("Warning: Intersection file is empty - no genes overlap with insertions")
        # Create empty output files with headers
        with open(gene_annotate_file, 'w') as f:
            f.write("chrom\tchromStart\tchromEnd\tseq\tsname\tgene_name\ttype\n")
        
        for output_file in [gtf_output_file, pop_vcf_sorted_file]:
            with open(output_file, 'w') as f:
                f.write("")
        return
    
    except Exception as e:
        print(f"Error processing gene annotations: {e}")
        # Create empty output files
        for output_file in [gene_annotate_file, gtf_output_file, pop_vcf_sorted_file]:
            with open(output_file, 'w') as f:
                f.write("")
        return

if __name__ == '__main__':
    main()