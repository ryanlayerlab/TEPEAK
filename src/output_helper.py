import pandas as pd
import os.path

def main():
    species = snakemake.params.species
    output_dir = snakemake.params.output_dir
    annotated_file = snakemake.input.annotated_file
    low, high = snakemake.params.low, snakemake.params.high

    peak_path = os.path.join(output_dir, f'peak_{low}-{high}')
    peak_species_filename = f'{species}_{low}-{high}'

    # Read the BED file into a DataFrame
    df = pd.read_csv(annotated_file, sep="\t", header=0, names=["Chrom", "start", "end", "sequence", "sampleID", "gene", "type"])

    # Sort the DataFrame
    df = df.sort_values(by=["Chrom", "start", "end"])

    # Merge overlapping rows
    merged_rows = []
    current_row = df.iloc[0].copy()
    max_seq_len = len(current_row["sequence"])
    sample_ids = [current_row["sampleID"]]

    for _, row in df.iloc[1:].iterrows():
        if row["Chrom"] == current_row["Chrom"] and row["start"] <= current_row["end"]:
            # Adjust end position based on sequence length
            end_pos = current_row["start"] + max_seq_len
            current_row["end"] = min(end_pos, row["end"])
            sample_ids.append(row["sampleID"])
            max_seq_len = max(max_seq_len, len(row["sequence"]))
        else:
            current_row["sampleID"] = ",".join(sample_ids)
            merged_rows.append(current_row)
            current_row = row
            max_seq_len = len(row["sequence"])
            sample_ids = [row["sampleID"]]

    # Handle the last row
    current_row["sampleID"] = ",".join(sample_ids)
    merged_rows.append(current_row)

    # Convert to a new DataFrame
    merged_df = pd.DataFrame(merged_rows)
    out_file = os.path.join(peak_path, f'{peak_species_filename}_merged_genes.txt')
    # Write to a new BED file
    merged_df.to_csv(out_file, sep="\t", header=False, index=False)

if __name__ == '__main__':
    main()