from argparse import ArgumentParser
import pandas as pd
import os.path

def parse_args():
    parser = ArgumentParser(description = "Process some arguments")
    parser.add_argument('-s', '--species', required = True, help = "species name")
    parser.add_argument('-l', '--lower', required = True, help = "lower bound")
    parser.add_argument('-u', '--upper', required = True, help = "upper bound")
    return parser.parse_args()

def main():
    args = parse_args()
    species = args.species
    low, high = args.lower, args.upper

    peak_path = os.path.join('output', species, f'peak_{low}-{high}')
    peak_species_filename = f'{species}_{low}-{high}'

    annotated_file = os.path.join(peak_path, f'{peak_species_filename}_gene_annotate.txt')
    # Read the BED file into a DataFrame
    df = pd.read_csv(annotated_file, sep="\t", header=None, names=["Chrom", "start", "end", "sequence", "sampleID", "gene", "type"])

    # Sort the DataFrame
    df = df.sort_values(by=["Chrom", "start", "end"])

    # Merge overlapping rows
    merged_rows = []
    current_row = df.iloc[0]
    max_seq_len = len(current_row["sequence"])
    sample_ids = [current_row["sampleID"]]

    for row in df.iloc[1:]:
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
    out_file = os.path.join(peak_path, f'{peak_species_filename}_genes_merged.txt')
    # Write to a new BED file
    merged_df.to_csv(out_file, sep="\t", header=False, index=False)

if __name__ == '__main__':
    main()