import argparse
from optparse import OptionParser
import pandas as pd
import os
import json
from subprocess import Popen, PIPE
import getpass
import yaml
import re

parser = OptionParser()

parser.add_option("-s",
	dest="species",
	help="species name")

parser.add_option("-l",
	dest="lower",
	help="")

parser.add_option("-u",
	dest="upper",
	help="")

(options, args) = parser.parse_args()

config_file = 'config_'+options.species + '.yaml'

with open(config_file, 'r') as stream:
    try:
        data = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

low = options.lower
high = options.upper
species =options.species
print(species)
data_dir = data['data_directory']

annotated_file="output/"+options.species+"/peak_"+low+"-"+high+"/"+species+"_"+low+"-"+high+"_gene_annotate.txt"
#loci_annotated_file="output/"+options.species+"/peak_"+low+"-"+high+"/"+species+"_"+low+"-"+high+"_gtf_loci.txt"
#pop_vcf_file="output/"+options.species+"/peak_"+low+"-"+high+"/"+species+"_"+low+"-"+high+"_pop_vcf.txt"


# Read the BED file into a DataFrame
df = pd.read_csv(annotated_file, sep="\t", header=None, names=["Chrom", "start", "end", "sequence", "sampleID", "gene", "type"])

# Sort the DataFrame
df = df.sort_values(by=["Chrom", "start", "end"])

# Merge overlapping rows
merged_rows = []
current_row = df.iloc[0]
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
out_file = "output/"+options.species+"/peak_"+low+"-"+high+"/"+species+"_"+low+"-"+high+"_genes_merged.txt"
# Write to a new BED file
merged_df.to_csv(out_file, sep="\t", header=False, index=False)