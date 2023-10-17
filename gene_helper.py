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

gcf_annotated_file="output/"+options.species+"/peak_"+low+"-"+high+"/"+species+"_"+low+"-"+high+"_gtf.txt"
loci_annotated_file="output/"+options.species+"/peak_"+low+"-"+high+"/"+species+"_"+low+"-"+high+"_gtf_loci.txt"
pop_vcf_file="output/"+options.species+"/peak_"+low+"-"+high+"/"+species+"_"+low+"-"+high+"_pop_vcf.txt"

with open(gcf_annotated_file) as f:
    with open(loci_annotated_file, 'w') as w:
        lines = f.readlines()
        for line in lines:
            #print(name = re.search('gene_id (.*); gene_version', line).group(1))
            try:
                #name = re.search('gene_name (.*); gene_source', line).group(1)
                name = re.search('gene_id (.*); gene_version', line).group(1)
                gene_id =  name.strip('"')

                line = line.split(';')
                gene_name = line[2].split()
                gene_name = gene_name[1].strip('"')
                
                s = line[-1:][0].split()
                
                s = s[0] + '\t' + s[1] + '\t' + s[2] + '\t' + s[4]
                int_type = line[0].split()[2]

                w.write(s + '\t' + gene_id + '\t' + gene_name +'\t' +  int_type +'\n')
            except:
                pass
            try:
                #name = re.search('gene_id (.*); gene_version', line).group(1)
                name = re.search('gene_id (.*); gene_version', line).group(1)
                gene_id =  name.strip('"')
                line = line.split(';')
                gene_name = line[2].split()
                gene_name = gene_name[1].strip('"')
                s = line[-1:][0].split()
        
                s = s[0] + '\t' + s[1] + '\t' + s[2] + '\t' + s[4]
            
                int_type = line[0].split()[2]
                
                
                w.write(s + '\t' + gene_id + '\t' + gene_name +'\t' +  int_type +'\n')
            except:
                continue
                
#rint(loci_annotated_file)            
df = pd.read_csv(loci_annotated_file, sep='\t', comment='t', header=None)
                
header = ['chrom', 'chromStart', 'chromEnd','sname', 'gene_id', 'gene_name', 'type']
df.columns = header[:len(df.columns)]
#print(df)

df_full = pd.read_csv(pop_vcf_file, sep='\t', lineterminator='\n')
#print('here')
#print(df_full)
df_full.columns = ['chrom','chromStart','chromEnd','seq','sname']

merged_df = pd.merge(df_full, df[['chrom', 'chromStart', 'chromEnd', 'sname', 'gene_name', 'type']], 
                     on=['chrom', 'chromStart', 'chromEnd', 'sname'], 
                     how='left')
merged_df['gene_name'].fillna('None', inplace=True)
merged_df['type'].fillna('None', inplace=True)
out_file = "output/"+options.species+"/peak_"+low+"-"+high+"/"+species+"_"+low+"-"+high+"_gene_annotate.txt"
merged_df.to_csv(out_file,  sep='\t',index=False)
