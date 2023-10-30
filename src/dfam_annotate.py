import pandas as pd
import matplotlib.pyplot as plt
import requests
import yaml, json
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.pairwise2 import format_alignment
from io import StringIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import configparser
import argparse
from optparse import OptionParser
import os.path

parser = OptionParser()
parser.add_option("-s",
	dest="species",
	help="species name")

(options,args) = parser.parse_args()

window_size = 50
min_window = 200
max_window = 6400

species = options.species

with open('configs/config_'+species+'.yaml', 'r') as file:
    config_data = yaml.safe_load(file)
data_directory = config_data['data_directory']
# sv_info_file = 'output/'+ species +'/' + species + '_global_vcf.txt'
sv_info_file = os.path.join('output', species, f'{species}_global_vcf.txt')

df = pd.read_csv(sv_info_file, sep='\t', lineterminator='\n')
df.columns = ['chrom','start','end','length','seq']
df = df[df['length']!='.']
df['length'] = df['length'].astype(int)

print(df.head())
peaks = []
peak_sizes = []
for i in range(min_window,max_window,window_size):
    min_i = i
    max_i = i + window_size
    t_rows = df.query('length >= ' + str(min_i))
    t_rows = t_rows.query('length <= ' + str(max_i))   
    t = t_rows['length'].value_counts()
    if len(t) > 0:
        mode = t_rows['length'].mode().values[0]
        mode_idx = t[mode]
        if mode_idx > 50:
            peaks.append(mode)
            peak_sizes.append(mode_idx)


    # check if the sequences are good match with each other
    # if good cluster then keep one sequence with peak number
    # will need error handle if the cluster is not a good cluster
# once all peaks are done then submit to dfam

peak_seqs = []
for peak in peaks: 
    print(peak)
    df_peak = df.loc[df['length'] == peak]
    samples = df_peak.sample(n=10)
    #print(samples['seq'].tolist())
    remaining_sequences = samples['seq'].tolist().copy()
    all_clusters = []
    for i,X in enumerate(remaining_sequences):
        cluster = []
        X_seq = X
        for j,Y in enumerate(remaining_sequences):
            Y_seq = Y
            alignments = pairwise2.align.globalxx(X_seq, Y_seq)
            if float(alignments[0].score) >= 0.75 * peak:
                cluster.append(Y)
                del remaining_sequences[j]
        if len(cluster) > 0:
            all_clusters.append(cluster)
    cluster_lens = [len(i) for i in all_clusters]
    largest_cluster = cluster_lens.index(max(cluster_lens))
    peak_seqs.append(all_clusters[largest_cluster][0])



url = "https://dfam.org/api/POST/searches"
params = {'sequence' : 'ATCGATCG', 'cutoff' : 'curated'}
data = "sequence=atcg&cutoff=curated"
#headers = {'Content-Type': 'application/x-www-form-urlencoded'}
#headers = {'Content-Type': 'application/json'}
annotations = []
for peak_s in peak_seqs:
    peak_s_annots = []
    
    response = requests.post('https://dfam.org/api/searches', data={'sequence': peak_s, 'organism' : "Equus caballus",'cutoff': 'curated'})
    print(response)
    print(response.json())
    resp_id = response.json()['id']
    print(resp_id)
    response = requests.get('https://dfam.org/api/searches/'+resp_id)
    results = response.json()
    #print(results)
    if results['duration'] == 'Not finished':
        finished = False
        i=0
        while not finished:
            print("query not finished - checking again in 5 seconds")
            time.sleep(5)
            response = requests.get('https://dfam.org/api/searches/'+resp_id)
            results = response.json()
            if 'seconds' in results['duration']:
                print("query finished")
                finished = True
                break
            i+=1
            if i == 10:
                print('Timed out')
                exit()
        #no_query_hits = len(response.json()['results'][0]['hits'])
    #print(results)
    no_query_hits = len(response.json()['results'][0]['hits'])
    #print(no_query_hits)
    if no_query_hits == 0:
        peak_s_annots.append('No matching annotation found')
    for hit in range(no_query_hits):
        peak_s_annots.append((response.json()['results'][0]['hits'][hit]['description'],
                            response.json()['results'][0]['hits'][hit]['type']))
        #print(response.json()['results'][0]['hits'][hit]['description'])
    annotations.append(peak_s_annots)


df_write = pd.DataFrame()
df_write['peak_number'] = peaks
df_write['peak_size'] = peak_sizes
df_write['annotations'] = annotations
df_write['peak_sequence'] = peak_seqs

out_file = 'output/' + species + '/dfam_annotate.csv'
df_write.to_csv(out_file)
