import pandas as pd
import requests
from Bio.Align import PairwiseAligner
from argparse import ArgumentParser
import os.path

def parse_args():
    parser = ArgumentParser(description = "Process some arguments")
    parser.add_argument('-s', '--species', required = True, help = "species name")
    return parser.parse_args()

def construct_peak_sizes(df, min_window, max_window, window_size):
    peaks = []
    peak_sizes = []
    for i in range(min_window, max_window, window_size):
        min_i = i
        max_i = i + window_size
        b_rows = df.query('length >= @min_i')
        t_rows = b_rows.query('length <= @max_i') 
        t = t_rows['length'].value_counts()
        if len(t) > 0:
            mode = t_rows['length'].mode().values[0]
            mode_idx = t[mode]
            if mode_idx > 50:
                peaks.append(mode)
                peak_sizes.append(mode_idx)
    return peaks, peak_sizes

# check if the sequences are good match with each other
# if good cluster then keep one sequence with peak number
# will need error handle if the cluster is not a good cluster
# once all peaks are done then submit to dfam
def construct_peak_seq(df, peaks): 
    peak_seqs = []
    for peak in peaks: 
        print(peak)
        df_peak = df.loc[df['length'] == peak]
        samples = df_peak.sample(n=10)
        #print(samples['seq'].tolist())
        remaining_sequences = samples['seq'].tolist().copy()
        all_clusters = []
        for X in remaining_sequences:
            cluster = []
            X_seq = X
            for j, Y in enumerate(remaining_sequences):
                Y_seq = Y
                alignments = PairwiseAligner.align.globalxx(X_seq, Y_seq)
                if float(alignments[0].score) >= 0.75 * peak:
                    cluster.append(Y)
                    del remaining_sequences[j]
            if len(cluster) > 0:
                all_clusters.append(cluster)
        cluster_lens = [len(i) for i in all_clusters]
        largest_cluster = cluster_lens.index(max(cluster_lens))
        peak_seqs.append(all_clusters[largest_cluster][0])
    return peak_seqs

def construct_annotations(peak_seqs, url, params):
    annotations = []
    for peak_s in peak_seqs:
        params['sequence'] = peak_s
        peak_s_annots = []
        
        response = requests.post(url, data = params)
        print(response)
        print(response.json())
        resp_id = response.json()['id']
        print(resp_id)
        response = requests.get(url + resp_id)
        results = response.json()
        #print(results)
        if results['duration'] == 'Not finished':
            finished = False
            i = 0
            while not finished:
                print("query not finished - checking again in 5 seconds")
                time.sleep(5)
                response = requests.get(url + resp_id)
                results = response.json()
                if 'seconds' in results['duration']:
                    print("query finished")
                    finished = True
                    break
                i += 1
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
    return annotations

def main(): 
    options = parse_args()

    window_size = 50
    min_window = 200
    max_window = 6400

    species = options.species
    sv_info_file = os.path.join('output', species, f'{species}_global_vcf.txt')

    df = pd.read_csv(sv_info_file, sep='\t', lineterminator='\n')
    df.columns = ['chrom','start','end','length','seq']
    df = df[df['length']!='.']
    df['length'] = df['length'].astype(int)

    print(df.head())

    url = 'https://dfam.org/api/searches/'
    params = {'sequence' : 'PEAK_S', 'organism': 'Escherichia coli', 'cutoff' : 'curated'}
    # data = "sequence=atcg&cutoff=curated"
    #headers = {'Content-Type': 'application/x-www-form-urlencoded'}
    #headers = {'Content-Type': 'application/json'}
    peaks, peak_sizes = construct_peak_sizes(df, min_window, max_window, window_size)
    # print("peaks", peaks)
    # print("peak_sizes", peak_sizes)
    peak_seqs = construct_peak_seq(df, peaks)
    # print("peak_seqs", peak_seqs)
    annotations = construct_annotations(peak_seqs, url, params)
    # print("annotations", annotations)

    df_write = pd.DataFrame()
    df_write['peak_number'] = peaks
    df_write['peak_size'] = peak_sizes
    df_write['annotations'] = annotations
    df_write['peak_sequence'] = peak_seqs

    out_file = os.path.join('output', species, 'dfam_annotate.csv')
    df_write.to_csv(out_file)

if __name__ == '__main__':
    main()