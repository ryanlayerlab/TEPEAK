import argparse
from optparse import OptionParser
import pandas as pd
import os
import json
from subprocess import Popen, PIPE
import getpass
import yaml


parser = OptionParser()

parser.add_option("-f",
	dest="sra_file",
	help="Path to sra runinfo file")

parser.add_option("-n",
	dest="max_n",
	help="Max number of sra numbers returned")

parser.add_option("-s",
	dest="species",
	help="species")

(options, args) = parser.parse_args()

config_file = 'config_'+options.species + '.yaml'

with open(config_file, 'r') as stream:
    try:
        data = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

data_dir = data['data_directory']


df = pd.read_csv (options.sra_file, low_memory = False)
sorted_df = df.sort_values(by='Bases', ascending=False)
top_N_runs = sorted_df.head(int(options.max_n))['Run']
#
with open(data_dir+'/'+options.species+'/'+options.species + '_samples.txt', 'w') as f:
    for run in top_N_runs:
        f.write(str(run) + '\n')


























