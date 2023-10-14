import argparse
from optparse import OptionParser
import pandas as pd
import os
import json
from subprocess import Popen, PIPE
import getpass

parser = OptionParser()

parser.add_option("-f",
	dest="sra_file",
	help="Path to sra runinfo file")

parser.add_option("-n",
	dest="max_n",
	help="Max number of sra numbers returned")

parser.add_option("-o",
	dest="out_file",
	help="Output file")

(options, args) = parser.parse_args()

df = pd.read_csv (options.sra_file, low_memory = False)
sorted_df = df.sort_values(by='bases', ascending=False)
top_N_runs = sorted_df.head(int(options.max_n))['Run']
#
with open(options.out_file, 'w') as f:
    for run in top_N_runs:
        f.write(str(run) + '\n')


#SRR8766967
#SRR8767000 

#SRR8766982 - vcf parse error
























