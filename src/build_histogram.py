from argparse import ArgumentParser
import pandas as pd
import matplotlib.pyplot as plt

def parse_args():
	parser = ArgumentParser(description = "Process some arguments.")
	parser.add_argument('-f', '--pop_vcf_file', required = True, help = "Population VCF file")
	parser.add_argument('-l', '--lower', required = False, help = "lower range")
	parser.add_argument('-u', '--upper', required = False, help = "upper range") 
	return parser.parse_args()

def main():
	args = parse_args()
	if args.lower is None: args.lower = 0
	if args.upper is None: args.upper = 10000

	min_size = str(args.lower)
	max_size = str(args.upper)
	sv_info_file = args.pop_vcf_file

	df = pd.read_csv(sv_info_file, sep='\t', lineterminator='\n')
	df.columns = ['chrom','start','end','length','seq']

	df['length'].value_counts()
	df = df[df['length']!='.']
	df.dropna(subset = ['length'], inplace = True)

	df['length'] = df['length'].astype(int)
	df['length'].value_counts()
	t_rows = df.query('length >= ' + min_size)
	t_rows = t_rows.query('length <=  ' + max_size)
	t = t_rows['length'].value_counts()

	### Use this to print the top N frequencies (helps when trying to narrow down a range)
	#print(t[0:20])
	#mean = t.sum()
	#print(mean)
	###

	plt.hist(t_rows['length'], density=False, bins=len(t))
	#plt.yscale('log')
	plt.ylabel('log(Frequency)')
	plt.xlabel('Insertion Size (bp)')
	plt.show()

if __name__ == '__main__':
	main()