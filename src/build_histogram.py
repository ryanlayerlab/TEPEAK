import pandas as pd
import matplotlib.pyplot as plt

def main():
	lower = 0 if snakemake.params.low is None else snakemake.params.low
	upper = 10000 if snakemake.params.high is None else snakemake.params.high

	min_size = str(lower)
	max_size = str(upper)
	sv_info_file = snakemake.input.glpbal_vcf_file

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

	plt.savefig(snakemake.output.plot)

if __name__ == '__main__':
	main()