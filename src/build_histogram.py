import pandas as pd
import matplotlib.pyplot as plt

def main():
	lower = 0 if snakemake.params.low is None else snakemake.params.low
	upper = 10000 if snakemake.params.high is None else snakemake.params.high

	min_size = str(lower)
	max_size = str(upper)
	sv_info_file = snakemake.input.global_vcf_file

	# Read file and handle varying column counts
	df = pd.read_csv(sv_info_file, sep='\t', lineterminator='\n', header=None)
	
	# Assign column names based on actual number of columns
	if df.shape[1] >= 5:
		df.columns = ['chrom','start','end','length','seq'] + [f'extra_{i}' for i in range(5, df.shape[1])]
	elif df.shape[1] == 4:
		df.columns = ['chrom','start','end','length']
		df['seq'] = ''  # Add empty sequence column
	else:
		print(f"Error: Expected at least 4 columns, got {df.shape[1]}")
		# Create empty plot
		fig, ax = plt.subplots(figsize=(8, 6))
		ax.text(0.5, 0.5, f'Invalid data format: {df.shape[1]} columns', 
				ha='center', va='center', transform=ax.transAxes)
		ax.set_title('Error: Invalid Input Data')
		plt.savefig(snakemake.output.plot)
		plt.close()
		return

	df = df[df['length']!='.']
	df.dropna(subset = ['length'], inplace = True)

	# Convert to numeric, handling invalid values
	df['length'] = pd.to_numeric(df['length'], errors='coerce')
	df = df.dropna(subset=['length'])
	
	if len(df) == 0:
		print("No valid length data found")
		# Create empty plot
		fig, ax = plt.subplots(figsize=(8, 6))
		ax.text(0.5, 0.5, 'No valid insertion length data', 
				ha='center', va='center', transform=ax.transAxes)
		ax.set_title('No Data Available')
		plt.savefig(snakemake.output.plot)
		plt.close()
		return

	df['length'] = df['length'].astype(int)
	t_rows = df.query('length >= ' + min_size)
	t_rows = t_rows.query('length <=  ' + max_size)
	
	if len(t_rows) == 0:
		print(f"No insertions found in size range {min_size}-{max_size}")
		# Create empty plot
		fig, ax = plt.subplots(figsize=(8, 6))
		ax.text(0.5, 0.5, f'No insertions in range {min_size}-{max_size} bp', 
				ha='center', va='center', transform=ax.transAxes)
		ax.set_title('No Data in Size Range')
		plt.savefig(snakemake.output.plot)
		plt.close()
		return
		
	t = t_rows['length'].value_counts()

	### Use this to print the top N frequencies (helps when trying to narrow down a range)
	#print(t[0:20])
	#mean = t.sum()
	#print(mean)
	### 

	plt.hist(t_rows['length'], density=False, bins=min(len(t), 100))  # Cap bins at 100
	plt.ylabel('Frequency')
	plt.xlabel('Insertion Size (bp)')
	plt.title(f'Insertion Size Distribution ({len(t_rows):,} insertions)')

	plt.savefig(snakemake.output.plot)
	plt.close()

if __name__ == '__main__':
	main()