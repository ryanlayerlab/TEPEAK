import os, json

def main():
    os.system(f'echo "n" | unzip {snakemake.input.ref}')
    species = snakemake.params.species
    species_dir = snakemake.params.species_dir

    with open('ncbi_dataset/data/dataset_catalog.json') as dataset_catalog: 
        parsed_dataset = json.load(dataset_catalog)
    
    fa_filepath = f"ncbi_dataset/data/{parsed_dataset['assemblies'][-1]['files'][0]['filePath']}"
    os.system(
        f"""
        cp {fa_filepath} {species_dir}/{species}.fa
        echo "File has been moved and renamed to {species}.fa"
        rm -r ncbi_dataset
        cd {species_dir}
        samtools faidx {species}.fa
        bwa index -p {species} {species}.fa
        """
    )

if __name__ == "__main__":
    main()