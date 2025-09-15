import os, json
from pathlib import Path

# noqa: F821

def main():
    species = snakemake.params.species
    species_dir = Path(snakemake.params.species_dir)
    zip_path = Path(snakemake.input.ref)

    species_dir.mkdir(parents=True, exist_ok=True)

    # Unzip directly into species_dir (creates species_dir/ncbi_dataset/…)
    catalog_probe = species_dir / "ncbi_dataset" / "data" / "dataset_catalog.json"
    if not catalog_probe.exists():
        os.system(f'unzip -o "{zip_path}" -d "{species_dir}" > /dev/null')
    else:
        print("Dataset already extracted, proceeding...")

    # Find dataset_catalog.json under species_dir/ncbi_dataset
    candidates = list((species_dir / "ncbi_dataset").glob("**/data/dataset_catalog.json"))
    if not candidates:
        raise FileNotFoundError(f"No dataset_catalog.json under {species_dir/'ncbi_dataset'}")
    catalog_path = candidates[0]

    with open(catalog_path) as f:
        parsed = json.load(f)

    fa_rel = parsed["assemblies"][-1]["files"][0]["filePath"]

    # Resolve FASTA under species_dir/ncbi_dataset/data/**
    fa_root = species_dir / "ncbi_dataset" / "data"
    fa_candidates = [p for p in fa_root.glob("**/*") if str(p).endswith(fa_rel)]
    fa_path = fa_candidates[0] if fa_candidates else (fa_root / fa_rel)

    target_fa = species_dir / f"{species}.fa"
    os.system(f'cp "{fa_path}" "{target_fa}"')
    print(f"File has been moved and renamed to {species}.fa")

    os.system(f'cd "{species_dir}" && samtools faidx "{species}.fa" && bwa index -p "{species}.fa" "{species}.fa"')

if __name__ == "__main__":
    main()
