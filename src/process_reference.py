# src/process_reference.py
import os, sys, zipfile, gzip, shutil, tempfile, subprocess
from pathlib import Path


ref_zip = str(snakemake.input.ref) 
species_dir = str(snakemake.params.species_dir).rstrip("/")
species = str(snakemake.params.species)

out_fa = f"{species_dir}/{species}.fa"

def run(cmd):
    subprocess.run(cmd, check=True)

def find_first_fasta(root: Path):
    exts = {".fa", ".fna", ".fasta"}
    gz_exts = {e + ".gz" for e in exts}
    # Prefer uncompressed, then gz
    cands_plain, cands_gz = [], []
    for p in root.rglob("*"):
        if not p.is_file():
            continue
        suf = "".join(p.suffixes[-2:]) if len(p.suffixes) >= 2 else p.suffix
        if p.suffix.lower() in exts:
            cands_plain.append(p)
        elif suf.lower() in gz_exts:
            cands_gz.append(p)
    if cands_plain:
        return cands_plain[0], False
    if cands_gz:
        return cands_gz[0], True
    return None, False

def main():
    os.makedirs(species_dir, exist_ok=True)

    # Check if reference and indexes already exist
    if all(Path(p).exists() and Path(p).stat().st_size > 0 for p in snakemake.output.ref):
        print("Reference and indexes already present; skipping rebuild.", file=sys.stderr)
        return
    
    # Use lock file to prevent multiple jobs from building reference simultaneously
    lock_file = Path(species_dir) / ".reference_building.lock"
    if lock_file.exists():
        print("Another process is building the reference. Waiting...", file=sys.stderr)
        import time
        while lock_file.exists() and time.time() - lock_file.stat().st_mtime < 3600:  # 1 hour timeout
            time.sleep(10)
        # Check again if reference is now complete
        if all(Path(p).exists() and Path(p).stat().st_size > 0 for p in snakemake.output.ref):
            print("Reference completed by another process.", file=sys.stderr)
            return
    
    # Create lock file
    try:
        lock_file.touch()
        print(f"Created lock file: {lock_file}", file=sys.stderr)

        # Check if ref_zip is actually a zip file or a plain FASTA
        ref_path = Path(ref_zip)
        
        if ref_path.suffix.lower() == '.zip':
            # Handle ZIP file case
            with tempfile.TemporaryDirectory() as tmpd:
                tmpd = Path(tmpd)
                with zipfile.ZipFile(ref_zip) as zf:
                    zf.extractall(tmpd)

                fasta_path, is_gz = find_first_fasta(tmpd)
                if not fasta_path:
                    print(f"[process_reference] ERROR: No FASTA-like file (.fa/.fna/.fasta[.gz]) found inside {ref_zip}", file=sys.stderr)
                    sys.exit(1)

                # Normalize to species.fa (decompress if needed)
                if is_gz:
                    with gzip.open(fasta_path, "rb") as inp, open(out_fa, "wb") as outp:
                        shutil.copyfileobj(inp, outp)
                else:
                    shutil.copyfile(fasta_path, out_fa)
        else:
            # Handle plain FASTA file case
            if ref_path.suffix.lower() in ['.gz']:
                # Compressed FASTA
                with gzip.open(ref_zip, "rb") as inp, open(out_fa, "wb") as outp:
                    shutil.copyfileobj(inp, outp)
            else:
                # Plain FASTA
                shutil.copyfile(ref_zip, out_fa)

        # Build indexes
        run(["samtools", "faidx", out_fa])
        run(["bwa", "index", "-p", out_fa, out_fa])
        
        # Remove lock file
        if lock_file.exists():
            lock_file.unlink()
            print(f"Removed lock file: {lock_file}", file=sys.stderr)
    except Exception as e:
        # Always try to clean up lock file on error
        if 'lock_file' in locals() and lock_file.exists():
            lock_file.unlink()
        raise e

if __name__ == "__main__":
    main()
