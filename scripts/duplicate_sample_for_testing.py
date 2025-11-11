#!/usr/bin/env python3
"""
Script to duplicate sample files for testing parallelization.
Creates a second sample by copying all files from the first sample.
"""
import os
import shutil
import sys
import argparse
from pathlib import Path

def duplicate_sample_files(species_dir, original_sample, new_sample, verbose=True):
    """
    Duplicate all files for a sample (BAM, BAI, etc.) with a new sample name.
    """
    species_path = Path(species_dir)
    
    if not species_path.exists():
        print(f"Error: Species directory {species_dir} does not exist")
        return False
    
    files_copied = 0
    
    # Look for all files that start with the original sample name
    for file_path in species_path.glob(f"{original_sample}*"):
        if file_path.is_file():
            # Create new filename by replacing original sample name with new sample name
            new_filename = file_path.name.replace(original_sample, new_sample)
            new_file_path = species_path / new_filename
            
            if verbose:
                print(f"Copying: {file_path.name} -> {new_filename}")
            
            try:
                shutil.copy2(file_path, new_file_path)
                files_copied += 1
            except Exception as e:
                print(f"Error copying {file_path}: {e}")
                return False
    
    if verbose:
        print(f"Successfully copied {files_copied} files for sample {new_sample}")
    
    return files_copied > 0

def update_sample_list(sample_file, new_sample, verbose=True):
    """
    Add the new sample to the sample list file.
    """
    try:
        # Read existing samples
        with open(sample_file, 'r') as f:
            existing_samples = [line.strip() for line in f if line.strip()]
        
        # Check if new sample already exists
        if new_sample in existing_samples:
            if verbose:
                print(f"Sample {new_sample} already exists in {sample_file}")
            return True
        
        # Add new sample
        with open(sample_file, 'a') as f:
            f.write(f"{new_sample}\n")
        
        if verbose:
            print(f"Added {new_sample} to {sample_file}")
        
        return True
        
    except Exception as e:
        print(f"Error updating sample list: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description='Duplicate sample files for testing parallelization')
    parser.add_argument('--species-dir', required=True, help='Species directory (e.g., data/eFus)')
    parser.add_argument('--original-sample', required=True, help='Original sample name to duplicate')
    parser.add_argument('--new-sample', required=True, help='New sample name')
    parser.add_argument('--sample-file', help='Sample list file to update (optional)')
    parser.add_argument('--quiet', action='store_true', help='Suppress verbose output')
    
    args = parser.parse_args()
    
    verbose = not args.quiet
    
    if verbose:
        print("="*60)
        print("TEPEAK: DUPLICATE SAMPLE FOR TESTING")
        print("="*60)
        print(f"Species directory: {args.species_dir}")
        print(f"Original sample: {args.original_sample}")
        print(f"New sample: {args.new_sample}")
        print()
    
    # Duplicate the sample files
    success = duplicate_sample_files(
        args.species_dir, 
        args.original_sample, 
        args.new_sample, 
        verbose
    )
    
    if not success:
        print("Failed to duplicate sample files")
        sys.exit(1)
    
    # Update sample list if provided
    if args.sample_file:
        if verbose:
            print()
        
        success = update_sample_list(args.sample_file, args.new_sample, verbose)
        
        if not success:
            print("Failed to update sample list")
            sys.exit(1)
    
    if verbose:
        print("\nDuplication complete! You can now test parallelization with multiple samples.")
        print(f"To test, delete VCF outputs and rerun pipeline:")
        print(f"  rm -rf output/*/*/out.pass.vcf.gz")
        print(f"  snakemake --configfile your_config.yaml --cores 64")

if __name__ == '__main__':
    main()
