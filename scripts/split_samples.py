#!/usr/bin/env python3
"""
TEPEAK Sample Splitter
Split sample files into smaller chunks for batch processing
"""

import argparse
import os
import math

def split_samples(input_file, output_prefix, chunk_size=None, num_chunks=None):
    """Split samples into chunks."""
    
    # Read all samples
    with open(input_file, 'r') as f:
        samples = [line.strip() for line in f if line.strip()]
    
    total_samples = len(samples)
    print(f"Total samples: {total_samples}")
    
    # Determine chunking strategy
    if chunk_size:
        num_chunks = math.ceil(total_samples / chunk_size)
        actual_chunk_size = chunk_size
    elif num_chunks:
        actual_chunk_size = math.ceil(total_samples / num_chunks)
    else:
        # Default: aim for ~10 samples per chunk
        actual_chunk_size = 10
        num_chunks = math.ceil(total_samples / actual_chunk_size)
    
    print(f"Creating {num_chunks} chunks with ~{actual_chunk_size} samples each")
    
    # Create output directory
    output_dir = os.path.dirname(output_prefix) or '.'
    os.makedirs(output_dir, exist_ok=True)
    
    chunk_files = []
    
    # Split samples into chunks
    for i in range(num_chunks):
        start_idx = i * actual_chunk_size
        end_idx = min((i + 1) * actual_chunk_size, total_samples)
        
        if start_idx >= total_samples:
            break
            
        chunk_samples = samples[start_idx:end_idx]
        chunk_file = f"{output_prefix}_chunk{i+1:02d}.txt"
        
        with open(chunk_file, 'w') as f:
            for sample in chunk_samples:
                f.write(f"{sample}\n")
        
        chunk_files.append(chunk_file)
        print(f"  Chunk {i+1}: {len(chunk_samples)} samples → {chunk_file}")
    
    return chunk_files

def main():
    parser = argparse.ArgumentParser(description='Split TEPEAK sample files into smaller chunks')
    parser.add_argument('input_file', help='Input sample file to split')
    parser.add_argument('output_prefix', help='Output prefix (e.g., "chunks/samples")')
    
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--chunk-size', type=int, 
                      help='Number of samples per chunk')
    group.add_argument('--num-chunks', type=int, 
                      help='Number of chunks to create')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input_file):
        print(f"Error: Input file {args.input_file} not found")
        return 1
    
    chunk_files = split_samples(
        args.input_file, 
        args.output_prefix,
        chunk_size=args.chunk_size,
        num_chunks=args.num_chunks
    )
    
    print(f"\nSample splitting complete!")
    print(f"Created {len(chunk_files)} chunk files")
    
    print(f"\nTo run Stage 1 on each chunk:")
    for i, chunk_file in enumerate(chunk_files, 1):
        print(f"  # Chunk {i}")
        print(f"  cp {chunk_file} temp_samples.txt")
        print(f"  snakemake -s Snakefile_stage1.py --configfile config.yaml --cores 64")
        print()
    
    print("After all Stage 1 jobs complete, run Stage 2:")
    print("  snakemake -s Snakefile_stage2.py --configfile config.yaml --cores 32")

if __name__ == '__main__':
    main()
