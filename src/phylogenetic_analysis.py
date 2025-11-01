"""
Phylogenetic analysis of representative TE sequences
Uses MAFFT for alignment and FastTree for phylogeny reconstruction
"""
import pandas as pd
import subprocess
import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import tempfile

def extract_representative_sequences(families_csv, min_sequences=3, max_sequences=50):
    """Extract representative sequences from TEPEAK families for phylogenetic analysis."""
    try:
        df = pd.read_csv(families_csv)
        
        # Filter for annotated families with sufficient data
        annotated = df[df['dfam_status'] == 'annotated'].copy()
        
        if len(annotated) < min_sequences:
            print(f"Warning: Only {len(annotated)} annotated families found, minimum {min_sequences} required")
            # Fall back to all families if we don't have enough annotated ones
            annotated = df[df['representative_seq'].str.len() >= 50].copy()
        
        if len(annotated) < min_sequences:
            print(f"Error: Only {len(annotated)} families with valid sequences, cannot build phylogeny")
            return []
        
        # Sort by estimated total count and take top sequences
        annotated = annotated.sort_values('estimated_total', ascending=False).head(max_sequences)
        
        print(f"Selected {len(annotated)} representative sequences for phylogenetic analysis")
        
        # Create sequence records
        sequences = []
        for idx, row in annotated.iterrows():
            # Create informative sequence ID
            seq_id = f"Fam{row['size_cluster_idx']}_{row['seq_family_idx']}"
            
            # Add annotation info if available
            if pd.notna(row['dfam_name']) and row['dfam_name']:
                seq_id += f"_{row['dfam_name'].replace(' ', '_')}"
            
            # Add size info
            seq_id += f"_{row['rep_length']}bp"
            
            # Truncate ID if too long
            if len(seq_id) > 50:
                seq_id = seq_id[:47] + "..."
            
            sequences.append(SeqRecord(
                Seq(row['representative_seq']),
                id=seq_id.replace('/', '_').replace(':', '_'),  # Clean ID for tree software
                description=f"size_range={row['size_range']} count={row['estimated_total']}"
            ))
        
        return sequences
        
    except Exception as e:
        print(f"Error extracting sequences: {e}")
        return []

def run_mafft_alignment(sequences, output_file, threads=4):
    """Run MAFFT multiple sequence alignment."""
    try:
        # Write sequences to temporary FASTA file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_fasta:
            SeqIO.write(sequences, temp_fasta, 'fasta')
            temp_fasta_path = temp_fasta.name
        
        # Run MAFFT
        cmd = [
            'mafft',
            '--auto',  # Automatic alignment strategy
            '--thread', str(threads),
            '--quiet',  # Reduce output
            temp_fasta_path
        ]
        
        print(f"Running MAFFT alignment with {len(sequences)} sequences...")
        
        with open(output_file, 'w') as outf:
            result = subprocess.run(cmd, stdout=outf, stderr=subprocess.PIPE, text=True)
        
        # Clean up temp file
        os.unlink(temp_fasta_path)
        
        if result.returncode != 0:
            print(f"MAFFT error: {result.stderr}")
            return False
        
        print(f"Alignment saved to: {output_file}")
        return True
        
    except FileNotFoundError:
        print("Error: MAFFT not found. Please install MAFFT: conda install -c bioconda mafft")
        return False
    except Exception as e:
        print(f"Error running MAFFT: {e}")
        return False

def run_fasttree_phylogeny(alignment_file, output_file, model='nt'):
    """Run FastTree phylogenetic reconstruction."""
    try:
        cmd = [
            'FastTree',
            '-nt' if model == 'nt' else '-lg',  # Nucleotide or protein model
            '-quiet',  # Reduce output
            alignment_file
        ]
        
        print("Running FastTree phylogenetic reconstruction...")
        
        with open(output_file, 'w') as outf:
            result = subprocess.run(cmd, stdout=outf, stderr=subprocess.PIPE, text=True)
        
        if result.returncode != 0:
            print(f"FastTree error: {result.stderr}")
            return False
        
        print(f"Phylogenetic tree saved to: {output_file}")
        return True
        
    except FileNotFoundError:
        print("Error: FastTree not found. Please install FastTree: conda install -c bioconda fasttree")
        return False
    except Exception as e:
        print(f"Error running FastTree: {e}")
        return False

def create_analysis_summary(families_csv, tree_file, summary_file):
    """Create a summary of the phylogenetic analysis."""
    try:
        df = pd.read_csv(families_csv)
        
        with open(summary_file, 'w') as f:
            f.write("=== TEPEAK PHYLOGENETIC ANALYSIS SUMMARY ===\n\n")
            
            # Basic statistics
            f.write(f"Total TE families analyzed: {len(df)}\n")
            annotated = df[df['dfam_status'] == 'annotated']
            f.write(f"Annotated families: {len(annotated)}\n")
            f.write(f"Unannotated families: {len(df) - len(annotated)}\n\n")
            
            # Size distribution
            f.write("Size distribution of families:\n")
            f.write(f"  Mean size: {df['rep_length'].mean():.1f} bp\n")
            f.write(f"  Size range: {df['rep_length'].min()}-{df['rep_length'].max()} bp\n\n")
            
            # TE classes found
            if len(annotated) > 0:
                f.write("TE classes identified:\n")
                class_counts = annotated['dfam_class'].value_counts()
                for te_class, count in class_counts.head(10).items():
                    if pd.notna(te_class) and te_class:
                        f.write(f"  {te_class}: {count} families\n")
                f.write("\n")
            
            # Tree file info
            if os.path.exists(tree_file):
                f.write(f"Phylogenetic tree: {os.path.basename(tree_file)}\n")
                f.write("Tree format: Newick\n")
                f.write("Method: FastTree (approximate maximum-likelihood)\n")
                f.write("Model: Nucleotide (Jukes-Cantor + CAT)\n\n")
            
            f.write("INTERPRETATION NOTES:\n")
            f.write("- Branch lengths represent evolutionary distance\n")
            f.write("- Closely related sequences cluster together\n")
            f.write("- Long branches may indicate divergent or novel TEs\n")
            f.write("- Use tree viewers like FigTree or iTOL for visualization\n")
        
        print(f"Analysis summary saved to: {summary_file}")
        
    except Exception as e:
        print(f"Error creating summary: {e}")

def main():
    """Main phylogenetic analysis function for Snakemake."""
    # Get inputs from Snakemake
    families_csv = snakemake.input.families_csv
    output_dir = snakemake.params.output_dir
    species = snakemake.params.species
    
    # Get phylogeny parameters from config
    config = snakemake.config
    threads = config.get('threads', 4)
    min_sequences = config.get('phylo_min_sequences', 3)
    max_sequences = config.get('phylo_max_sequences', 50)
    
    # Output files
    alignment_file = os.path.join(output_dir, f'{species}_te_alignment.fasta')
    tree_file = os.path.join(output_dir, f'{species}_te_phylogeny.newick')
    summary_file = os.path.join(output_dir, f'{species}_phylogeny_summary.txt')
    
    print("="*60)
    print("TEPEAK: PHYLOGENETIC ANALYSIS")
    print("="*60)
    
    # Step 1: Extract representative sequences
    sequences = extract_representative_sequences(
        families_csv, 
        min_sequences=min_sequences,
        max_sequences=max_sequences
    )
    
    if len(sequences) < min_sequences:
        print(f"Insufficient sequences ({len(sequences)}) for phylogenetic analysis")
        # Create empty output files
        open(alignment_file, 'w').close()
        open(tree_file, 'w').close()
        with open(summary_file, 'w') as f:
            f.write("Phylogenetic analysis skipped: insufficient sequences\n")
        return
    
    # Step 2: Multiple sequence alignment
    alignment_success = run_mafft_alignment(sequences, alignment_file, threads)
    
    if not alignment_success:
        print("Alignment failed, skipping phylogenetic reconstruction")
        open(tree_file, 'w').close()
        with open(summary_file, 'w') as f:
            f.write("Phylogenetic analysis failed: alignment error\n")
        return
    
    # Step 3: Phylogenetic reconstruction
    tree_success = run_fasttree_phylogeny(alignment_file, tree_file, model='nt')
    
    if not tree_success:
        print("Phylogenetic reconstruction failed")
        with open(summary_file, 'w') as f:
            f.write("Phylogenetic analysis failed: tree reconstruction error\n")
        return
    
    # Step 4: Create summary
    create_analysis_summary(families_csv, tree_file, summary_file)
    
    print("\nPhylogenetic analysis complete!")
    print(f"Alignment: {alignment_file}")
    print(f"Tree: {tree_file}")
    print(f"Summary: {summary_file}")

if __name__ == '__main__':
    main()
