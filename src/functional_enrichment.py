"""
Functional enrichment analysis of genes intersecting with TE insertions
Uses g:Profiler API for GO and pathway enrichment
"""
import pandas as pd
import requests
import json
import os
import sys
from collections import defaultdict

def extract_gene_names_from_annotation(annotation_file):
    """Extract unique gene names from the gene annotation file."""
    try:
        df = pd.read_csv(annotation_file, sep='\t')
        
        # Handle different possible column names
        gene_col = None
        for col in ['gene_name', 'gene_id', 'gene', 'Gene', 'GENE']:
            if col in df.columns:
                gene_col = col
                break
        
        if gene_col is None:
            print("Warning: No gene name column found in annotation file")
            return []
        
        # Filter out 'None' and empty values, get unique genes
        genes = df[gene_col].dropna()
        genes = genes[genes != 'None']
        genes = genes[genes != '']
        unique_genes = list(set(genes.astype(str)))
        
        print(f"Extracted {len(unique_genes)} unique genes for enrichment analysis")
        return unique_genes
        
    except Exception as e:
        print(f"Error reading annotation file: {e}")
        return []

def run_gprofiler_enrichment(gene_list, organism='hsapiens', sources=None, max_pvalue=0.05):
    """
    Run functional enrichment analysis using g:Profiler API.
    
    Parameters:
    -----------
    gene_list : list
        List of gene names/IDs
    organism : str
        Organism code (e.g., 'hsapiens', 'mmusculus', 'dmelanogaster')
    sources : list
        Data sources to query (default: GO, KEGG, Reactome)
    max_pvalue : float
        Maximum p-value threshold
    """
    if sources is None:
        sources = ['GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC']
    
    if len(gene_list) == 0:
        print("No genes provided for enrichment analysis")
        return pd.DataFrame()
    
    print(f"Running g:Profiler enrichment analysis for {len(gene_list)} genes")
    print(f"Organism: {organism}")
    print(f"Sources: {', '.join(sources)}")
    
    # g:Profiler API endpoint
    url = 'https://biit.cs.ut.ee/gprofiler/api/gost/profile/'
    
    # Prepare request data
    data = {
        'organism': organism,
        'query': gene_list,
        'sources': sources,
        'user_threshold': max_pvalue,
        'significance_threshold_method': 'fdr',
        'domain_scope': 'annotated',
        'measure_underrepresentation': False,
        'evcodes': True,
        'combined': True,
        'ordered': False
    }
    
    try:
        # Make request
        response = requests.post(url, json=data, timeout=120)
        response.raise_for_status()
        
        result = response.json()
        
        if 'result' not in result or not result['result']:
            print("No significant enrichment results found")
            return pd.DataFrame()
        
        # Parse results
        enrichment_results = []
        for term in result['result']:
            enrichment_results.append({
                'term_id': term.get('native', ''),
                'term_name': term.get('name', ''),
                'source': term.get('source', ''),
                'p_value': term.get('p_value', 1.0),
                'adjusted_p_value': term.get('p_value', 1.0),  # g:Profiler returns adjusted p-values
                'term_size': term.get('term_size', 0),
                'query_size': term.get('query_size', 0),
                'intersection_size': term.get('intersection_size', 0),
                'intersection_genes': ','.join(term.get('intersections', [[]])[0]),
                'effective_domain_size': term.get('effective_domain_size', 0),
                'precision': term.get('precision', 0.0),
                'recall': term.get('recall', 0.0)
            })
        
        df_results = pd.DataFrame(enrichment_results)
        df_results = df_results.sort_values('adjusted_p_value')
        
        print(f"Found {len(df_results)} significant enrichment terms")
        return df_results
        
    except requests.exceptions.RequestException as e:
        print(f"Error connecting to g:Profiler API: {e}")
        return pd.DataFrame()
    except Exception as e:
        print(f"Error processing enrichment results: {e}")
        return pd.DataFrame()

def format_enrichment_summary(df_results, max_terms=20):
    """Create a summary of top enrichment results."""
    if df_results.empty:
        return "No significant enrichment terms found.\n"
    
    summary = []
    summary.append("=== FUNCTIONAL ENRICHMENT ANALYSIS SUMMARY ===\n")
    summary.append(f"Total significant terms: {len(df_results)}")
    summary.append(f"Showing top {min(max_terms, len(df_results))} results:\n")
    
    # Group by source
    sources = df_results['source'].unique()
    for source in sorted(sources):
        source_results = df_results[df_results['source'] == source].head(10)
        if len(source_results) > 0:
            summary.append(f"\n--- {source} ---")
            for _, row in source_results.iterrows():
                summary.append(f"{row['term_name']}")
                summary.append(f"  ID: {row['term_id']}")
                summary.append(f"  P-value: {row['adjusted_p_value']:.2e}")
                summary.append(f"  Genes: {row['intersection_size']}/{row['query_size']}")
                if len(row['intersection_genes']) < 200:  # Don't show if too long
                    summary.append(f"  Intersecting genes: {row['intersection_genes']}")
                summary.append("")
    
    return "\n".join(summary)

def detect_organism_from_config(species_name):
    """
    Attempt to map species name to g:Profiler organism code.
    Returns best guess or 'hsapiens' as default.
    """
    species_lower = species_name.lower()
    
    organism_mapping = {
        'human': 'hsapiens',
        'homo': 'hsapiens',
        'sapiens': 'hsapiens',
        'mouse': 'mmusculus',
        'mus': 'mmusculus',
        'musculus': 'mmusculus',
        'rat': 'rnorvegicus',
        'rattus': 'rnorvegicus',
        'fly': 'dmelanogaster',
        'drosophila': 'dmelanogaster',
        'melanogaster': 'dmelanogaster',
        'worm': 'celegans',
        'caenorhabditis': 'celegans',
        'elegans': 'celegans',
        'yeast': 'scerevisiae',
        'saccharomyces': 'scerevisiae',
        'cerevisiae': 'scerevisiae',
        'zebrafish': 'drerio',
        'danio': 'drerio',
        'arabidopsis': 'athaliana',
        'thaliana': 'athaliana',
        'rice': 'osativa',
        'oryza': 'osativa'
    }
    
    for key, organism in organism_mapping.items():
        if key in species_lower:
            print(f"Detected organism: {organism} (from species: {species_name})")
            return organism
    
    print(f"Could not detect organism from species '{species_name}', using human (hsapiens) as default")
    print("Set 'enrichment_organism' in config to override")
    return 'hsapiens'

def main():
    """Main enrichment analysis function for Snakemake."""
    # Get inputs from Snakemake
    annotation_file = snakemake.input.annotation_file
    output_dir = snakemake.params.output_dir
    species = snakemake.params.species
    
    # Get enrichment parameters from config
    config = snakemake.config
    organism = config.get('enrichment_organism', None)
    if organism is None:
        organism = detect_organism_from_config(species)
    
    sources = config.get('enrichment_sources', ['GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC'])
    max_pvalue = config.get('enrichment_max_pvalue', 0.05)
    max_terms = config.get('enrichment_max_terms', 50)
    
    # Output files
    results_file = os.path.join(output_dir, f'{species}_enrichment_results.csv')
    summary_file = os.path.join(output_dir, f'{species}_enrichment_summary.txt')
    
    print("="*60)
    print("FUNCTIONAL ENRICHMENT ANALYSIS")
    print("="*60)
    
    # Extract genes from annotation file
    genes = extract_gene_names_from_annotation(annotation_file)
    
    if len(genes) == 0:
        print("No genes found for enrichment analysis")
        # Create empty output files
        pd.DataFrame().to_csv(results_file, index=False)
        with open(summary_file, 'w') as f:
            f.write("No genes found for enrichment analysis.\n")
        return
    
    if len(genes) < 5:
        print(f"Too few genes ({len(genes)}) for meaningful enrichment analysis")
        # Create empty output files with explanation
        pd.DataFrame().to_csv(results_file, index=False)
        with open(summary_file, 'w') as f:
            f.write(f"Too few genes ({len(genes)}) for meaningful enrichment analysis.\n")
            f.write("Minimum 5 genes recommended.\n")
        return
    
    # Run enrichment analysis
    results_df = run_gprofiler_enrichment(
        gene_list=genes,
        organism=organism,
        sources=sources,
        max_pvalue=max_pvalue
    )
    
    # Save detailed results
    if not results_df.empty:
        # Limit to top results
        results_df_top = results_df.head(max_terms)
        results_df_top.to_csv(results_file, index=False)
        print(f"Detailed results saved to: {results_file}")
    else:
        pd.DataFrame().to_csv(results_file, index=False)
    
    # Create and save summary
    summary = format_enrichment_summary(results_df, max_terms=20)
    with open(summary_file, 'w') as f:
        f.write(summary)
    
    print(f"Summary saved to: {summary_file}")
    print("\n" + summary)

if __name__ == '__main__':
    main()
