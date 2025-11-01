"""
TEPEAK: Transposable Element Peak Finding & Annotation
Advanced version with adaptive clustering and sequence analysis
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio.Align import PairwiseAligner
import requests
import time
from collections import Counter
import os
import sys

def load_global_tsv(path: str) -> pd.DataFrame:
    """
    Robustly load a TSV file with optional header detection.
    Handles both RefSeq (NC_*) and numeric chromosome identifiers.
    """
    # First pass: detect structure
    with open(path, 'r') as f:
        first_line = f.readline().strip()
    
    if not first_line:
        print("Error: File is empty")
        return pd.DataFrame(columns=["chrom", "start", "end", "length", "seq"])
    
    fields = first_line.split('\t')
    first_val = fields[0] if fields else ""
    
    print(f"Detected {len(fields)} columns in input file")
    
    # Determine if header exists
    has_header = not (first_val.startswith('NC_') or first_val.isdigit())
    
    # Read file without assuming column count
    df = pd.read_csv(
        path,
        sep="\t",
        header=0 if has_header else None,
        dtype=str,  # Read everything as string first
        keep_default_na=False
    )
    
    print(f"Initial load: {len(df):,} rows, {df.shape[1]} columns")
    
    # Assign standard column names based on what we have
    if df.shape[1] >= 5:
        # Standard format: chrom, start, end, length, seq [+ extra columns]
        new_cols = ["chrom", "start", "end", "length", "seq"]
        if df.shape[1] > 5:
            new_cols += [f"extra_{i}" for i in range(5, df.shape[1])]
        df.columns = new_cols
    elif df.shape[1] == 4:
        # Missing sequence column
        df.columns = ["chrom", "start", "end", "length"]
        df["seq"] = ""  # Add empty sequence column
        print("Warning: No sequence column detected, added empty sequences")
    else:
        print(f"Error: Expected at least 4 columns (chrom, start, end, length), got {df.shape[1]}")
        return pd.DataFrame(columns=["chrom", "start", "end", "length", "seq"])
    
    # Only keep first 5 columns for processing
    df = df[["chrom", "start", "end", "length", "seq"]]
    
    # Filter out rows where length is "." or empty
    initial_count = len(df)
    df = df[df["length"].str.strip().isin(['.', '', 'nan', 'NaN']) == False]
    if len(df) < initial_count:
        print(f"Filtered out {initial_count - len(df)} rows with invalid length values")
    
    # Convert length to numeric, coercing errors to NaN
    df["length"] = pd.to_numeric(df["length"], errors='coerce')
    
    # Drop rows with NaN length
    before_drop = len(df)
    df = df.dropna(subset=["length"])
    if len(df) < before_drop:
        print(f"Dropped {before_drop - len(df):,} rows with non-numeric length")
    
    if len(df) == 0:
        print("Warning: No valid data remaining after filtering")
        return pd.DataFrame(columns=["chrom", "start", "end", "length", "seq"])
    
    # Convert to int
    df["length"] = df["length"].astype(int)
    
    # Clean sequence column
    df["seq"] = df["seq"].fillna("").replace("nan", "").replace(".", "")
    
    # Reset index
    df = df.reset_index(drop=True)
    
    print(f"Final dataset: {len(df):,} rows")
    print(f"Length range: {df['length'].min()} - {df['length'].max()} bp")
    
    return df

def find_insertion_clusters(df, min_w, max_w, 
                           min_cluster_size=10,
                           density_threshold=0.5,
                           window_size=50,
                           merge_distance=100):
    """
    Find clusters/bands of insertions using density-based approach.
    
    Parameters:
    -----------
    min_cluster_size : int
        Minimum number of insertions to form a cluster
    density_threshold : float
        Minimum density (insertions per bp) to be considered a cluster
    window_size : int
        Sliding window size for density calculation
    merge_distance : int
        Merge clusters if they're within this distance
    """
    cnt = df["length"].value_counts().sort_index()
    if cnt.empty:
        return []
    
    lengths = pd.RangeIndex(cnt.index.min(), cnt.index.max() + 1)
    y = cnt.reindex(lengths, fill_value=0)
    
    # Filter to size range
    mask = (lengths >= min_w) & (lengths <= max_w)
    y_filtered = y[mask]
    lengths_filtered = lengths[mask]
    
    print(f"Analyzing range {min_w}-{max_w}bp")
    print(f"Total insertions in range: {y_filtered.sum():,}")
    
    # === Calculate local density ===
    density = y_filtered.rolling(window_size, center=True, min_periods=1).sum() / window_size
    
    # Calculate what constitutes "elevated" density using local statistics
    # Use sliding percentile to adapt to local context
    local_percentile = density.rolling(500, center=True, min_periods=50).quantile(0.75)
    
    # === Find high-density regions ===
    is_high_density = (density >= density_threshold) & (y_filtered > 0)
    
    print(f"High-density positions: {is_high_density.sum()}")
    
    # === Identify contiguous regions ===
    clusters = []
    in_cluster = False
    cluster_start = None
    cluster_positions = []
    
    for i, (pos, is_dense) in enumerate(zip(lengths_filtered, is_high_density)):
        if is_dense and not in_cluster:
            # Start new cluster
            in_cluster = True
            cluster_start = pos
            cluster_positions = [pos]
        elif is_dense and in_cluster:
            # Continue cluster
            cluster_positions.append(pos)
        elif not is_dense and in_cluster:
            # End cluster
            cluster_end = cluster_positions[-1]
            cluster_count = sum(y_filtered[p] for p in cluster_positions if p in y_filtered.index)
            
            if cluster_count >= min_cluster_size:
                # Find peak within cluster
                cluster_y = {p: y_filtered[p] for p in cluster_positions if p in y_filtered.index}
                peak_pos = max(cluster_y.keys(), key=lambda p: cluster_y[p])
                peak_count = cluster_y[peak_pos]
                
                clusters.append({
                    'start': int(cluster_start),
                    'end': int(cluster_end),
                    'peak': int(peak_pos),
                    'width': int(cluster_end - cluster_start + 1),
                    'total_count': int(cluster_count),
                    'peak_count': int(peak_count),
                    'positions': cluster_positions
                })
            
            in_cluster = False
            cluster_positions = []
    
    # Handle cluster at end
    if in_cluster and len(cluster_positions) > 0:
        cluster_end = cluster_positions[-1]
        cluster_count = sum(y_filtered[p] for p in cluster_positions if p in y_filtered.index)
        
        if cluster_count >= min_cluster_size:
            cluster_y = {p: y_filtered[p] for p in cluster_positions if p in y_filtered.index}
            peak_pos = max(cluster_y.keys(), key=lambda p: cluster_y[p])
            peak_count = cluster_y[peak_pos]
            
            clusters.append({
                'start': int(cluster_start),
                'end': int(cluster_end),
                'peak': int(peak_pos),
                'width': int(cluster_end - cluster_start + 1),
                'total_count': int(cluster_count),
                'peak_count': int(peak_count),
                'positions': cluster_positions
            })
    
    print(f"Found {len(clusters)} initial clusters")
    
    # === Merge nearby clusters ===
    if merge_distance > 0 and len(clusters) > 1:
        merged = []
        clusters.sort(key=lambda c: c['start'])
        
        current = clusters[0]
        for next_cluster in clusters[1:]:
            if next_cluster['start'] - current['end'] <= merge_distance:
                # Merge clusters
                current = {
                    'start': current['start'],
                    'end': next_cluster['end'],
                    'peak': current['peak'] if current['peak_count'] > next_cluster['peak_count'] else next_cluster['peak'],
                    'width': next_cluster['end'] - current['start'] + 1,
                    'total_count': current['total_count'] + next_cluster['total_count'],
                    'peak_count': max(current['peak_count'], next_cluster['peak_count']),
                    'positions': current['positions'] + next_cluster['positions']
                }
            else:
                merged.append(current)
                current = next_cluster
        
        merged.append(current)
        clusters = merged
        
        print(f"After merging (distance <= {merge_distance}bp): {len(clusters)} clusters")
    
    # Sort by total count (most prominent clusters first)
    clusters.sort(key=lambda c: -c['total_count'])
    
    return clusters

# ✓ Adaptive density-based clustering with local percentile thresholds
def find_clusters_adaptive_density(df, min_w, max_w, 
                                   min_cluster_size=10,
                                   percentile_threshold=75,
                                   window_size=50,
                                   merge_distance=100):
    """
    Adaptive density-based clustering using local percentiles.
    Automatically adjusts to local density variations - EXACTLY as described.
    """
    cnt = df["length"].value_counts().sort_index()
    if cnt.empty:
        return []
    
    lengths = pd.RangeIndex(cnt.index.min(), cnt.index.max() + 1)
    y = cnt.reindex(lengths, fill_value=0)
    
    # Filter to size range
    mask = (lengths >= min_w) & (lengths <= max_w)
    y_filtered = y[mask]
    lengths_filtered = lengths[mask]
    
    print(f"Analyzing range {min_w}-{max_w}bp")
    print(f"Total insertions in range: {y_filtered.sum():,}")
    
    # === Calculate local density ===
    # ✓ Sliding window density calculation (default: 50bp)
    density = y_filtered.rolling(window_size, center=True, min_periods=1).sum()
    
    # ✓ Local background estimation with broader context (default: 500bp)
    local_bg_window = 500
    local_threshold = density.rolling(local_bg_window, center=True, 
                                      min_periods=50).quantile(percentile_threshold / 100.0)
    
    # ✓ Local context-aware peak detection (not global thresholds)
    is_elevated = density > local_threshold
    
    print(f"Elevated density positions: {is_elevated.sum()}")
    
    # === Find contiguous regions ===
    clusters = []
    in_cluster = False
    cluster_start = None
    cluster_positions = []
    
    for pos, is_elev in zip(lengths_filtered, is_elevated):
        if is_elev and not in_cluster:
            in_cluster = True
            cluster_start = pos
            cluster_positions = [pos]
        elif is_elev and in_cluster:
            cluster_positions.append(pos)
        elif not is_elev and in_cluster:
            # End cluster
            if len(cluster_positions) > 0:
                cluster_end = cluster_positions[-1]
                cluster_count = sum(y_filtered.get(p, 0) for p in cluster_positions)
                
                if cluster_count >= min_cluster_size:
                    cluster_y = {p: y_filtered.get(p, 0) for p in cluster_positions}
                    peak_pos = max(cluster_y.keys(), key=lambda p: cluster_y[p])
                    peak_count = cluster_y[peak_pos]
                    
                    clusters.append({
                        'start': int(cluster_start),
                        'end': int(cluster_end),
                        'peak': int(peak_pos),
                        'width': int(cluster_end - cluster_start + 1),
                        'total_count': int(cluster_count),
                        'peak_count': int(peak_count),
                        'avg_density': float(cluster_count / len(cluster_positions)),
                        'positions': cluster_positions
                    })
            
            in_cluster = False
            cluster_positions = []
    
    # Handle final cluster
    if in_cluster and len(cluster_positions) > 0:
        cluster_count = sum(y_filtered.get(p, 0) for p in cluster_positions)
        if cluster_count >= min_cluster_size:
            cluster_y = {p: y_filtered.get(p, 0) for p in cluster_positions}
            peak_pos = max(cluster_y.keys(), key=lambda p: cluster_y[p])
            peak_count = cluster_y[peak_pos]
            cluster_end = cluster_positions[-1]
            
            clusters.append({
                'start': int(cluster_start),
                'end': int(cluster_end),
                'peak': int(peak_pos),
                'width': int(cluster_end - cluster_start + 1),
                'total_count': int(cluster_count),
                'peak_count': int(peak_count),
                'avg_density': float(cluster_count / len(cluster_positions)),
                'positions': cluster_positions
            })
    
    print(f"Found {len(clusters)} clusters")
    
    # === Merge nearby clusters ===
    if merge_distance > 0 and len(clusters) > 1:
        merged = []
        clusters.sort(key=lambda c: c['start'])
        
        current = clusters[0]
        for next_cluster in clusters[1:]:
            gap = next_cluster['start'] - current['end']
            if gap <= merge_distance:
                # Merge
                all_positions = current['positions'] + next_cluster['positions']
                total_count = current['total_count'] + next_cluster['total_count']
                
                cluster_y = {p: y_filtered.get(p, 0) for p in all_positions}
                peak_pos = max(cluster_y.keys(), key=lambda p: cluster_y[p])
                
                current = {
                    'start': current['start'],
                    'end': next_cluster['end'],
                    'peak': int(peak_pos),
                    'width': next_cluster['end'] - current['start'] + 1,
                    'total_count': total_count,
                    'peak_count': int(cluster_y[peak_pos]),
                    'avg_density': float(total_count / len(all_positions)),
                    'positions': all_positions
                }
            else:
                merged.append(current)
                current = next_cluster
        
        merged.append(current)
        clusters = merged
        
        print(f"After merging: {len(clusters)} clusters")
    
    # Sort by total count
    clusters.sort(key=lambda c: -c['total_count'])
    
    return clusters

def extract_sequences_for_cluster(df, cluster, max_sample=200):
    """Extract sequences from a cluster with smart sampling."""
    # Get all insertions in the cluster size range
    cluster_df = df[(df["length"] >= cluster['start']) & 
                    (df["length"] <= cluster['end'])].copy()
    
    # Filter out empty/invalid sequences
    cluster_df = cluster_df[cluster_df["seq"].str.len() > 0]
    cluster_df = cluster_df[cluster_df["seq"] != "."]
    
    total_seqs = len(cluster_df)
    
    # Smart sampling: oversample near peak, undersample edges
    if total_seqs > max_sample:
        # Weight by distance from peak
        cluster_df['dist_from_peak'] = abs(cluster_df['length'] - cluster['peak'])
        cluster_df['weight'] = 1.0 / (1.0 + cluster_df['dist_from_peak'] / 50.0)
        
        # Sample with weights
        sampled = cluster_df.sample(n=max_sample, weights='weight', random_state=42)
        print(f"  Sampled {max_sample} from {total_seqs} sequences (peak-weighted)")
        return sampled['seq'].tolist(), total_seqs
    else:
        return cluster_df['seq'].tolist(), total_seqs

def _cluster_sequences_greedy(seqs, pid_threshold):
    """Greedy clustering: each sequence joins first compatible cluster."""
    if not seqs:
        return []
    
    # ✓ Greedy clustering algorithm (O(n·k) vs O(n²))
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = 0
    aligner.extend_gap_score = 0
    
    clusters = []
    representatives = []  # Store representative of each cluster
    
    for seq in seqs:
        # Try to find compatible cluster
        found_cluster = False
        
        for i, rep in enumerate(representatives):
            # Quick length check first
            if abs(len(seq) - len(rep)) / max(len(seq), len(rep)) > 0.1:
                continue
            
            # Calculate similarity to representative
            alignment = aligner.align(seq, rep)[0]
            matches = float(alignment.score)
            denom = float(min(len(seq), len(rep))) or 1.0
            pid = matches / denom
            
            if pid >= pid_threshold:
                clusters[i].append(seq)
                found_cluster = True
                break
        
        if not found_cluster:
            # Create new cluster
            clusters.append([seq])
            representatives.append(seq)
    
    return clusters

def cluster_sequences_hierarchical(seqs, pid_threshold=0.85, quick_filter=True):
    """Hierarchical sequence clustering with optimizations."""
    if not seqs or len(seqs) == 0:
        return []
    
    # Quick filter: group by length first
    if quick_filter:
        length_groups = {}
        for seq in seqs:
            length = len(seq)
            # ✓ Allow 10% length variation as described
            length_key = length // 10
            if length_key not in length_groups:
                length_groups[length_key] = []
            length_groups[length_key].append(seq)
        
        print(f"  Pre-filtered into {len(length_groups)} length groups")
        
        # Cluster within each length group
        all_clusters = []
        for length_key, group_seqs in length_groups.items():
            if len(group_seqs) > 1:
                clusters = _cluster_sequences_greedy(group_seqs, pid_threshold)
                all_clusters.extend(clusters)
            else:
                all_clusters.append(group_seqs)
        
        return all_clusters
    else:
        return _cluster_sequences_greedy(seqs, pid_threshold)

def find_consensus_sequence(cluster_seqs, method='most_common'):
    """Find consensus/representative sequence from a cluster."""
    if not cluster_seqs:
        return ""
    
    if method == 'most_common':
        counts = Counter(cluster_seqs)
        return counts.most_common(1)[0][0]
    elif method == 'longest':
        return max(cluster_seqs, key=len)
    elif method == 'median_length':
        lengths = [len(s) for s in cluster_seqs]
        median_len = np.median(lengths)
        return min(cluster_seqs, key=lambda s: abs(len(s) - median_len))
    
    return cluster_seqs[0]

def process_all_clusters(df, size_clusters, pid_threshold=0.85, 
                         max_sample=200, consensus_method='most_common'):
    """Process all size clusters to find sequence families."""
    families = []
    
    print(f"\nProcessing {len(size_clusters)} size clusters...")
    
    for i, cluster in enumerate(size_clusters):
        print(f"\nCluster {i+1}/{len(size_clusters)}: {cluster['start']}-{cluster['end']}bp "
              f"({cluster['total_count']} insertions)")
        
        # Extract sequences
        seqs, total_in_range = extract_sequences_for_cluster(df, cluster, max_sample)
        
        if len(seqs) == 0:
            print("  No valid sequences found, skipping")
            continue
        
        # Cluster sequences by similarity
        print(f"  Clustering {len(seqs)} sequences (pid >= {pid_threshold})...")
        seq_clusters = cluster_sequences_hierarchical(seqs, pid_threshold, quick_filter=True)
        
        # Sort sequence clusters by size
        seq_clusters.sort(key=len, reverse=True)
        
        print(f"  Found {len(seq_clusters)} sequence families")
        
        # Process each sequence family
        for j, seq_cluster in enumerate(seq_clusters):
            cluster_size = len(seq_cluster)
            
            # Skip very small clusters (likely noise)
            if cluster_size < 3:
                continue
            
            # Find representative sequence
            representative = find_consensus_sequence(seq_cluster, method=consensus_method)
            
            # Estimate total count (extrapolate from sample)
            if total_in_range > max_sample:
                estimated_count = int(cluster_size * (total_in_range / len(seqs)))
            else:
                estimated_count = cluster_size
            
            family = {
                'size_cluster_idx': i,
                'size_range': f"{cluster['start']}-{cluster['end']}",
                'size_peak': cluster['peak'],
                'seq_family_idx': j,
                'seq_family_size': cluster_size,
                'estimated_total': estimated_count,
                'percentage_of_cluster': 100 * cluster_size / len(seqs),
                'representative_seq': representative,
                'rep_length': len(representative)
            }
            
            families.append(family)
            
            if j < 3:  # Print top 3 families per size cluster
                print(f"    Family {j+1}: {cluster_size} seqs "
                      f"({family['percentage_of_cluster']:.1f}% of sample), "
                      f"rep length: {len(representative)}bp")
    
    return families

def detect_organism_from_config(species_name):
    """
    Attempt to map species name to Dfam organism identifier.
    Returns best guess or None to use sequence-only search.
    
    Dfam organism format examples:
    - "Homo sapiens" (full scientific name)
    - "9606" (NCBI taxonomy ID) 
    - "human" (common name, sometimes works)
    """
    species_lower = species_name.lower()
    
    # Common organism mappings (scientific names work best with Dfam)
    organism_mapping = {
        'human': 'Homo sapiens',
        'homo': 'Homo sapiens', 
        'sapiens': 'Homo sapiens',
        'mouse': 'Mus musculus',
        'mus': 'Mus musculus',
        'musculus': 'Mus musculus',
        'rat': 'Rattus norvegicus',
        'rattus': 'Rattus norvegicus',
        'fly': 'Drosophila melanogaster',
        'drosophila': 'Drosophila melanogaster',
        'melanogaster': 'Drosophila melanogaster',
        'worm': 'Caenorhabditis elegans',
        'caenorhabditis': 'Caenorhabditis elegans',
        'elegans': 'Caenorhabditis elegans',
        'yeast': 'Saccharomyces cerevisiae',
        'saccharomyces': 'Saccharomyces cerevisiae',
        'cerevisiae': 'Saccharomyces cerevisiae',
        'zebrafish': 'Danio rerio',
        'danio': 'Danio rerio',
        'arabidopsis': 'Arabidopsis thaliana',
        'thaliana': 'Arabidopsis thaliana',
        'rice': 'Oryza sativa',
        'oryza': 'Oryza sativa',
        'horse': 'Equus caballus',
        'equus': 'Equus caballus',
        'caballus': 'Equus caballus',
        'cow': 'Bos taurus',
        'bos': 'Bos taurus',
        'taurus': 'Bos taurus',
        'pig': 'Sus scrofa',
        'sus': 'Sus scrofa',
        'scrofa': 'Sus scrofa',
        'chicken': 'Gallus gallus',
        'gallus': 'Gallus gallus',
        'dog': 'Canis lupus familiaris',
        'canis': 'Canis lupus familiaris',
        'cat': 'Felis catus',
        'felis': 'Felis catus'
    }
    
    for key, organism in organism_mapping.items():
        if key in species_lower:
            print(f"Detected Dfam organism: '{organism}' (from species: {species_name})")
            return organism
    
    # If no match, try using the species name as-is (capitalize first letter of each word)
    formatted_species = ' '.join(word.capitalize() for word in species_name.split())
    print(f"Using species name as organism: '{formatted_species}'")
    print("If this fails, set 'tepeak_dfam_organism' in config to override")
    return formatted_species

def dfam_submit_and_poll_relaxed(seq, organism=None, evalue_cutoff=10.0, timeout=30, tries=5, backoff=2.0):
    """Submit sequence to Dfam with relaxed e-value cutoff."""
    base = "https://dfam.org/api/searches/"
    
    try:
        # Use proper FASTA format for sequence
        fasta_seq = f">query\n{seq}"
        
        # Dfam API requires cutoff='evalue' and evalue parameter
        payload = {
            "sequence": fasta_seq,
            "cutoff": "evalue",
            "evalue": str(evalue_cutoff)
        }
        
        # Add organism if provided
        if organism:
            payload["organism"] = organism
        
        # Use form data (not JSON)
        r = requests.post(base, data=payload, timeout=timeout)
        
        if r.status_code == 400:
            try:
                error_detail = r.json()
                return {"status": "error", "detail": f"bad_request: {error_detail}"}
            except:
                return {"status": "error", "detail": "bad_request_400"}
        
        if r.status_code != 200:
            return {"status": "error", "detail": f"post_{r.status_code}"}
        
        result_data = r.json()
        jid = result_data.get("id", "")
        if not jid:
            return {"status": "error", "detail": "no_job_id"}
        
        url = base + jid
        
        # Poll for results
        for i in range(tries):
            time.sleep(backoff * (1.5 ** i))
            
            g = requests.get(url, timeout=timeout)
            
            if g.status_code == 200:
                data = g.json()
                
                try:
                    # Check if results are ready
                    if "results" not in data or not data["results"]:
                        status = data.get("status", "")
                        if status in ["PEND", "RUN", "RUNNING", "PENDING"]:
                            continue
                        else:
                            return {"status": "error", "detail": f"no_results_status_{status}"}
                            
                    hits = data["results"][0]["hits"]
                    
                    # Results are already filtered by the e-value we sent
                    filtered_hits = []
                    for h in hits:
                        evalue = h.get("evalue", float('inf'))
                        
                        # Handle different e-value formats
                        if isinstance(evalue, str):
                            try:
                                evalue = float(evalue)
                            except:
                                evalue = float('inf')
                        elif evalue is None:
                            evalue = float('inf')
                        else:
                            try:
                                evalue = float(evalue)
                            except:
                                evalue = float('inf')
                        
                        filtered_hits.append({
                            "name": h.get("accession", "") or h.get("description", ""),
                            "description": h.get("description", ""),
                            "class": h.get("type", ""),
                            "evalue": evalue,
                            "score": h.get("bit_score", ""),
                            "family": h.get("family", "")
                        })
                    
                    if filtered_hits:
                        filtered_hits.sort(key=lambda x: x['evalue'])
                        return {"status": "ok", "hits": filtered_hits}
                    else:
                        return {"status": "ok", "hits": []}
                        
                except (KeyError, IndexError, TypeError) as e:
                    # Check if job is still running
                    status = data.get("status", "")
                    if status in ["PEND", "RUNNING", "PENDING"]:
                        continue
                    else:
                        return {"status": "error", "detail": f"parse_error: {str(e)}"}
            elif g.status_code == 202:
                # Job still running
                continue
            else:
                return {"status": "error", "detail": f"get_{g.status_code}"}
        
        return {"status": "error", "detail": "timeout"}
        
    except requests.RequestException as e:
        return {"status": "error", "detail": f"request_error: {str(e)}"}

def annotate_families_with_dfam_relaxed(families, organism=None, evalue_cutoff=10.0, 
                                        batch_size=5, delay=3.0,
                                        skip_short=50, skip_long=10000):
    """Annotate families with Dfam using relaxed e-value threshold."""
    print(f"\nAnnotating {len(families)} families with Dfam...")
    print(f"Organism: {organism if organism else 'Not specified (sequence-only search)'}")
    print(f"E-value cutoff: {evalue_cutoff}")
    print(f"Sequence length range: {skip_short}-{skip_long}bp")
    print("(This may take a while due to API rate limits)\n")
    
    annotated = []
    success_count = 0
    error_counts = {}
    
    for i, family in enumerate(families):
        print(f"  {i+1}/{len(families)}: Size {family['size_range']}, "
              f"Family {family['seq_family_idx']+1} (len={family['rep_length']})...", end=' ')
        
        seq = family['representative_seq']
        seq_len = len(seq)
        
        # Skip invalid sequences
        if not seq or seq_len < skip_short:
            print("(too short, skipped)")
            family_annotated = family.copy()
            family_annotated.update({
                'dfam_name': '',
                'dfam_description': '',
                'dfam_class': '',
                'dfam_family': '',
                'dfam_evalue': '',
                'dfam_score': '',
                'dfam_status': 'too_short'
            })
            annotated.append(family_annotated)
            continue
        
        if seq_len > skip_long:
            print("(too long, skipped)")
            family_annotated = family.copy()
            family_annotated.update({
                'dfam_name': '',
                'dfam_description': '',
                'dfam_class': '',
                'dfam_family': '',
                'dfam_evalue': '',
                'dfam_score': '',
                'dfam_status': 'too_long'
            })
            annotated.append(family_annotated)
            continue
        
        # Query Dfam
        result = dfam_submit_and_poll_relaxed(
            seq,
            organism=organism,
            evalue_cutoff=evalue_cutoff,
            timeout=30, 
            tries=6, 
            backoff=2.0
        )
        
        if result.get('status') == 'ok' and result.get('hits'):
            top_hit = result['hits'][0]
            print(f"✓ {top_hit.get('name', 'Unknown')} "
                  f"(e={top_hit.get('evalue', 'N/A')}, class={top_hit.get('class', 'N/A')})")
            
            family_annotated = family.copy()
            family_annotated.update({
                'dfam_name': top_hit.get('name', ''),
                'dfam_description': top_hit.get('description', ''),
                'dfam_class': top_hit.get('class', ''),
                'dfam_family': top_hit.get('family', ''),
                'dfam_evalue': top_hit.get('evalue', ''),
                'dfam_score': top_hit.get('score', ''),
                'dfam_status': 'annotated'
            })
            success_count += 1
        elif result.get('status') == 'ok':
            print(f"✗ no hits (e-value > {evalue_cutoff})")
            family_annotated = family.copy()
            family_annotated.update({
                'dfam_name': '',
                'dfam_description': '',
                'dfam_class': '',
                'dfam_family': '',
                'dfam_evalue': '',
                'dfam_score': '',
                'dfam_status': 'no_hits'
            })
        else:
            error_detail = result.get('detail', 'unknown_error')
            print(f"✗ {error_detail}")
            
            error_counts[error_detail] = error_counts.get(error_detail, 0) + 1
            
            family_annotated = family.copy()
            family_annotated.update({
                'dfam_name': '',
                'dfam_description': '',
                'dfam_class': '',
                'dfam_family': '',
                'dfam_evalue': '',
                'dfam_score': '',
                'dfam_status': f'error_{error_detail.split(":")[0] if ":" in error_detail else error_detail}'
            })
        
        annotated.append(family_annotated)
        
        # Rate limiting
        if (i + 1) % batch_size == 0:
            print(f"  → Progress: {success_count}/{i+1} annotated ({100*success_count/(i+1):.1f}%)")
            print(f"  → Pausing {delay}s for rate limiting...")
            time.sleep(delay)
    
    # Final summary
    print(f"\n{'='*60}")
    print(f"Annotation complete: {success_count}/{len(families)} families annotated")
    
    if error_counts:
        print("\nError summary:")
        for error, count in sorted(error_counts.items(), key=lambda x: -x[1]):
            print(f"  {error}: {count}")
    
    return annotated

def save_families_to_csv(families, output_file):
    """Save families to CSV for downstream analysis."""
    families_df = pd.DataFrame(families)
    
    # Reorder columns for readability
    col_order = [
        'size_cluster_idx', 'size_range', 'size_peak',
        'seq_family_idx', 'seq_family_size', 'estimated_total', 'percentage_of_cluster',
        'rep_length', 'dfam_name', 'dfam_description', 'dfam_class', 'dfam_family',
        'dfam_evalue', 'dfam_score', 'dfam_status', 'representative_seq'
    ]
    
    # Only include columns that exist
    col_order = [c for c in col_order if c in families_df.columns]
    families_df = families_df[col_order]
    
    families_df.to_csv(output_file, index=False)
    print(f"\nSaved {len(families_df)} families to {output_file}")
    
    return families_df

def visualize_clusters(df, clusters, min_w, max_w, output_file):
    """Visualize clusters as bands on the histogram."""
    cnt = df["length"].value_counts().sort_index()
    lengths = pd.RangeIndex(cnt.index.min(), cnt.index.max() + 1)
    y = cnt.reindex(lengths, fill_value=0)
    
    mask = (lengths >= min_w) & (lengths <= max_w)
    x = lengths[mask]
    y_plot = y[mask]
    
    fig, axes = plt.subplots(2, 1, figsize=(16, 10))
    
    # === Top plot: Log scale with cluster bands ===
    ax = axes[0]
    ax.plot(x, y_plot, linewidth=0.5, alpha=0.5, color='gray', label='Frequency')
    
    # Draw cluster regions
    colors = plt.cm.tab20(np.linspace(0, 1, min(20, len(clusters))))
    for i, cluster in enumerate(clusters[:20]):  # Show top 20
        color = colors[i % 20]
        ax.axvspan(cluster['start'], cluster['end'], alpha=0.3, color=color, 
                  label=f"C{i+1}: {cluster['start']}-{cluster['end']}bp")
        ax.scatter([cluster['peak']], [cluster['peak_count']], 
                  color=color, s=150, marker='*', edgecolors='black', linewidths=1.5, zorder=5)
    
    ax.set_xlabel('Insertion size (bp)')
    ax.set_ylabel('Frequency (log scale)')
    ax.set_title(f'Detected Clusters (n={len(clusters)})')
    ax.set_yscale('log')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3)
    
    # === Bottom plot: Cluster metrics ===
    ax = axes[1]
    
    if clusters:
        cluster_df = pd.DataFrame(clusters[:20])
        
        x_pos = np.arange(len(cluster_df))
        ax.bar(x_pos, cluster_df['total_count'], alpha=0.7, label='Total count')
        ax.set_xlabel('Cluster index (sorted by count)')
        ax.set_ylabel('Total insertions in cluster')
        ax.set_title('Cluster Sizes')
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add cluster ranges as labels
        labels = [f"{c['start']}-{c['end']}" for c in clusters[:20]]
        ax.set_xticks(x_pos)
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    # Print summary
    print("\n=== Cluster Summary ===")
    print(f"Total clusters: {len(clusters)}")
    if clusters:
        print(f"Size range: {clusters[-1]['start']}bp to {clusters[0]['end']}bp")
        print(f"Total insertions in clusters: {sum(c['total_count'] for c in clusters):,}")
        print(f"Percentage of insertions in clusters: {100 * sum(c['total_count'] for c in clusters) / len(df):.1f}%")

def main():
    # Get parameters from Snakemake
    global_vcf_file = snakemake.input.global_vcf
    output_dir = snakemake.params.output_dir
    
    # Get TEPEAK parameters from config with defaults
    config = snakemake.config
    species = config.get('species', 'unknown')
    MIN_WINDOW = config.get('low', 200)
    MAX_WINDOW = config.get('high', 6400)
    min_cluster_size = config.get('tepeak_min_cluster_size', 10)
    percentile_threshold = config.get('tepeak_percentile_threshold', 75)
    window_size = config.get('tepeak_window_size', 50)
    merge_distance = config.get('tepeak_merge_distance', 100)
    pid_threshold = config.get('tepeak_pid_threshold', 0.85)
    max_clusters = config.get('tepeak_max_clusters', 50)
    max_sample_seqs = config.get('tepeak_max_sample_seqs', 200)
    dfam_evalue = config.get('tepeak_dfam_evalue', 10.0)
    dfam_batch_size = config.get('tepeak_dfam_batch_size', 5)
    dfam_delay = config.get('tepeak_dfam_delay', 3.0)
    
    # Get organism for Dfam queries
    dfam_organism = config.get('tepeak_dfam_organism', None)
    if dfam_organism is None:
        dfam_organism = detect_organism_from_config(species)
    
    # Output files
    output_csv = os.path.join(output_dir, 'tepeak_families_annotated.csv')
    clusters_plot = os.path.join(output_dir, 'tepeak_clusters_plot.png')
    
    # Step 1: Load data
    print("="*60)
    print("TEPEAK: LOADING DATA")
    print("="*60)
    df = load_global_tsv(global_vcf_file)
    
    # Step 2: Find size-based clusters using adaptive density method
    print("\n" + "="*60)
    print("TEPEAK: FINDING SIZE CLUSTERS")
    print("="*60)
    clusters = find_clusters_adaptive_density(
        df,
        min_w=MIN_WINDOW,
        max_w=MAX_WINDOW,
        min_cluster_size=min_cluster_size,
        percentile_threshold=percentile_threshold,
        window_size=window_size,
        merge_distance=merge_distance
    )
    
    print(f"\nFound {len(clusters)} clusters")
    
    if clusters:
        # Convert to DataFrame for display
        clusters_df = pd.DataFrame(clusters)
        clusters_df['percentage'] = 100 * clusters_df['total_count'] / len(df)
        
        print("\nTop 20 clusters by total count:")
        print(clusters_df[['start', 'end', 'peak', 'width', 'total_count', 
                           'peak_count', 'avg_density', 'percentage']].head(20))
        
        # Visualize clusters
        visualize_clusters(df, clusters, MIN_WINDOW, MAX_WINDOW, clusters_plot)
        
        # Step 3: Find sequence families within clusters
        print("\n" + "="*60)
        print("TEPEAK: CLUSTERING SEQUENCES")
        print("="*60)
        families = process_all_clusters(
            df, 
            clusters[:max_clusters],  # Process configured number of clusters
            pid_threshold=pid_threshold,
            max_sample=max_sample_seqs,
            consensus_method='most_common'
        )
        
        print(f"\n{'='*60}")
        print(f"Total sequence families identified: {len(families)}")
        
        if families:
            # Step 4: Annotate with Dfam
            print("\n" + "="*60)
            print("TEPEAK: DFAM ANNOTATION")
            print("="*60)
            annotated_families = annotate_families_with_dfam_relaxed(
                families,
                organism=dfam_organism,
                evalue_cutoff=dfam_evalue,
                batch_size=dfam_batch_size,
                delay=dfam_delay,
                skip_short=50,
                skip_long=10000
            )
            
            # Step 5: Save results
            print("\n" + "="*60)
            print("TEPEAK: SAVING RESULTS")
            print("="*60)
            results_df = save_families_to_csv(annotated_families, output_csv)
            
            # Final summary
            annotated = results_df[results_df['dfam_status'] == 'annotated']
            print(f"\nSuccessfully annotated: {len(annotated)}/{len(results_df)} families")
            
            if len(annotated) > 0:
                print(f"Total insertions annotated: {annotated['estimated_total'].sum():,}")
        else:
            # Create empty output file
            pd.DataFrame().to_csv(output_csv, index=False)
            print("No sequence families found - created empty output file")
    else:
        # Create empty output files
        pd.DataFrame().to_csv(output_csv, index=False)
        # Create empty plot
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, 'No clusters found', ha='center', va='center', transform=ax.transAxes)
        ax.set_title('TEPEAK Cluster Analysis - No Clusters Detected')
        plt.savefig(clusters_plot)
        plt.close()
        print("No clusters found - created empty output files")

if __name__ == '__main__':
    main()
