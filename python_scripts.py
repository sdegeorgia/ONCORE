# =============================================================================
# scripts/01_prepare_cell_lines.py
# =============================================================================

#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import os

def main():
    # Read DepMap data
    depmap_file = snakemake.input["depmap_data"]
    
    if os.path.exists(depmap_file):
        depmap = pd.read_csv(depmap_file, index_col=0)
    else:
        print("Warning: DepMap data not found, creating mock metadata")
        depmap = pd.DataFrame()
    
    # Parameters
    cancer_type = snakemake.params["cancer_type"]
    target_gene = snakemake.params["target_gene"]
    cell_lines = snakemake.params["cell_lines"]
    
    # Create metadata dataframe
    metadata = pd.DataFrame({
        'cell_line': cell_lines,
        'cancer_type': [cancer_type] * len(cell_lines),
        'target_gene': [target_gene] * len(cell_lines)
    })
    
    # Add additional metadata if DepMap data available
    if not depmap.empty and target_gene in depmap.columns:
        gene_effects = []
        for cell_line in cell_lines:
            matching_rows = depmap.index[depmap.index.str.contains(cell_line, case=False, na=False)]
            if len(matching_rows) > 0:
                gene_effects.append(depmap.loc[matching_rows[0], target_gene])
            else:
                gene_effects.append(np.nan)
        
        metadata['gene_effect_score'] = gene_effects
        metadata['essential'] = metadata['gene_effect_score'] < -0.5
    else:
        metadata['gene_effect_score'] = np.nan
        metadata['essential'] = False
    
    # Add data availability flags (would be implemented with real API calls)
    metadata['chip_seq_available'] = True  # Mock data
    metadata['rna_seq_available'] = True   # Mock data
    metadata['atac_seq_available'] = [True, True, False, True, False][:len(cell_lines)]
    
    # Save metadata
    metadata.to_csv(snakemake.output["metadata"], sep='\t', index=False)
    
    print(f"Prepared metadata for {len(cell_lines)} {cancer_type} cell lines")
    print(f"Target gene: {target_gene}")

if __name__ == "__main__":
    main()

# =============================================================================
# scripts/02_consensus_peaks.py
# =============================================================================

#!/usr/bin/env python3

import pandas as pd
import numpy as np
from collections import defaultdict
import pybedtools

def main():
    peak_files = snakemake.input["peaks"]
    min_overlap = snakemake.params["min_overlap"]
    output_file = snakemake.output["consensus"]
    
    # Read all peak files
    all_peaks = []
    for i, peak_file in enumerate(peak_files):
        if os.path.getsize(peak_file) > 0:  # Check if file is not empty
            peaks = pd.read_csv(peak_file, sep='\t', header=None,
                              names=['chr', 'start', 'end', 'name', 'score'])
            peaks['cell_line'] = f"cell_{i}"
            all_peaks.append(peaks)
    
    if not all_peaks:
        # Create empty output file
        open(output_file, 'w').close()
        return
    
    # Combine all peaks
    combined_peaks = pd.concat(all_peaks, ignore_index=True)
    
    # Use bedtools to find overlapping regions
    peak_bed = pybedtools.BedTool.from_dataframe(
        combined_peaks[['chr', 'start', 'end', 'name', 'score']]
    )
    
    # Merge overlapping peaks
    merged_peaks = peak_bed.merge(c=4, o='count')
    
    # Filter by minimum overlap requirement
    min_cell_lines = max(1, int(len(peak_files) * min_overlap))
    
    consensus_regions = []
    for interval in merged_peaks:
        if int(interval[3]) >= min_cell_lines:
            consensus_regions.append([
                interval[0],  # chr
                int(interval[1]),  # start
                int(interval[2]),  # end
                f"consensus_peak_{len(consensus_regions) + 1}",  # name
                int(interval[3]) * 100,  # score (number of cell lines * 100)
                "."  # strand
            ])
    
    # Save consensus peaks
    if consensus_regions:
        consensus_df = pd.DataFrame(consensus_regions,
                                  columns=['chr', 'start', 'end', 'name', 'score', 'strand'])
        consensus_df.to_csv(output_file, sep='\t', header=False, index=False)
        print(f"Identified {len(consensus_regions)} consensus peaks")
    else:
        # Create empty file
        open(output_file, 'w').close()
        print("No consensus peaks found")

if __name__ == "__main__":
    main()

# =============================================================================
# scripts/04_identify_cancer_tfs.py
# =============================================================================

#!/usr/bin/env python3

import pandas as pd
import requests
import json

def get_cancer_tfs(cancer_type):
    """
    Get cancer-specific transcription factors from various databases
    """
    
    # Define cancer-specific TF databases (simplified version)
    cancer_tf_db = {
        'melanoma': ['SOX10', 'MITF', 'PAX3', 'TFAP2A', 'LEF1', 'FOXD3'],
        'breast': ['ESR1', 'FOXA1', 'GATA3', 'PGR', 'AP1', 'TFAP2C'],
        'lung': ['NKX2-1', 'FOXA2', 'GATA6', 'HOXA5', 'TBX4', 'FOXP2'],
        'colorectal': ['CDX2', 'FOXF1', 'GATA4', 'HNF4A', 'TCF7L2', 'MSX1'],
        'prostate': ['AR', 'FOXA1', 'HOXB13', 'NKX3-1', 'ERG', 'ETV1'],
        'pancreatic': ['PDX1', 'FOXA2', 'GATA6', 'HNF1A', 'KLF4', 'MEIS1']
    }
    
    # Oncogenes and tumor suppressors (common across cancers)
    common_oncogenes = ['MYC', 'JUN', 'FOS', 'ETS1', 'STAT3', 'NFKB1']
    tumor_suppressors = ['TP53', 'RB1', 'CDKN2A', 'APC', 'BRCA1', 'PTEN']
    
    # Get cancer-specific TFs
    cancer_specific = cancer_tf_db.get(cancer_type.lower(), [])
    
    # Combine all TFs
    all_tfs = list(set(cancer_specific + common_oncogenes + tumor_suppressors))
    
    return all_tfs

def create_tf_metadata(tf_list, cancer_type, target_gene):
    """
    Create metadata for transcription factors
    """
    
    tf_metadata = []
    for tf in tf_list:
        # Determine TF category
        if tf in ['MYC', 'JUN', 'FOS', 'ETS1', 'STAT3', 'NFKB1']:
            category = 'oncogene'
        elif tf in ['TP53', 'RB1', 'CDKN2A', 'APC', 'BRCA1', 'PTEN']:
            category = 'tumor_suppressor'
        else:
            category = 'tissue_specific'
        
        # Priority scoring (simplified)
        if tf == target_gene:
            priority = 1
        elif category == 'tissue_specific':
            priority = 2
        elif category == 'oncogene':
            priority = 3
        else:
            priority = 4
        
        tf_metadata.append({
            'tf_name': tf,
            'category': category,
            'priority': priority,
            'cancer_type': cancer_type,
            'druggable': tf in ['ESR1', 'AR', 'MYC', 'STAT3']  # Simplified
        })
    
    return pd.DataFrame(tf_metadata)

def main():
    cancer_type = snakemake.params["cancer_type"]
    target_gene = snakemake.params["target_gene"]
    
    # Get cancer-specific TFs
    tf_list = get_cancer_tfs(cancer_type)
    
    # Add target gene if not already present
    if target_gene not in tf_list:
        tf_list.append(target_gene)
    
    # Create metadata
    tf_metadata = create_tf_metadata(tf_list, cancer_type, target_gene)
    
    # Sort by priority
    tf_metadata = tf_metadata.sort_values('priority')
    
    # Save results
    tf_metadata.to_csv(snakemake.output["tf_list"], sep='\t', index=False)
    
    print(f"Identified {len(tf_list)} transcription factors for {cancer_type}")
    print("Top 5 priority TFs:")
    print(tf_metadata.head()['tf_name'].tolist())

if __name__ == "__main__":
    main()

# =============================================================================
# scripts/06_therapeutic_targets.py
# =============================================================================

#!/usr/bin/env python3

import pandas as pd
import numpy as np

def score_therapeutic_potential(tf_binding, expression, mutations):
    """
    Score therapeutic potential based on multiple criteria
    """
    
    # Read input data
    if os.path.getsize(tf_binding) > 0:
        binding_df = pd.read_csv(tf_binding, sep='\t', header=None,
                                names=['chr', 'start', 'end', 'name', 'score', 'strand'])
    else:
        binding_df = pd.DataFrame()
    
    expr_df = pd.read_csv(expression, sep='\t')
    mut_df = pd.read_csv(mutations, sep='\t')
    
    # Extract TF names from binding sites
    if not binding_df.empty:
        binding_df['tf'] = binding_df['name'].str.split('_').str[0]
        tf_counts = binding_df['tf'].value_counts()
    else:
        tf_counts = pd.Series()
    
    # Define druggable TFs (simplified database)
    druggable_tfs = {
        'ESR1': {'drug': 'Tamoxifen', 'mechanism': 'Antagonist', 'approved': True},
        'AR': {'drug': 'Enzalutamide', 'mechanism': 'Antagonist', 'approved': True},
        'MYC': {'drug': 'BET_inhibitors', 'mechanism': 'Indirect', 'approved': False},
        'STAT3': {'drug': 'Stattic', 'mechanism': 'Direct', 'approved': False},
        'NFKB1': {'drug': 'Bortezomib', 'mechanism': 'Indirect', 'approved': True},
        'SOX10': {'drug': 'Sox10_inhibitor', 'mechanism': 'Direct', 'approved': False}
    }
    
    # Score each TF
    therapeutic_targets = []
    for tf in tf_counts.index:
        if tf in druggable_tfs:
            # Binding score (normalized)
            binding_score = min(tf_counts[tf] / tf_counts.max(), 1.0) if tf_counts.max() > 0 else 0
            
            # Expression correlation score (mock)
            expr_score = np.random.uniform(0.3, 0.9)  # Would use real correlation data
            
            # Mutation frequency score (mock)
            mut_score = np.random.uniform(0.1, 0.8)   # Would use real mutation data
            
            # Clinical readiness score
            clinical_score = 1.0 if druggable_tfs[tf]['approved'] else 0.5
            
            # Combined score
            combined_score = (binding_score * 0.3 + expr_score * 0.3 + 
                            mut_score * 0.2 + clinical_score * 0.2)
            
            therapeutic_targets.append({
                'tf_name': tf,
                'binding_sites': tf_counts[tf],
                'binding_score': binding_score,
                'expression_score': expr_score,
                'mutation_score': mut_score,
                'clinical_score': clinical_score,
                'combined_score': combined_score,
                'drug_name': druggable_tfs[tf]['drug'],
                'mechanism': druggable_tfs[tf]['mechanism'],
                'approved': druggable_tfs[tf]['approved'],
                'priority': 'High' if combined_score > 0.7 else 'Medium' if combined_score > 0.4 else 'Low'
            })
    
    return pd.DataFrame(therapeutic_targets)

def main():
    # Score therapeutic potential
    targets_df = score_therapeutic_potential(
        snakemake.input["tf_binding"],
        snakemake.input["expression"],
        snakemake.input["mutations"]
    )
    
    if not targets_df.empty:
        # Sort by combined score
        targets_df = targets_df.sort_values('combined_score', ascending=False)
        
        # Save results
        targets_df.to_csv(snakemake.output["targets"], sep='\t', index=False)
        
        print(f"Identified {len(targets_df)} potential therapeutic targets")
        print("Top 3 targets:")
        for _, row in targets_df.head(3).iterrows():
            print(f"  {row['tf_name']}: {row['drug_name']} (Score: {row['combined_score']:.2f})")
    else:
        # Create empty file
        pd.DataFrame().to_csv(snakemake.output["targets"], sep='\t', index=False)
        print("No therapeutic targets identified")

if __name__ == "__main__":
    main()