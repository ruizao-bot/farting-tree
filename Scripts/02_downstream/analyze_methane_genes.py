#!/usr/bin/env python3
"""
Analyze DIAMOND methane gene search results
Extract species abundance and functional gene abundance
"""

import pandas as pd
import sys
import os
from collections import Counter
import re

def parse_diamond_results(input_file):
    """
    Parse DIAMOND blastx results
    Format: qseqid sseqid pident length evalue bitscore stitle
    """
    columns = ['query_id', 'subject_id', 'pident', 'length', 'evalue', 'bitscore', 'description']
    
    df = pd.read_csv(input_file, sep='\t', names=columns, header=None)
    
    return df

def extract_gene_info(description):
    """Extract gene name and species from description"""
    # Pattern: protein_name [species_name]
    match = re.match(r'^(.+?)\s+\[(.+?)\]', description)
    if match:
        gene_name = match.group(1).strip()
        species = match.group(2).strip()
        return gene_name, species
    return description, "Unknown"

def parse_singlem_otu_table(singlem_file):
    """
    Parse SingleM OTU table for taxonomic abundance
    Returns DataFrame with taxonomy and abundance information
    """
    try:
        df = pd.read_csv(singlem_file, sep='\t')
        # Extract species-level taxonomy
        df['full_taxonomy'] = df['taxonomy']
        
        # Parse taxonomy string to extract different levels
        def extract_taxonomy_levels(tax_string):
            levels = {'domain': 'Unknown', 'phylum': 'Unknown', 'class': 'Unknown', 
                     'order': 'Unknown', 'family': 'Unknown', 'genus': 'Unknown', 'species': 'Unknown'}
            
            if pd.isna(tax_string) or tax_string == 'Root':
                return pd.Series(levels)
            
            parts = [p.strip() for p in tax_string.split(';')]
            for part in parts:
                if part.startswith('d__'):
                    levels['domain'] = part.replace('d__', '')
                elif part.startswith('p__'):
                    levels['phylum'] = part.replace('p__', '')
                elif part.startswith('c__'):
                    levels['class'] = part.replace('c__', '')
                elif part.startswith('o__'):
                    levels['order'] = part.replace('o__', '')
                elif part.startswith('f__'):
                    levels['family'] = part.replace('f__', '')
                elif part.startswith('g__'):
                    levels['genus'] = part.replace('g__', '')
                elif part.startswith('s__'):
                    levels['species'] = part.replace('s__', '')
            
            return pd.Series(levels)
        
        # Apply taxonomy parsing
        tax_df = df['taxonomy'].apply(extract_taxonomy_levels)
        df = pd.concat([df, tax_df], axis=1)
        
        return df
    except FileNotFoundError:
        print(f"WARNING: SingleM file not found: {singlem_file}")
        return None
    except Exception as e:
        print(f"WARNING: Error parsing SingleM file: {e}")
        return None

def classify_methane_gene(description):
    """Classify genes into functional categories"""
    desc_lower = description.lower()
    
    # Methane oxidation (methanotrophs)
    if any(x in desc_lower for x in ['particulate methane monooxygenase', 'pmoa']):
        return 'pMMO (Particulate methane monooxygenase)'
    elif any(x in desc_lower for x in ['soluble methane monooxygenase', 'mmox', 'mmoy', 'mmoz']):
        return 'sMMO (Soluble methane monooxygenase)'
    elif 'mmor' in desc_lower:
        return 'MmoR (MMO regulatory protein)'
    elif 'mmog' in desc_lower:
        return 'MmoG (MMO chaperone)'
    
    # Methane production (methanogens)
    elif any(x in desc_lower for x in ['methyl-coenzyme m reductase', 'mcra', 'mcrb', 'mcrg']):
        return 'MCR (Methyl-coenzyme M reductase)'
    
    # Related metabolism
    elif 'methanol dehydrogenase' in desc_lower:
        return 'MDH (Methanol dehydrogenase)'
    elif 'formaldehyde' in desc_lower:
        return 'Formaldehyde metabolism'
    elif 'formate' in desc_lower:
        return 'Formate metabolism'
    
    # Supporting enzymes
    elif any(x in desc_lower for x in ['ammonia', 'ammonium']):
        return 'Nitrogen metabolism (related)'
    
    return 'Other'

def analyze_gene_abundance(df):
    """Analyze functional gene abundance"""
    print("\n" + "="*80)
    print("FUNCTIONAL GENE ABUNDANCE ANALYSIS")
    print("="*80 + "\n")
    
    # Extract gene info
    df['gene_name'], df['species'] = zip(*df['description'].apply(extract_gene_info))
    df['gene_category'] = df['description'].apply(classify_methane_gene)
    
    # Filter for methane-specific genes
    methane_genes = df[df['gene_category'] != 'Other'].copy()
    
    print(f"Total DIAMOND hits: {len(df):,}")
    print(f"Methane-related hits: {len(methane_genes):,}\n")
    
    # Gene category abundance (unique reads)
    print("1. GENE CATEGORY ABUNDANCE (by unique reads)")
    print("-" * 60)
    category_counts = methane_genes.groupby('gene_category')['query_id'].nunique()
    category_counts_sorted = category_counts.sort_values(ascending=False)
    
    for category, count in category_counts_sorted.items():
        pct = (count / len(methane_genes['query_id'].unique())) * 100
        print(f"  {category:<45} {count:>6} reads ({pct:>5.2f}%)")
    
    # Top genes by read count
    print("\n2. TOP 20 METHANE-RELATED GENES (by read count)")
    print("-" * 60)
    gene_counts = methane_genes.groupby('gene_name')['query_id'].nunique().sort_values(ascending=False).head(20)
    
    for gene, count in gene_counts.items():
        print(f"  {gene[:60]:<62} {count:>6} reads")
    
    # Gene presence/absence
    print("\n3. KEY METHANE METABOLISM GENES DETECTED")
    print("-" * 60)
    
    key_genes = {
        'pmoA': ['pmoa', 'particulate methane monooxygenase'],
        'mmoX': ['mmox', 'soluble methane monooxygenase'],
        'mcrA': ['mcra', 'methyl-coenzyme m reductase'],
        'mxaF': ['mxaf', 'methanol dehydrogenase'],
    }
    
    for gene_key, keywords in key_genes.items():
        matches = df[df['description'].str.lower().str.contains('|'.join(keywords))]
        if len(matches) > 0:
            unique_reads = matches['query_id'].nunique()
            avg_pident = matches['pident'].mean()
            print(f"  ✓ {gene_key:<8} Detected: {unique_reads:>5} reads (avg identity: {avg_pident:.1f}%)")
        else:
            print(f"  ✗ {gene_key:<8} Not detected")
    
    return methane_genes

def analyze_species_abundance(df, singlem_df=None):
    """Analyze species/organism abundance using SingleM data if available"""
    print("\n" + "="*80)
    print("SPECIES/ORGANISM ABUNDANCE ANALYSIS")
    print("="*80 + "\n")
    
    if singlem_df is not None:
        print("Using SingleM taxonomic classification\n")
        
        # Aggregate by species
        species_abundance = singlem_df.groupby('species').agg({
            'num_hits': 'sum',
            'coverage': 'sum'
        }).sort_values('num_hits', ascending=False)
        
        # Also get genus-level for unclassified species
        genus_abundance = singlem_df[singlem_df['species'] == 'Unknown'].groupby('genus').agg({
            'num_hits': 'sum',
            'coverage': 'sum'
        }).sort_values('num_hits', ascending=False)
        
        print(f"Total species detected: {len(species_abundance[species_abundance.index != 'Unknown'])}")
        print(f"Total genera detected: {singlem_df['genus'].nunique()}\n")
        
        print("1. TOP 30 SPECIES BY READ ABUNDANCE (from SingleM)")
        print("-" * 80)
        
        total_hits = singlem_df['num_hits'].sum()
        
        for i, (species, row) in enumerate(species_abundance.head(30).iterrows(), 1):
            if species != 'Unknown':
                pct = (row['num_hits'] / total_hits) * 100
                print(f"  {i:>2}. {species[:60]:<62} {int(row['num_hits']):>6} reads ({pct:>5.2f}%)")
        
        # Show top genera for unclassified species
        if len(genus_abundance) > 0:
            print("\n   TOP GENERA (species-level unclassified)")
            print("   " + "-" * 76)
            for i, (genus, row) in enumerate(genus_abundance.head(10).iterrows(), 1):
                if genus != 'Unknown':
                    pct = (row['num_hits'] / total_hits) * 100
                    print(f"      {genus[:60]:<62} {int(row['num_hits']):>6} reads ({pct:>5.2f}%)")
        
        species_reads = species_abundance
    else:
        print("Using DIAMOND protein annotation (less accurate for taxonomy)\n")
        
        # Extract species from DIAMOND results
        df['gene_name'], df['species'] = zip(*df['description'].apply(extract_gene_info))
        
        # Count unique reads per species
        species_reads = df.groupby('species')['query_id'].nunique().sort_values(ascending=False)
        
        print(f"Total species/organisms detected: {len(species_reads)}\n")
        
        print("1. TOP 30 SPECIES BY READ ABUNDANCE (from DIAMOND)")
        print("-" * 80)
        
        total_unique_reads = df['query_id'].nunique()
        
        for i, (species, count) in enumerate(species_reads.head(30).items(), 1):
            pct = (count / total_unique_reads) * 100
            print(f"  {i:>2}. {species[:60]:<62} {count:>6} reads ({pct:>5.2f}%)")
    
    # Methanotroph species
    print("\n2. KNOWN METHANOTROPH SPECIES")
    print("-" * 80)
    
    methanotroph_keywords = [
        'methylosinus', 'methylocystis', 'methylomonas', 'methylobacter',
        'methylococcus', 'methylomicrobium', 'methylocaldum', 'methylocapsa'
    ]
    
    if singlem_df is not None:
        # Check both species and genus columns
        methanotrophs = singlem_df[
            singlem_df['species'].str.lower().str.contains('|'.join(methanotroph_keywords), na=False) |
            singlem_df['genus'].str.lower().str.contains('|'.join(methanotroph_keywords), na=False)
        ]
        
        if len(methanotrophs) > 0:
            methano_abundance = methanotrophs.groupby(['genus', 'species'])['num_hits'].sum().sort_values(ascending=False)
            for (genus, species), count in methano_abundance.items():
                display_name = f"{genus} {species}" if species != 'Unknown' else genus
                print(f"  ✓ {display_name:<60} {int(count):>6} reads")
        else:
            print("  No known methanotroph species detected")
    else:
        methanotrophs = df[df['species'].str.lower().str.contains('|'.join(methanotroph_keywords), na=False)]
        
        if len(methanotrophs) > 0:
            methano_species = methanotrophs.groupby('species')['query_id'].nunique().sort_values(ascending=False)
            for species, count in methano_species.items():
                print(f"  ✓ {species:<60} {count:>6} reads")
        else:
            print("  No known methanotroph species detected")
    
    # Methanogen species
    print("\n3. KNOWN METHANOGEN SPECIES")
    print("-" * 80)
    
    methanogen_keywords = [
        'methanobacterium', 'methanobrevibacter', 'methanosarcina',
        'methanococcus', 'methanothermobacter', 'methanospirillum'
    ]
    
    if singlem_df is not None:
        methanogens = singlem_df[
            singlem_df['species'].str.lower().str.contains('|'.join(methanogen_keywords), na=False) |
            singlem_df['genus'].str.lower().str.contains('|'.join(methanogen_keywords), na=False)
        ]
        
        if len(methanogens) > 0:
            methano_abundance = methanogens.groupby(['genus', 'species'])['num_hits'].sum().sort_values(ascending=False)
            for (genus, species), count in methano_abundance.items():
                display_name = f"{genus} {species}" if species != 'Unknown' else genus
                print(f"  ✓ {display_name:<60} {int(count):>6} reads")
        else:
            print("  No known methanogen species detected")
    else:
        methanogens = df[df['species'].str.lower().str.contains('|'.join(methanogen_keywords), na=False)]
        
        if len(methanogens) > 0:
            methano_species = methanogens.groupby('species')['query_id'].nunique().sort_values(ascending=False)
            for species, count in methano_species.items():
                print(f"  ✓ {species:<60} {count:>6} reads")
        else:
            print("  No known methanogen species detected")
    
    return species_reads

def calculate_rpkm(df, total_reads_millions):
    """Calculate RPKM for genes"""
    print("\n" + "="*80)
    print("GENE ABUNDANCE NORMALIZATION (RPKM)")
    print("="*80 + "\n")
    
    df['gene_name'], df['species'] = zip(*df['description'].apply(extract_gene_info))
    df['gene_category'] = df['description'].apply(classify_methane_gene)
    
    methane_genes = df[df['gene_category'] != 'Other'].copy()
    
    # Group by gene name
    gene_stats = methane_genes.groupby('gene_name').agg({
        'query_id': 'nunique',  # Unique reads
        'length': 'mean',        # Average alignment length
        'pident': 'mean'         # Average identity
    }).rename(columns={'query_id': 'read_count'})
    
    # Calculate RPKM: (read_count * 1000) / (avg_length * total_reads_millions)
    gene_stats['rpkm'] = (gene_stats['read_count'] * 1000) / (gene_stats['length'] * total_reads_millions)
    gene_stats = gene_stats.sort_values('rpkm', ascending=False)
    
    print(f"Total reads (millions): {total_reads_millions:.2f}M\n")
    print("TOP 20 GENES BY RPKM")
    print("-" * 100)
    print(f"{'Gene':<50} {'Reads':<10} {'Avg Len':<10} {'Avg ID%':<10} {'RPKM':<10}")
    print("-" * 100)
    
    for gene, row in gene_stats.head(20).iterrows():
        print(f"{gene[:48]:<50} {row['read_count']:<10} {row['length']:<10.1f} {row['pident']:<10.1f} {row['rpkm']:<10.4f}")
    
    return gene_stats

def main():
    if len(sys.argv) < 3:
        print("Usage: python analyze_methane_genes.py <diamond_results.txt> <singlem_otu_table.csv> [total_reads_millions]")
        print("\nExample:")
        print("  python analyze_methane_genes.py 53394_combined_methane_hits.txt Data/processed_data/singlem_otu_table.csv")
        print("  python analyze_methane_genes.py 53394_combined_methane_hits.txt Data/processed_data/singlem_otu_table.csv 8.79")
        print("\nNote: If total_reads_millions is not provided, it will be auto-detected from bbduk_stats.txt")
        sys.exit(1)
    
    input_file = sys.argv[1]
    sample_id = input_file.split('/')[-1].split('_')[0]
    singlem_file = sys.argv[2]
    
    # Auto-detect or use provided total reads
    if len(sys.argv) > 3:
        total_reads_millions = float(sys.argv[3])
        print(f"Using provided total reads: {total_reads_millions:.2f}M")
    else:
        # Try to read from bbduk stats
        base_dir = '/'.join(input_file.split('/')[:-3])
        bbduk_stats = f"{base_dir}/processed_data/bbduk_cleaned/{sample_id}_bbduk_stats.txt"
        if os.path.exists(bbduk_stats):
            with open(bbduk_stats, 'r') as f:
                for line in f:
                    if line.startswith('#Total'):
                        total_reads = int(line.split()[1])
                        total_reads_millions = total_reads / 1_000_000
                        print(f"Auto-detected from {bbduk_stats}: {total_reads_millions:.2f}M reads")
                        break
                else:
                    print(f"ERROR: Could not find #Total line in {bbduk_stats}")
                    sys.exit(1)
        else:
            print(f"ERROR: BBduk stats file not found: {bbduk_stats}")
            print("Please provide total_reads_millions as third argument")
            sys.exit(1)
    
    print("\n" + "="*80)
    print(f"METHANE GENE ANALYSIS - Sample {sample_id}")
    print("="*80)
    print(f"DIAMOND results: {input_file}")
    print(f"Total reads: {total_reads_millions:.2f} million")
    print(f"SingleM OTU table: {singlem_file}")
    print("="*80)
    
    # Parse DIAMOND results
    df = parse_diamond_results(input_file)
    
    # Parse SingleM results
    singlem_df = parse_singlem_otu_table(singlem_file)
    if singlem_df is None:
        print("ERROR: Failed to load SingleM OTU table")
        sys.exit(1)
    print(f"✓ Loaded {len(singlem_df)} SingleM OTU records\n")
    
    # Analyze gene abundance
    methane_genes = analyze_gene_abundance(df)
    
    # Analyze species abundance using SingleM
    species_reads = analyze_species_abundance(df, singlem_df)
    
    # Calculate RPKM
    gene_stats = calculate_rpkm(df, total_reads_millions)
    
    # Set output directory to Results/quick_search
    output_dir = "Results/quick_search"
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Export detailed tables
    methane_genes_file = f"{output_dir}/{sample_id}_methane_genes_detailed.csv"
    methane_genes.to_csv(methane_genes_file, index=False)
    
    species_file = f"{output_dir}/{sample_id}_species_abundance.csv"
    species_reads.to_csv(species_file)
    
    gene_stats_file = f"{output_dir}/{sample_id}_gene_rpkm.csv"
    gene_stats.to_csv(gene_stats_file)
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print(f"\nOutput files:")
    print(f"  1. Methane genes detail: {methane_genes_file}")
    print(f"  2. Species abundance: {species_file}")
    print(f"  3. Gene RPKM: {gene_stats_file}")
    print("\n")

if __name__ == "__main__":
    main()
