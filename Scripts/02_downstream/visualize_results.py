#!/usr/bin/env python3
"""
Visualize methane gene analysis results
Creates contribution plots and comparative analyses
"""

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.patches import FancyBboxPatch
from scipy.spatial import procrustes

# Set style
sns.set_style("whitegrid")
plt.rcParams['font.size'] = 10
plt.rcParams['figure.dpi'] = 300

def plot_gene_category_contribution(detailed_df, sample_id, output_dir):
    """Plot stacked bar chart of gene category contributions"""
    
    # Count reads per gene category
    category_counts = detailed_df['gene_category'].value_counts()
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Color palette
    colors = sns.color_palette("husl", len(category_counts))
    
    # Create horizontal stacked bar
    left = 0
    for i, (category, count) in enumerate(category_counts.items()):
        percentage = (count / len(detailed_df)) * 100
        ax.barh(0, count, left=left, height=0.5, 
                label=f'{category} ({percentage:.1f}%)', 
                color=colors[i], edgecolor='white', linewidth=2)
        
        # Add count label in the middle of the bar
        if percentage > 5:  # Only show label if segment is large enough
            ax.text(left + count/2, 0, f'{count}', 
                   ha='center', va='center', fontweight='bold', fontsize=9)
        
        left += count
    
    ax.set_xlim(0, len(detailed_df))
    ax.set_ylim(-0.5, 0.5)
    ax.set_xlabel('Number of Reads', fontsize=12, fontweight='bold')
    ax.set_title(f'Gene Category Contribution - Sample {sample_id}', 
                fontsize=14, fontweight='bold', pad=20)
    ax.set_yticks([])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=True, fontsize=10)
    ax.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/{sample_id}_gene_category_contribution.png', 
                bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"✓ Created: {sample_id}_gene_category_contribution.png")


def plot_species_contribution(species_df, detailed_df, sample_id, output_dir, top_n=15):
    """Plot stacked bar chart of top species contributions to methane genes"""
    
    # Get species from detailed genes (from DIAMOND annotations)
    diamond_species = detailed_df['species'].value_counts().head(top_n)
    
    # Get species from SingleM (actual taxonomy)
    singlem_species = species_df.nlargest(top_n, 'num_hits')
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Plot 1: DIAMOND annotations (gene-associated species)
    colors1 = sns.color_palette("Set3", len(diamond_species))
    left = 0
    for i, (species, count) in enumerate(diamond_species.items()):
        percentage = (count / len(detailed_df)) * 100
        ax1.barh(0, count, left=left, height=0.6,
                label=f'{species[:30]}... ({percentage:.1f}%)', 
                color=colors1[i], edgecolor='white', linewidth=1.5)
        if percentage > 3:
            ax1.text(left + count/2, 0, f'{count}', 
                    ha='center', va='center', fontsize=8, fontweight='bold')
        left += count
    
    ax1.set_xlim(0, len(detailed_df))
    ax1.set_ylim(-0.5, 0.5)
    ax1.set_xlabel('Number of Gene Hits', fontsize=11, fontweight='bold')
    ax1.set_title('DIAMOND Gene Annotations\n(Species from methane gene hits)', 
                 fontsize=12, fontweight='bold', pad=15)
    ax1.set_yticks([])
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), 
              ncol=2, frameon=True, fontsize=8)
    ax1.grid(axis='x', alpha=0.3)
    
    # Plot 2: SingleM taxonomy (actual community composition)
    colors2 = sns.color_palette("Set2", len(singlem_species))
    left = 0
    total_hits = singlem_species['num_hits'].sum()
    for i, (idx, row) in enumerate(singlem_species.iterrows()):
        species = row['species']
        count = row['num_hits']
        percentage = (count / total_hits) * 100
        ax2.barh(0, count, left=left, height=0.6,
                label=f'{species[:30]}... ({percentage:.1f}%)', 
                color=colors2[i], edgecolor='white', linewidth=1.5)
        if percentage > 3:
            ax2.text(left + count/2, 0, f'{count}', 
                    ha='center', va='center', fontsize=8, fontweight='bold')
        left += count
    
    ax2.set_xlim(0, total_hits)
    ax2.set_ylim(-0.5, 0.5)
    ax2.set_xlabel('Number of SingleM Marker Hits', fontsize=11, fontweight='bold')
    ax2.set_title('SingleM Taxonomy\n(Actual community composition)', 
                 fontsize=12, fontweight='bold', pad=15)
    ax2.set_yticks([])
    ax2.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), 
              ncol=2, frameon=True, fontsize=8)
    ax2.grid(axis='x', alpha=0.3)
    
    plt.suptitle(f'Species Contribution Comparison - Sample {sample_id}', 
                fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/{sample_id}_species_contribution.png', 
                bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"✓ Created: {sample_id}_species_contribution.png")


def plot_procrustes_analysis(species_df, detailed_df, sample_id, output_dir):
    """
    Procrustes analysis comparing DIAMOND vs SingleM species profiles
    Note: For single sample, shows composition comparison
    """
    
    # Get species abundance from both methods
    diamond_species = detailed_df['species'].value_counts()
    singlem_species = species_df.set_index('species')['num_hits']
    
    # Get common species
    common_species = list(set(diamond_species.index) & set(singlem_species.index))
    
    if len(common_species) < 3:
        print(f"⚠ Warning: Only {len(common_species)} common species found. Skipping Procrustes analysis.")
        print("   Procrustes analysis requires at least 3 common species for meaningful comparison.")
        return
    
    # Create abundance matrix for common species
    diamond_abund = np.array([diamond_species.get(sp, 0) for sp in common_species])
    singlem_abund = np.array([singlem_species.get(sp, 0) for sp in common_species])
    
    # Normalize
    diamond_abund = diamond_abund / diamond_abund.sum()
    singlem_abund = singlem_abund / singlem_abund.sum()
    
    # Create 2D coordinates for visualization (using log abundance and rank)
    diamond_coords = np.column_stack([
        np.log10(diamond_abund + 1e-10),
        np.arange(len(common_species))
    ])
    
    singlem_coords = np.column_stack([
        np.log10(singlem_abund + 1e-10),
        np.arange(len(common_species))
    ])
    
    # Perform Procrustes analysis
    mtx1, mtx2, disparity = procrustes(diamond_coords, singlem_coords)
    
    # Calculate correlation
    correlation = np.corrcoef(diamond_abund, singlem_abund)[0, 1]
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: Procrustes superimposition
    ax1.scatter(mtx1[:, 0], mtx1[:, 1], c='blue', s=100, alpha=0.6, 
               label='DIAMOND', edgecolors='black', linewidth=1)
    ax1.scatter(mtx2[:, 0], mtx2[:, 1], c='red', s=100, alpha=0.6, 
               label='SingleM', edgecolors='black', linewidth=1)
    
    # Draw arrows connecting corresponding points
    for i in range(len(mtx1)):
        ax1.arrow(mtx1[i, 0], mtx1[i, 1], 
                 mtx2[i, 0] - mtx1[i, 0], mtx2[i, 1] - mtx1[i, 1],
                 alpha=0.3, head_width=0.05, head_length=0.05, 
                 fc='gray', ec='gray', linestyle='--', linewidth=0.5)
    
    ax1.set_xlabel('Dimension 1 (log abundance)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Dimension 2 (rank)', fontsize=11, fontweight='bold')
    ax1.set_title(f'Procrustes Analysis\nDisparity: {disparity:.4f}', 
                 fontsize=12, fontweight='bold')
    ax1.legend(fontsize=10, frameon=True)
    ax1.grid(alpha=0.3)
    
    # Add text box with statistics
    textstr = f'M² statistic: {disparity:.4f}\nCorrelation: {correlation:.3f}\nCommon species: {len(common_species)}'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax1.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=9,
            verticalalignment='top', bbox=props)
    
    # Plot 2: Abundance correlation
    ax2.scatter(diamond_abund, singlem_abund, s=100, alpha=0.6, 
               c=np.arange(len(common_species)), cmap='viridis',
               edgecolors='black', linewidth=1)
    
    # Add 1:1 line
    max_val = max(diamond_abund.max(), singlem_abund.max())
    ax2.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, label='1:1 line')
    
    # Add regression line
    z = np.polyfit(diamond_abund, singlem_abund, 1)
    p = np.poly1d(z)
    x_line = np.linspace(0, max_val, 100)
    ax2.plot(x_line, p(x_line), 'r-', alpha=0.7, linewidth=2, 
            label=f'Fit: y={z[0]:.2f}x+{z[1]:.2e}')
    
    ax2.set_xlabel('DIAMOND Relative Abundance', fontsize=11, fontweight='bold')
    ax2.set_ylabel('SingleM Relative Abundance', fontsize=11, fontweight='bold')
    ax2.set_title(f'Abundance Correlation\nR² = {correlation**2:.3f}', 
                 fontsize=12, fontweight='bold')
    ax2.legend(fontsize=9, frameon=True)
    ax2.grid(alpha=0.3)
    
    plt.suptitle(f'Method Comparison: DIAMOND vs SingleM - Sample {sample_id}', 
                fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/{sample_id}_procrustes_analysis.png', 
                bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"✓ Created: {sample_id}_procrustes_analysis.png")
    print(f"  Procrustes disparity (M²): {disparity:.4f} (lower is better)")
    print(f"  Abundance correlation: {correlation:.3f}")
    print(f"  Common species analyzed: {len(common_species)}")


def plot_gene_rpkm_heatmap(rpkm_df, sample_id, output_dir, top_n=20):
    """Plot heatmap of top genes by RPKM"""
    
    top_genes = rpkm_df.nlargest(top_n, 'rpkm')
    
    # Create data for heatmap
    data = top_genes[['read_count', 'length', 'pident', 'rpkm']].T
    data.columns = [idx[:40] + '...' if len(idx) > 40 else idx 
                    for idx in top_genes.index]
    
    # Normalize each row for better visualization
    data_norm = data.div(data.max(axis=1), axis=0)
    
    fig, ax = plt.subplots(figsize=(14, 6))
    sns.heatmap(data_norm, annot=False, cmap='YlOrRd', 
                cbar_kws={'label': 'Normalized Value'}, 
                linewidths=0.5, linecolor='white', ax=ax)
    
    ax.set_xlabel('Genes', fontsize=11, fontweight='bold')
    ax.set_ylabel('Metrics', fontsize=11, fontweight='bold')
    ax.set_title(f'Top {top_n} Genes by RPKM - Sample {sample_id}', 
                fontsize=12, fontweight='bold', pad=15)
    ax.set_yticklabels(['Read Count', 'Avg Length', 'Avg Identity%', 'RPKM'], 
                       rotation=0)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/{sample_id}_gene_rpkm_heatmap.png', 
                bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"✓ Created: {sample_id}_gene_rpkm_heatmap.png")


def main():
    if len(sys.argv) < 2:
        print("Usage: python visualize_results.py <sample_id> [output_dir]")
        print("\nExample:")
        print("  python visualize_results.py 53394")
        print("  python visualize_results.py 53394 Results/quick_search")
        sys.exit(1)
    
    sample_id = sys.argv[1]
    
    # Determine output directory
    if len(sys.argv) > 2:
        output_dir = sys.argv[2]
    else:
        # Auto-detect from quick_search results
        output_dir = "Results/quick_search"
    
    # Check if output directory exists
    if not os.path.exists(output_dir):
        print(f"ERROR: Output directory not found: {output_dir}")
        sys.exit(1)
    
    # Load data
    detailed_file = f"{output_dir}/{sample_id}_methane_genes_detailed.csv"
    species_file = f"{output_dir}/{sample_id}_species_abundance.csv"
    rpkm_file = f"{output_dir}/{sample_id}_gene_rpkm.csv"
    
    print("\n" + "="*80)
    print(f"VISUALIZING RESULTS - Sample {sample_id}")
    print("="*80)
    print(f"Loading data from: {output_dir}")
    
    # Load dataframes
    if not os.path.exists(detailed_file):
        print(f"ERROR: {detailed_file} not found")
        sys.exit(1)
    if not os.path.exists(species_file):
        print(f"ERROR: {species_file} not found")
        sys.exit(1)
    if not os.path.exists(rpkm_file):
        print(f"ERROR: {rpkm_file} not found")
        sys.exit(1)
    
    detailed_df = pd.read_csv(detailed_file)
    species_df = pd.read_csv(species_file)
    rpkm_df = pd.read_csv(rpkm_file, index_col=0)
    
    print(f"✓ Loaded {len(detailed_df)} gene records")
    print(f"✓ Loaded {len(species_df)} species records")
    print(f"✓ Loaded {len(rpkm_df)} gene RPKM values\n")
    
    print("Generating plots...")
    print("-" * 80)
    
    # Generate all plots
    plot_gene_category_contribution(detailed_df, sample_id, output_dir)
    plot_species_contribution(species_df, detailed_df, sample_id, output_dir)
    plot_procrustes_analysis(species_df, detailed_df, sample_id, output_dir)
    plot_gene_rpkm_heatmap(rpkm_df, sample_id, output_dir)
    
    print("\n" + "="*80)
    print("VISUALIZATION COMPLETE")
    print("="*80)
    print(f"\nAll plots saved to: {output_dir}/")
    print()


if __name__ == "__main__":
    main()
