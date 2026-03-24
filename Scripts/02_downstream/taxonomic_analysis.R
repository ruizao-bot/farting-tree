#!/usr/bin/env Rscript

# Minimal downstream analysis: heatmap, PCoA, PERMANOVA

# Load required packages
required_packages <- c("tidyverse", "vegan", "ggplot2", "RColorBrewer", "pheatmap")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

setwd("/Users/jiayi/Desktop/metagenomic_pipeline/QIIME")

# Import data
feature_table <- read.table(
  "Results/denoise_mode/exported-table/feature-table.tsv",
  header = TRUE,
  sep = "\t",
  skip = 1,
  row.names = 1,
  comment.char = ""
)
feature_table <- feature_table[, !(colnames(feature_table) %in% c("A1", "A2")), drop = FALSE]

taxonomy <- read.table(
  "Results/denoise_mode/exported-taxonomy/taxonomy.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  quote = ""
)

metadata <- read.table(
  "Data/metadata/metadata.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  comment.char = ""
)
# Do not subset metadata by feature table columns. 
# Keep full metadata to support PCoA plotting of all samples in the ordination file.
# metadata <- metadata[colnames(feature_table), ]

# Parse taxonomy
parse_taxonomy <- function(tax_string) {
  tax_clean <- gsub("\\s*\\([^\\)]+\\)", "", tax_string)
  levels <- strsplit(tax_clean, ";")[[1]]
  levels <- trimws(levels)
  tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  result <- rep(NA, 7)
  names(result) <- tax_levels
  result[1:min(length(levels), 7)] <- levels[1:min(length(levels), 7)]
  return(result)
}
tax_table <- t(sapply(taxonomy$Taxon, parse_taxonomy))
rownames(tax_table) <- rownames(taxonomy)

# Remove Cyanobacteria and Chloroplast
cat("Filtering out Cyanobacteria and Chloroplast...\n")
cat("  Features before filtering:", nrow(feature_table), "\n")

# Identify sequences to remove
remove_idx <- (
  grepl("Cyanobacteria", tax_table[, "Phylum"], ignore.case = TRUE) |
  grepl("Cyanobacteriia", tax_table[, "Class"], ignore.case = TRUE) |
  grepl("Chloroplast", tax_table[, "Order"], ignore.case = TRUE) |
  grepl("Chloroplast", tax_table[, "Family"], ignore.case = TRUE) |
  grepl("Chloroplast", tax_table[, "Genus"], ignore.case = TRUE)
)

cat("  Cyanobacteria/Chloroplast features found:", sum(remove_idx), "\n")

# Remove from both feature table and taxonomy table
feature_table <- feature_table[!remove_idx, ]
tax_table <- tax_table[!remove_idx, ]

cat("  Features after filtering:", nrow(feature_table), "\n\n")

# Filter datasets by domain
cat("Filtering datasets by domain...\n")

# Get Kingdom/Domain column
kingdoms <- tax_table[, "Kingdom"]

# Dataset 1: Bacteria only
bacteria_idx <- grepl("Bacteria", kingdoms, ignore.case = TRUE)
feature_table_bacteria <- feature_table[bacteria_idx, ]
tax_table_bacteria <- tax_table[bacteria_idx, ]

# Dataset 2: Bacteria and Archaea
bacteria_archaea_idx <- grepl("Bacteria|Archaea", kingdoms, ignore.case = TRUE)
feature_table_bacteria_archaea <- feature_table[bacteria_archaea_idx, ]
tax_table_bacteria_archaea <- tax_table[bacteria_archaea_idx, ]

cat("Dataset summary:\n")
cat("  Total features:", nrow(feature_table), "\n")
cat("  Bacteria only:", nrow(feature_table_bacteria), "\n")
cat("  Bacteria + Archaea:", nrow(feature_table_bacteria_archaea), "\n\n")

# Create output directories for each dataset
dir.create("Results/denoise_mode/bacteria_only", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/denoise_mode/bacteria_archaea", showWarnings = FALSE, recursive = TRUE)

# Function to run analysis on a dataset
run_analysis <- function(ft, tt, meta, output_prefix, dataset_name) {
  cat("\n=== Processing", dataset_name, "===\n")
  
  # Relative abundance
  rel_abundance <- sweep(ft, 2, colSums(ft), "/") * 100
  
  # Aggregate by genus
  aggregate_by_level <- function(feat_table, tax_table, level) {
    # Ensure feature tables and taxonomy tables are aligned
    common_features <- intersect(rownames(feat_table), rownames(tax_table))
    feat_table <- feat_table[common_features, , drop = FALSE]
    tax_table <- tax_table[common_features, , drop = FALSE]
    
    taxa <- tax_table[, level]
    taxa[is.na(taxa)] <- "Unclassified"
    agg_table <- aggregate(feat_table, by = list(Taxa = taxa), FUN = sum)
    rownames(agg_table) <- agg_table$Taxa
    agg_table$Taxa <- NULL
    return(as.data.frame(agg_table))
  }
  genus_abund <- aggregate_by_level(rel_abundance, tt, "Genus")
  
  # Heatmap
  plot_heatmap <- function(abund_table, metadata, top_n = 30, cluster_samples = TRUE, group_var = "sample_type") {
    mean_abund <- rowMeans(abund_table)
    top_taxa <- names(sort(mean_abund, decreasing = TRUE)[1:top_n])
    heat_data <- abund_table[top_taxa, ]
    heat_data_log <- log10(heat_data + 0.01)
    annotation_col <- metadata[, group_var, drop = FALSE]
    
    # Custom colors for sample_type annotation
    sample_types <- unique(metadata[, group_var])
    annotation_colors <- list()
    if (length(sample_types) == 3) {
      annotation_colors[[group_var]] <- c("green3", "purple3", "gold")
      names(annotation_colors[[group_var]]) <- sort(sample_types)
    }
    
    pheatmap::pheatmap(
      heat_data_log,
      cluster_rows = TRUE,
      cluster_cols = cluster_samples,
      annotation_col = annotation_col,
      annotation_colors = annotation_colors,
      color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
      fontsize = 8,
      fontsize_row = 7,
      fontsize_col = 8,
      main = paste("Top", top_n, "Taxa Heatmap -", dataset_name, "\n[log10(Relative Abundance %)]"),
      filename = paste0(output_prefix, "_taxa_heatmap.pdf"),
      width = 10,
      height = 12
    )
  }
  plot_heatmap(genus_abund, meta, top_n = 30, group_var = "sample_type")
  cat("Heatmap saved\n")
  
  # Alpha diversity
  alpha_div <- data.frame(
    Sample = colnames(ft),
    ASV_abundance = colSums(ft > 0),
    Shannon = vegan::diversity(t(ft), index = "shannon"),
    Simpson = vegan::diversity(t(ft), index = "simpson")
  )
  alpha_div <- cbind(alpha_div, meta[alpha_div$Sample, ])
  
  alpha_long <- alpha_div %>%
    pivot_longer(cols = c(ASV_abundance, Shannon, Simpson), names_to = "Metric", values_to = "Value")
  
  p_alpha <- ggplot(alpha_long, aes(x = sample_type, y = Value, fill = sample_type)) +
    geom_boxplot() +
    facet_wrap(~ Metric, scales = "free_y") +
    theme_bw() +
    labs(title = paste("Alpha Diversity -", dataset_name), x = "Sample Type", y = "Value") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p_alpha)
  ggsave(paste0(output_prefix, "_alpha_diversity_boxplot.pdf"), p_alpha, width = 10, height = 4)
  cat("Alpha diversity boxplot saved\n")

  # --- PCoA and PERMANOVA (Global Unweighted UniFrac) ---
  # Note: This uses the global ordination/distance matrix, not specific to the subset (Bacteria/Archaea)
  # unless those files are updated. We run this block to visualize the provided QIIME2 results.
  
  cat("\n=== PCoA with Unweighted UniFrac (QIIME2) ===\n")
  
  unweighted_pcoa_file <- "Results/denoise_mode/diversity/exported-unweighted_unifrac-pcoa/ordination.txt"
  unweighted_dist_file <- "Results/denoise_mode/diversity/unweighted_unifrac_exported/distance-matrix.tsv"
  
  if (file.exists(unweighted_pcoa_file)) {
    # Robust Parsing of QIIME2 ordination file
    lines <- readLines(unweighted_pcoa_file)
    prop_idx_all <- grep("^Proportion", lines)
    site_idx_all <- grep("^Site", lines)
    
    if (length(prop_idx_all) > 0 && length(site_idx_all) > 0) {
      prop_idx <- prop_idx_all[1]
      site_idx <- site_idx_all[1]
      
      # Parse Proportion Explained
      prop_line_idx <- prop_idx + 1
      while(prop_line_idx <= length(lines) && nchar(trimws(lines[prop_line_idx])) == 0) prop_line_idx <- prop_line_idx + 1
      prop_tokens <- unlist(strsplit(lines[prop_line_idx], "[\t ]+"))
      prop_tokens <- prop_tokens[prop_tokens != ""]
      percent_var <- round(as.numeric(prop_tokens) * 100, 2)
      
      # Parse Site Coordinates
      coord_start <- site_idx + 1
      # Find end of coordinate block (first blank line after start)
      blank_lines <- which(nchar(trimws(lines)) == 0)
      coord_end <- if (any(blank_lines > coord_start)) min(blank_lines[blank_lines > coord_start]) - 1 else length(lines)
      
      coord_lines <- lines[coord_start:coord_end]
      coord_lines <- coord_lines[nzchar(trimws(coord_lines))]
      
      if (length(coord_lines) > 0) {
        coord_tab <- read.table(text = paste(coord_lines, collapse = "\n"), header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
        # Expecting SampleID, PC1, PC2, ...
        pcoa_scores <- data.frame(
          Sample = coord_tab[, 1],
          PC1 = as.numeric(coord_tab[, 2]),
          PC2 = as.numeric(coord_tab[, 3]),
          stringsAsFactors = FALSE
        )
        
        # Merge with Metadata (Left Join to keep all ordination samples)
        # We use the full 'metadata' object (which is now not subsetted)
        meta_df <- data.frame(Sample = rownames(metadata), metadata, stringsAsFactors = FALSE)
        plot_data <- merge(pcoa_scores, meta_df, by = "Sample", all.x = TRUE)
        
        # Handle missing metadata for plotting
        if (!"sample_type" %in% colnames(plot_data)) plot_data$sample_type <- NA
        if (!"species" %in% colnames(plot_data)) plot_data$species <- NA
        
        plot_data$sample_type[is.na(plot_data$sample_type)] <- "Unknown"
        plot_data$species[is.na(plot_data$species)] <- "Unknown"
        
        # Plot PCoA
        p_pcoa <- ggplot(plot_data, aes(x = PC1, y = PC2, color = sample_type, shape = species)) +
          geom_point(size = 4, alpha = 0.8) +
          theme_bw() +
          labs(
            title = paste("PCoA - Unweighted UniFrac (QIIME2)"),
            subtitle = paste("Dataset context:", dataset_name),
            x = paste0("PC1 (", percent_var[1], "%)"),
            y = paste0("PC2 (", percent_var[2], "%)"),
            color = "Sample Type",
            shape = "Species"
          ) +
          theme(plot.title = element_text(hjust = 0.5))
        
        ggsave(paste0(output_prefix, "_pcoa_unweighted_unifrac.pdf"), p_pcoa, width = 10, height = 7)
        cat("PCoA plot saved to:", paste0(output_prefix, "_pcoa_unweighted_unifrac.pdf\n"))
        
        # Save Scores
        write.table(plot_data, paste0(output_prefix, "_pcoa_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
      }
    }
  } else {
    cat("Ordination file not found:", unweighted_pcoa_file, "\n")
  }
  
  # --- PERMANOVA ---
  cat("\n=== PERMANOVA Analysis ===\n")
  if (file.exists(unweighted_dist_file)) {
    dist_mat <- read.table(unweighted_dist_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    
    # Intersect samples between distance matrix and metadata
    # We need valid metadata for PERMANOVA factors
    common_samples <- intersect(rownames(dist_mat), rownames(metadata))
    
    # Filter out samples with NA in key columns if necessary, or just use common ones
    # Here we just use samples present in both.
    
    if (length(common_samples) > 0) {
      dist_subset <- as.dist(dist_mat[common_samples, common_samples])
      meta_subset <- metadata[common_samples, , drop = FALSE]
      
      sink(paste0(output_prefix, "_PERMANOVA_results.txt"))
      cat("PERMANOVA Results (Unweighted UniFrac)\n")
      cat("Samples used:", length(common_samples), "\n\n")
      
      # Test Species
      if ("species" %in% colnames(meta_subset) && length(unique(meta_subset$species)) > 1) {
        cat("--- Effect of Species ---\n")
        tryCatch({
          print(adonis2(dist_subset ~ species, data = meta_subset, permutations = 999))
        }, error = function(e) cat("Error running PERMANOVA on species:", e$message, "\n"))
        cat("\n")
      }
      
      # Test Sample Type
      if ("sample_type" %in% colnames(meta_subset) && length(unique(meta_subset$sample_type)) > 1) {
        cat("--- Effect of Sample Type ---\n")
        tryCatch({
          print(adonis2(dist_subset ~ sample_type, data = meta_subset, permutations = 999))
        }, error = function(e) cat("Error running PERMANOVA on sample_type:", e$message, "\n"))
        cat("\n")
      }
      
      sink()
      cat("PERMANOVA results saved to:", paste0(output_prefix, "_PERMANOVA_results.txt\n"))
    } else {
      cat("No common samples between distance matrix and metadata for PERMANOVA.\n")
    }
  } else {
    cat("Distance matrix file not found:", unweighted_dist_file, "\n")
  }
  
  cat("\n", dataset_name, "analysis complete\n")
}

# Run analysis for both datasets
run_analysis(feature_table_bacteria, tax_table_bacteria, metadata, 
             "Results/denoise_mode/bacteria_only/bacteria_only", "Bacteria Only")
run_analysis(feature_table_bacteria_archaea, tax_table_bacteria_archaea, metadata, 
             "Results/denoise_mode/bacteria_archaea/bacteria_archaea", "Bacteria + Archaea")
