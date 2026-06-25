# =============================================================================
# Supplemental Heterozygosity Plots
# =============================================================================
# Plots windowed heterozygosity (H) for control vs rotenone samples at
# generation 60, per cage, for chromosomes 2L, 2R, and 3R. Significant SNP
# positions are marked with vertical lines.
#
# Author: Leah Darwin, Camille Brown
#
# Output files:
#   output/hetero_2L.pdf : Heterozygosity tracks for chromosome 2L, all cages
#   output/hetero_2R.pdf : Heterozygosity tracks for chromosome 2R, all cages
#   output/hetero_3R.pdf : Heterozygosity tracks for chromosome 3R, all cages
# =============================================================================

# ============================================================================
# LOAD PACKAGES
# ============================================================================

library(dplyr)     # For data manipulation
library(ggplot2)   # For creating plots
library(patchwork) # For combining multiple plots

# ============================================================================
# DATA IMPORT AND PREPROCESSING
# ============================================================================

# Read sample labels and create file paths for heterozygosity data
labels = read.csv("data/pca/sync_labs.csv") %>%
  mutate(pop_file = paste("data/hetero/f", row_number(), "_H.tsv",sep = "")) %>%
  # Filter to generation 60 only
  filter(gen == 60)

# Add index column for matching data frames
labels$idx = row_number(labels)

# Read all heterozygosity files into list of data frames
dfs = lapply(labels$pop_file,read.csv, sep="\t")

# ============================================================================
# FUNCTION: Plot Heterozygosity
# ============================================================================

#' Plot heterozygosity (H) for a specific chromosome
#'
#' @param df_C Data frame with heterozygosity for control samples
#' @param df_R Data frame with heterozygosity for rotenone samples
#' @param chr Chromosome name (e.g., "2R", "3R")
#' @param title Plot title (optional)
#' @return ggplot object
plot_H = function(df_C, df_R, chr, title){
  
  # Filter to specified chromosome
  df_C = df_C %>%
    filter(chrom == chr)
  df_R = df_R %>%
    filter(chrom == chr)
  
  # Create heterozygosity plot
  # Black line = control, red line = rotenone
  p = ggplot() +
    geom_line(data = df_C, aes(x = (window_start+window_end)/2000000, y = H_mean), color="black", linewidth=0.5)+
    geom_line(data = df_R, aes(x = (window_start+window_end)/2000000, y = H_mean), color = "red", linewidth=0.5)+
    theme_classic() +
    labs(x ="", y="")+
    theme(
      plot.margin = margin(1, 1, 1, 1)
    ) +
    scale_y_continuous(
      breaks = c(0.1, 0.2, 0.3),
      limits = c(0, 0.3),
      expand = c(0, 0)
    )
  
  # Add vertical line at significant SNP position
  if(chr=="3R"){
    # Position on 3R: 21111506 bp
    p = p+geom_vline(xintercept = 21111506/1000000, color = "#3182bd", alpha=0.6)
    p = p+geom_vline(xintercept = 7987407/1000000, color = "#3182bd", alpha=0.6)
    p = p+geom_vline(xintercept = 30211593/1000000, color = "#3182bd", alpha=0.6)
  }else if(chr=="2R"){
    # Position on 2R: 14854642 bp
    p = p+geom_vline(xintercept = 14854642/1000000, color = "#3182bd", alpha=0.6)
    p = p+geom_vline(xintercept = 11173099/1000000, color = "#3182bd", alpha=0.6)
    p = p+geom_vline(xintercept = 7013536/1000000, color = "#3182bd", alpha=0.6)
  }else{
    # Position on 2L: 
    p = p+geom_vline(xintercept =20820035/1000000,  color = "#3182bd", alpha=0.6)
  }
  
  return(p)
  
}

# ============================================================================
# FUNCTION: Style Multiple Heterozygosity Plots
# ============================================================================

#' Style and arrange a list of heterozygosity plots
#'
#' @param plots List of ggplot objects (12 plots for 6 cages x 2 columns)
#' @return Combined patchwork plot with styled axes
style_plots = function(plots){
  
  # Interleave plots: columns 1-6 in col1, columns 7-12 in col2
  #col1 <- 1:6
  #col2 <- 7:12
  #new_order <- as.vector(rbind(col1, col2))
  #plots_list <- plots[new_order]
  
  #original order
  plots_list = plots
  
  # Arrange in 2 columns and remove all axis elements by default
  p = wrap_plots(plots_list ,ncol=3)&
    theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank()
    )
  
  # Add y-axis to left-side plots (odd indices)
  left_idx <- seq(1, 12, by = 2)
  
  ##for 4x3 arrangement
  left_idx = seq(1,12, by = 3)
  
  for (i in left_idx) {
    p[[i]] <- p[[i]] +
      theme(
        axis.text.y  = element_text(size=12),
        axis.ticks.y = element_line()
      )
  }
  
  # Add x-axis to bottom plots (last row)
  bottom_idx <- c(11,12)
  
  ##for 4x3 arrangement
  bottom_idx = c(10,11,12)
  
  for (i in bottom_idx) {
    p[[i]] <- p[[i]] +
      theme(
        axis.text.x  = element_text(size=12),
        axis.ticks.x = element_line()
      )
  }
  
  return(p)
  
}

# ============================================================================
# CREATE HETEROZYGOSITY PLOTS FOR ALL CAGES
# ============================================================================

cages = unique(labels$cage)

# Plot heterozygosity for chromosome 2L across all cages
# by() applies function to each subset defined by labels$cage
p_2L = by(labels, labels$cage, function(x) {
  plot_H(dfs[[x$idx[x$treatment == "C"]]],
         dfs[[x$idx[x$treatment == "R"]]], "2L", x$cage)
})

# Style and arrange 2R plots
p4 = style_plots(p_2L)
ggsave("output/hetero_2L.pdf", plot=p4, width = 8, height = 4)

# Plot heterozygosity for chromosome 2R across all cages
# by() applies function to each subset defined by labels$cage
p_2R = by(labels, labels$cage, function(x) {
  plot_H(dfs[[x$idx[x$treatment == "C"]]],
         dfs[[x$idx[x$treatment == "R"]]], "2R", x$cage)
})

# Style and arrange 2R plots
p2 = style_plots(p_2R)
ggsave("output/hetero_2R.pdf", plot=p2, width = 8, height = 4)

# Plot heterozygosity for chromosome 3R across all cages
p_3R = by(labels, labels$cage, function(x) {
  plot_H(dfs[[x$idx[x$treatment == "C"]]],
         dfs[[x$idx[x$treatment == "R"]]], "3R")
})

# Style and arrange 3R plots
p3 = style_plots(p_3R)
ggsave("output/hetero_3R.pdf", plot=p3, width = 8, height = 4)
