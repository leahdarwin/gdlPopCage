# =============================================================================
# RNA-seq MDS Analysis
# =============================================================================
# Performs multidimensional scaling (MDS) on VST-transformed RNA-seq count data
# using DESeq2, then plots samples by treatment and population type for all
# samples combined and for individual cage populations.
#
# Author: Leah Darwin, Yevgeniy Raynes
#
# Output files:
#   output/mdsplot_all.png              : MDS plot for all samples combined
#   output/mdsplot_1A3.png              : MDS plot for cage population 1A3
#   output/mdsplot_2A2.png              : MDS plot for cage population 2A2
#   output/mdsplot_2A3.png              : MDS plot for cage population 2A3
#   output/mdsplot_population_panel.png : Combined panel of per-population MDS
#                                         plots (1A3, 2A2, 2A3)
# =============================================================================

# Setup
knitr::opts_chunk$set(echo = TRUE)

# Parameters
set.seed(1234)
FDR.Thresh <- 0.05
LFC.Thresh <- 0.0

# Set to TRUE to exclude the 2A2C.R.1 sample from all plots
# (this sample clusters with untreated samples instead of treated samples)
EXCLUDE_2A2C_R1 <- TRUE

# Load packages
library(tidyverse)
library(DESeq2)
library(ggrepel)
library(dplyr)
library(ggforce)
library(patchwork)

# ============================================================================
# FUNCTION DEFINITIONS
# ============================================================================

#' Perform MDS analysis on count data
#'
#' @param count_data Matrix of count data
#' @param col_data Data frame with sample metadata
#' @return List containing mds data frame and dimension percentages
perform_mds_analysis <- function(count_data, col_data) {
  # Create DESeq2 dataset and perform VST
  dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData = col_data,
    design = ~ Group
  )
  vst <- vst(dds, blind = TRUE)

  # Calculate distances and perform MDS
  sampleDists <- dist(t(assay(vst)))
  sampleDistMatrix <- as.matrix(sampleDists)
  mds_dist <- cmdscale(sampleDistMatrix, k = 3, eig = TRUE)

  # Combine MDS results with sample metadata
  mds <- cbind(as.data.frame(colData(vst)), mds_dist$points)

  # Calculate percentage variation explained by each dimension
  Dim1 <- signif(mds_dist$eig[1] / sum(mds_dist$eig) * 100, digits = 4)
  Dim2 <- signif(mds_dist$eig[2] / sum(mds_dist$eig) * 100, digits = 4)

  # Add name column
  mds$name <- paste0(mds$Group, mds$Replicate)

  return(list(
    mds = mds,
    Dim1 = Dim1,
    Dim2 = Dim2
  ))
}

#' Create MDS plot
#'
#' @param mds_result List returned from perform_mds_analysis()
#' @param label_col Column name to use for text labels (quoted or as expression)
#' @return ggplot object
create_mds_plot <- function(mds_result, label_col = "Name", show_labels = TRUE) {
  mds <- mds_result$mds
  Dim1 <- mds_result$Dim1
  Dim2 <- mds_result$Dim2

  # Create label aesthetic based on label_col
  if (is.character(label_col) && length(label_col) == 1 && label_col %in% colnames(mds)) {
    label_aes <- aes_string(label = label_col)
  } else {
    # label_col is a vector of labels - use directly
    mds$plot_label <- label_col
    label_aes <- aes(label = plot_label)
  }

  plot <- ggplot(mds, aes(x = `1`, y = `2`, color = Treatment, shape = Type)) +
    geom_point(size = 4, stroke = 1.7) +
   # coord_fixed() +
    labs(
      x = paste0("Dimension 1 (", Dim1, "%)"),
      y = paste0("Dimension 2 (", Dim2, "%)"),
      color = "Population genotype",
      shape = "Treatment exposure"
    ) +
    scale_color_manual(
      values = c(
        "R" = "red",
        "C" = "black"
      )
    ) +
    scale_fill_manual(
      values = c(
        "R" = "red",
        "C" = "black"
      )
      )+
    theme_minimal() +
    theme(
      legend.position = "none",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 14)
    )

  # Add labels if requested
  if (show_labels) {
    plot <- plot + geom_text_repel(
      label_aes,
      size = 4,
      max.overlaps = Inf,
      show.legend = FALSE
    )
  }

  return(plot)
}


#' Filter count and coldata by sample name pattern
#'
#' @param cts Count matrix
#' @param coldata Sample metadata
#' @param pattern Substring pattern to match at start of column names
#' @param nchar_match Number of characters to match from start
#' @return List with filtered cts and coldata
filter_samples <- function(cts, coldata, pattern, nchar_match = NULL) {
  if (is.null(nchar_match)) {
    nchar_match <- nchar(pattern)
  }

  cts_filtered <- cts[, substring(colnames(cts), 1, nchar_match) == pattern]
  coldata_filtered <- coldata[colnames(cts_filtered), ]

  return(list(
    cts = cts_filtered,
    coldata = coldata_filtered
  ))
}

#' Run complete MDS analysis and save plot
#'
#' @param count_data Matrix of count data
#' @param col_data Data frame with sample metadata
#' @param output_file Path for output plot file
#' @param label_col Column name for plot labels
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Plot resolution
analyze_and_save_mds <- function(count_data, col_data, output_file,
                                  label_col = "Name", show_labels = TRUE,
                                  width = 6, height = 4, dpi = 300) {
  # Perform MDS analysis
  mds_result <- perform_mds_analysis(count_data, col_data)

  # Create plot
  plot <- create_mds_plot(mds_result, label_col, show_labels)

  # Display and save plot
  print(plot)
  ggsave(output_file, plot = plot, width = width, height = height, dpi = dpi)

  return(plot)
}

# ============================================================================
# MAIN ANALYSIS
# ============================================================================

# Read count tables/group info into R
cts <- as.matrix(read.csv("data/GDL_count_table.csv", row.names = 1, check.names = FALSE))
colnames(cts)
gene_universe <- row.names(cts[rowSums(cts) > 0, ])
coldata <- read.csv("data/GDL_group_data.csv", row.names = 1) 

# Convert columns to factors
coldata$Name <- as.factor(coldata$Name)
coldata$Type <- as.factor(coldata$Type)
coldata$Treatment <- as.factor(coldata$Treatment)
coldata$Group <- as.factor(coldata$Group)
coldata$History <- as.factor(coldata$History)
coldata$TypeTreat <- as.factor(paste(coldata$Type, coldata$Treatment,sep=""))

#

# Optionally exclude the 2A2C.R.1 sample
if (EXCLUDE_2A2C_R1) {
  cat("Excluding outlier samples from all analyses...\n")
  samples_to_keep <- colnames(cts) != "2A2.R.C.2" & colnames(cts) != "2A2.C.R.1" & colnames(cts) != "2A3.R.C.1"
  cts <- cts[, samples_to_keep]
  coldata <- coldata[colnames(cts), ]
}

# ============================================================================
# MDS Analysis - All Samples
# ============================================================================

cat("Analyzing all samples...\n")
analyze_and_save_mds(
  count_data = cts,
  col_data = coldata,
  output_file = "output/mdsplot_all.png",
  label_col = "Name"
)

# ============================================================================
# MDS Analysis by History
# ============================================================================

# 1A3
cat("\nAnalyzing 1A3 samples...\n")
filtered_1A3 <- filter_samples(cts, coldata, pattern = "1A3", nchar_match = 3)
mds_result_1A3 <- perform_mds_analysis(filtered_1A3$cts, filtered_1A3$coldata)
plot_1A3 <- create_mds_plot(mds_result_1A3, show_labels = FALSE)
print(plot_1A3)
ggsave("output/mdsplot_1A3.png", plot = plot_1A3, width = 6, height = 4, dpi = 300)

# 2A2
cat("\nAnalyzing 2A2 samples...\n")
filtered_2A2 <- filter_samples(cts, coldata, pattern = "2A2", nchar_match = 3)
mds_result_2A2 <- perform_mds_analysis(filtered_2A2$cts, filtered_2A2$coldata)
plot_2A2 <- create_mds_plot(mds_result_2A2, show_labels = FALSE)
print(plot_2A2)
ggsave("output/mdsplot_2A2.png", plot = plot_2A2, width = 6, height = 4, dpi = 300)



# 2A3
cat("\nAnalyzing 2A3 samples...\n")
filtered_2A3 <- filter_samples(cts, coldata, pattern = "2A3", nchar_match = 3)
mds_result_2A3 <- perform_mds_analysis(filtered_2A3$cts, filtered_2A3$coldata)
plot_2A3 <- create_mds_plot(mds_result_2A3, show_labels = FALSE)
print(plot_2A3)
ggsave("output/mdsplot_2A3.png", plot = plot_2A3, width = 6, height = 4, dpi = 300)


final = plot_1A3 + plot_spacer()  + plot_2A2 + plot_2A3 + plot_layout(nrow=2) &
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')',
                  theme = theme(plot.margin = margin(-1, -1, -1, -1)))
ggsave("output/mdsplot_population_panel.png", final, width = 10, height = 8, dpi = 300)
