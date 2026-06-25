# =============================================================================
# SNP Cluster and Haplotype Frequency Analysis
# =============================================================================
# Clusters significant SNPs on each chromosome using k-means, subsets haplotype
# frequency data to those cluster regions, and plots founder frequency
# trajectories and inter-cage Jaccard similarity heatmaps.
#
# Author: Leah Darwin
#
# Output files:
#   data/clust_bounds.csv          : Cluster genomic boundaries and founders that
#                                    increased in frequency (>0.1 change F20->F60)
#   output/jaccard_[treatment].pdf : Jaccard similarity heatmap between cage
#                                    populations based on shared increased founders
# =============================================================================

# ============================================================================
# LOAD PACKAGES
# ============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)

# ============================================================================
# DEFINE PARAMETERS
# ============================================================================

# Bonferroni-corrected significance threshold (0.05 / total tests)
bf = 0.05 / 3433901

# Small value to avoid zero frequencies in visualization
eps = 1e-6

# Generation to analyze
gen_i = 60

# Treatment to analyze ("C" = Control, "R" = Rotenone)
treatment_i = "R"

# Chromosome order for plotting
chr_order = c("2L", "2R", "3L", "3R", "4", "X", "Y")

# Color palette for founder groups
color_palette <- c(
  "Beijing"   = "#18458e",
  "Zimbabwe"  = "#4fc2ba",
  "OreR"  = "#ec2328",
  "DGRP375"= "#fbbd67",
  "w1118"    = "#984ea3"
)

# Valid line types for ggplot
valid_linetypes = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")

# All possible cages (ensures 12 facets always appear)
all_cages = c("1A1", "1A2", "1A3", "1B1", "1B2", "1B3",
              "2A1", "2A2", "2A3", "2B1", "2B2", "2B3")

# ============================================================================
# FUNCTIONS
# ============================================================================

#' Cluster SNPs on a chromosome using k-means
#'
#' @param chr Chromosome name (e.g., "2R")
#' @param k Number of clusters for k-means
#' @return Data frame of SNPs with cluster assignments (format: "CHR.cluster_num")
cluster_snps = function(chr, k) {
  chr_snps = snps %>% filter(CHR == chr)
  km = kmeans(chr_snps %>% pull(BP), k)
  chr_snps$cluster = paste(chr, km[["cluster"]], sep = ".")
  return(chr_snps)
}

#' Subset a data frame to windows overlapping cluster boundaries
#'
#' Filters genomic windows to only those falling within defined cluster regions
#' and adds a column indicating which cluster each window belongs to.
#'
#' @param df Data frame with chr and pos columns (haplotype frequency data)
#' @param bounds Data frame with cluster, CHR, start, end columns
#' @return Data frame filtered to cluster regions with cluster column added
subset_by_clusters = function(df, bounds) {
  df %>%
    rowwise() %>%
    mutate(
      cluster = {
        # Find cluster where chromosome matches and position falls within boundaries
        match_idx = which(bounds$CHR == chr & pos >= bounds$start & pos <= bounds$end)
        if (length(match_idx) > 0) bounds$cluster[match_idx[1]] else NA_character_
      }
    ) %>%
    ungroup() %>%
    # Remove rows not overlapping any cluster
    filter(!is.na(cluster))
}

#' Compute element-wise average of frequency vectors grouped by cluster
#'
#' For each cluster, parses the semicolon-separated frequency strings into
#' numeric vectors and computes the element-wise mean across all windows.
#'
#' @param df Data frame with frequencies column (semicolon-separated) and cluster column
#' @return Data frame with one row per cluster containing averaged frequency vector
average_freqs_by_cluster = function(df) {
  df %>%
    group_by(cluster,population) %>%
    summarise(
      frequencies = {
        # Parse each frequency string into a numeric vector
        freq_matrix = sapply(frequencies, function(x) {
          as.numeric(strsplit(x, ";")[[1]])
        })
        # Compute element-wise mean across columns (each column is one window)
        avg_vec = rowMeans(freq_matrix, na.rm = TRUE)
        # Collapse back to semicolon-separated string
        paste(avg_vec, collapse = ";")
      },
      .groups = "drop"
    )
}

#' Convert haplotype frequency data to long format for a specific chromosome
#'
#' @param df Data frame with haplotype frequencies
#' @param chr_i Chromosome name (e.g., "2R", "3R")
#' @return Data frame in long format with one row per founder per position
make_long = function(df){
  long = df %>%
    # Separate semicolon-delimited frequencies into multiple rows
    separate_rows(frequencies, sep = ";") %>%
    mutate(
      # Convert to numeric and add small epsilon
      freq = as.numeric(frequencies) + eps,
    ) %>%
    group_by(cluster, cage, gen) %>%
    # Assign founder names in order (matches frequency order)
    mutate(founder = founders$name[seq_len(n())]) %>%
    ungroup() %>%
    # Assign founder to group based on name prefix
    mutate(founder_group = case_when(grepl("^B",founder) ~ "Beijing",
                                     grepl("^Z", founder) ~ "Zimbabwe",
                                     grepl("^O", founder) ~ "OreR",
                                     grepl("^D", founder) ~ "DGRP375",
                                     grepl("^w", founder) ~ "w1118"))
}

#' Plot founder frequency trajectories for a given cluster
#'
#' Creates a faceted plot showing founder frequencies over generations,
#' colored by founder group with different line types for multiple founders
#' within the same group. Always displays all 12 cages.
#'
#' @param cluster_id Cluster to filter on (e.g., "2R.1", "3R.3")
#' @return ggplot object
plot_cluster_freqs = function(cluster_id) {

  # Filter and prepare plot data
  plot_data = long_df_increased %>%
    filter(cluster == cluster_id) %>%
    group_by(founder_group) %>%
    mutate(n_founders = n_distinct(founder)) %>%
    ungroup() %>%
    # Ensure cage is a factor with all levels for complete faceting
    mutate(cage = factor(cage, levels = all_cages))

  # Return empty plot if no data for this cluster
  if (nrow(plot_data) == 0) {
    return(
      ggplot() +
        facet_wrap(~factor(all_cages, levels = all_cages), nrow = 4) +
        labs(title = paste("Cluster:", cluster_id), x = "Generation", y = "Frequency") +
        theme_minimal()
    )
  }

  # Create linetype mapping: solid for single-founder groups, varied for multi-founder groups
  linetype_map = plot_data %>%
    distinct(founder_group, founder, n_founders) %>%
    group_by(founder_group) %>%
    mutate(
      lt_idx = row_number(),
      linetype = if_else(n_founders == 1, "solid", valid_linetypes[lt_idx])
    ) %>%
    ungroup()

  # Join linetype back to plot data
  plot_data = plot_data %>%
    left_join(linetype_map %>% select(founder, linetype), by = "founder")

  # Create named vector for scale_linetype_manual
  linetype_vals = linetype_map %>%
    distinct(founder, linetype) %>%
    pull(linetype, name = founder)

  # Build plot
  p = plot_data %>%
    ggplot(aes(
      x = gen,
      y = as.numeric(frequencies),
      color = founder_group,
      linetype = founder,
      group = interaction(founder, cage)
    )) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = color_palette) +
    scale_linetype_manual(values = linetype_vals) +
    facet_wrap(~cage, nrow = 4, drop = FALSE) +
    ylim(0, 1) +
    labs(
      title = paste("Cluster:", cluster_id),
      x = "Generation",
      y = "Frequency",
      color = "Founder Group",
      linetype = "Founder"
    ) +
    theme_linedraw() +
    guides(linetype = guide_legend(override.aes = list(color = "black")))

  return(p)
}

#' Calculate Jaccard similarity between two sets
#'
#' @param a First set (vector)
#' @param b Second set (vector)
#' @return Jaccard similarity (intersection / union)
jaccard = function(a, b) {
  intersection = length(intersect(a, b))
  union = length(union(a, b))
  if (union == 0) return(0)
  return(intersection / union)
}

#' Create Jaccard similarity heatmap between cages based on shared founders
#'
#' Calculates Jaccard similarity within each cluster first, then aggregates
#' (averages) across clusters for the final heatmap.
#'
#' @param df Data frame with cluster, cage, and founder columns (e.g., long_df_increased)
#' @return ggplot heatmap object
plot_jaccard_heatmap = function(df) {

  # Get unique founders per cage within each cluster
  cage_cluster_founders = df %>%
    group_by(cluster, cage) %>%
    summarise(founders = list(unique(founder)), .groups = "drop")

  clusters = unique(cage_cluster_founders$cluster)
  cages = sort(unique(df$cage))  # All cages in the data
  n_cages = length(cages)

  # Expand to include all cage-cluster combinations (empty list for missing)
  all_combinations = expand.grid(cluster = clusters, cage = cages, stringsAsFactors = FALSE)
  cage_cluster_founders = all_combinations %>%
    left_join(cage_cluster_founders, by = c("cluster", "cage")) %>%
    mutate(founders = lapply(founders, function(x) if (is.null(x)) character(0) else x))

  # Initialize list to store per-cluster Jaccard matrices
  cluster_matrices = list()

  # Calculate Jaccard similarity within each cluster for ALL cage pairs
  for (clust in clusters) {
    clust_data = cage_cluster_founders %>% filter(cluster == clust)

    # Initialize matrix for this cluster
    clust_matrix = matrix(NA, nrow = n_cages, ncol = n_cages,
                          dimnames = list(cages, cages))

    for (i in seq_along(cages)) {
      for (j in seq_along(cages)) {
        cage_i = cages[i]
        cage_j = cages[j]
        founders_i = clust_data$founders[clust_data$cage == cage_i][[1]]
        founders_j = clust_data$founders[clust_data$cage == cage_j][[1]]

        # If both cages have no founders in this cluster, set NA (skip in average)
        # If one has founders and other doesn't, Jaccard = 0
        # If both have founders, calculate normally
        if (length(founders_i) == 0 && length(founders_j) == 0) {
          clust_matrix[cage_i, cage_j] = NA
        } else {
          clust_matrix[cage_i, cage_j] = jaccard(founders_i, founders_j)
        }
      }
    }

    cluster_matrices[[clust]] = clust_matrix
  }

  # Aggregate: average Jaccard across clusters
  # Includes 0s where one cage has data and other doesn't
  jaccard_matrix = matrix(0, nrow = n_cages, ncol = n_cages,
                          dimnames = list(cages, cages))

  for (i in cages) {
    for (j in cages) {
      values = sapply(cluster_matrices, function(m) m[i, j])
      jaccard_matrix[i, j] = mean(values, na.rm = TRUE)
    }
  }

  # Hierarchical clustering to reorder cages
  dist_matrix = as.dist(1 - jaccard_matrix)  # Convert similarity to distance
  hc = hclust(dist_matrix, method = "complete")
  cage_order = hc$labels[hc$order]

  # Convert to long format for ggplot
  jaccard_df = as.data.frame(as.table(jaccard_matrix)) %>%
    rename(cage1 = Var1, cage2 = Var2, similarity = Freq) %>%
    mutate(
      cage1 = factor(cage1, levels = cage_order),
      cage2 = factor(cage2, levels = cage_order)
    ) %>%
    # Keep only upper triangle (where cage2 index >= cage1 index)
    filter(as.numeric(cage2) >= as.numeric(cage1))

  # Create heatmap
  p = ggplot(jaccard_df, aes(x = cage1, y = cage2, fill = similarity)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(similarity, 2)), size = 2.5) +
    scale_fill_gradient(low = "white", high = "#3182bd", limits = c(0, 1)) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    labs(
      title = "",
      x = "Population",
      y = "Population",
      fill = "Mean Jaccard\nSimilarity"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    ) +
    coord_fixed()

  return(p)
}

# ============================================================================
# IMPORT SIGNIFICANCE DATA
# ============================================================================

# Read significance test results for treatment-time interaction at F60 (100kb windows)
sig = read.csv("data/treatment_time_repl_F60_100000.csv") %>%
  # Keep only main chromosomes
  dplyr::filter(CHR %in% chr_order) %>%
  # Set chromosome order as factor for consistent plotting
  dplyr::mutate(CHR = factor(CHR, levels = chr_order))

# ============================================================================
# IMPORT HAPLOTYPE FREQUENCY DATA
# ============================================================================

# Read sample labels and create file paths for haplotype frequency data
labels = read.csv("data/pca/sync_labs.csv") %>%
  mutate(pop_file = paste0("data/haplo_frqs/f", row_number(), ".tsv"),
         pop_id = paste0("f",row_number())) %>%
  # Filter to selected treatment
  filter(treatment == treatment_i)

# Read founder genotype names
founders = read.csv("data/haplo_frqs/foundergt.names", header = FALSE, col.names = "name")

# Add index for data frame matching
labels$idx = row_number(labels)

# Read all haplotype frequency files into list
dfs = lapply(labels$pop_file, read.csv, sep = "\t")

# ============================================================================
# CLUSTER SIGNIFICANT SNPS
# ============================================================================

# Extract SNPs passing Bonferroni threshold
snps = sig %>% filter(P < bf)

# Apply clustering to each chromosome with specified number of clusters
clust = rbind(
  cluster_snps("2L", 1),
  cluster_snps("2R", 3),
  cluster_snps("3L", 1),
  cluster_snps("3R", 3)
)

# Calculate genomic boundaries (start/end positions) for each cluster
clust_bounds = clust %>%
  group_by(cluster, CHR) %>%
  summarise(
    start = min(BP, na.rm = TRUE),
    end = max(BP, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================================
# SUBSET HAPLOTYPE DATA BY CLUSTER BOUNDARIES
# ============================================================================

# Apply subsetting to all haplotype frequency data frames
dfs = lapply(dfs, subset_by_clusters, bounds = clust_bounds)

# ============================================================================
# AVERAGE FREQUENCIES BY CLUSTER
# ============================================================================

# Apply averaging to all data frames
dfs_avg = lapply(dfs, average_freqs_by_cluster)

# ============================================================================
# CONVERT TO LONG FORMAT
# ============================================================================

# Combine all averaged data frames into a single data frame
df_avg = bind_rows(dfs_avg) %>%
  left_join(labels, by = join_by(population==pop_id), relationship="many-to-one")

long_df = make_long(df_avg)

# ============================================================================
# PREPARE MANHATTAN PLOT DATA
# ============================================================================

# Calculate cumulative chromosome sizes for continuous x-axis
chr_sizes = sig %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP)) %>%
  arrange(CHR) %>%
  # Each chromosome starts where the previous one ended
  mutate(chr_start = lag(cumsum(chr_len), default = 0))

# Add cumulative genomic position to significance data
sig = sig %>%
  left_join(chr_sizes, by = "CHR") %>%
  mutate(pos_cum = BP + chr_start)

# Calculate center position of each chromosome for axis labels
axis_df = chr_sizes %>%
  mutate(center = chr_start + chr_len / 2)

# Add cumulative positions to cluster bounds for plotting
clust_bounds = clust_bounds %>%
  left_join(chr_sizes, by = "CHR") %>%
  mutate(
    start_cum = start + chr_start,
    end_cum = end + chr_start,
    center_cum = (start_cum + end_cum) / 2
  )

# Assign point colors based on significance and chromosome
sig = sig %>%
  mutate(
    chr_index = as.numeric(CHR),
    # Blue if significant, alternating grey shades otherwise
    chr_color = if_else(
      P < bf,
      "#3182bd",
      if_else(chr_index %% 2 == 0, "grey60", "grey30")
    )
  )

# ============================================================================
# CREATE MANHATTAN PLOT
# ============================================================================

# Position cluster boxes above the highest data point
y_max = max(-log10(sig$P), na.rm = TRUE)
box_ymin = y_max + 0.5
box_ymax = y_max + 1.5

p1 = ggplot(sig, aes(x = pos_cum, y = -log10(P))) +
  # Cluster boundary boxes above the plot
  geom_rect(
    data = clust_bounds,
    aes(xmin = start_cum, xmax = end_cum, ymin = box_ymin, ymax = box_ymax, fill = cluster),
    inherit.aes = FALSE,
    alpha = 0.7
  ) +
  geom_text(
    data = clust_bounds,
    aes(x = center_cum, y = (box_ymin + box_ymax) / 2, label = cluster),
    inherit.aes = FALSE,
    size = 2.5,
    fontface = "bold"
  ) +
  # SNP points
  geom_point(aes(color = chr_color), size = 0.7) +
  scale_color_identity() +
  # Chromosome labels at center positions
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$CHR
  ) +
  # Expand y-axis to accommodate cluster boxes
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = "Chromosome", y = expression(-log[10](P))) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  # Bonferroni significance threshold line
  geom_hline(yintercept = -log10(bf), linetype = "dashed", color = "black")

p1

# ============================================================================
# FILTER INCREASED FOUNDERS AND ANNOTATE CLUSTER BOUNDS
# ============================================================================

# Subset to founders that increased in frequency by >0.1 between gen 20 and 60
long_df_increased = long_df %>%
  filter(gen %in% c(20, 60)) %>%
  group_by(cluster, cage, founder) %>%
  summarise(
    freq_change = as.numeric(frequencies[gen == 60]) - as.numeric(frequencies[gen == 20]),
    .groups = "drop"
  ) %>%
  filter(freq_change > 0.1) %>%
  # Join back to get all generations for the selected founders
  select(cluster, cage, founder) %>%
  inner_join(long_df, by = c("cluster", "cage", "founder"))

# Add column to clust_bounds listing founders that passed freq_change > 0.1 filter
founders_increased = long_df_increased %>%
  distinct(cluster, founder) %>%
  group_by(cluster) %>%
  summarise(increased_founders = list(sort(unique(founder))), .groups = "drop")

clust_bounds = clust_bounds %>%
  left_join(founders_increased, by = "cluster") %>%
  select(cluster, CHR, start, end, increased_founders) %>%
  mutate(cluster = paste0(CHR, "(", start, ";", end, ")"))

clust_bounds %>%
  mutate(increased_founders = sapply(increased_founders, paste, collapse = ";")) %>%
  write.csv("data/clust_bounds.csv", row.names = F, quote = F)

# ============================================================================
# PLOT CLUSTER FREQUENCY TRAJECTORIES
# ============================================================================

print_clust = "3R.1"
p2 = plot_cluster_freqs(print_clust)
p2
print = p2 +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position="none",
    plot.title = element_blank()
  )

#ggsave(paste0("/Users/leahdarwin/Documents/drand/presentations/postdoc_talk/",print_clust,".pdf"), print,width = 5, height=4)

# ============================================================================
# JACCARD SIMILARITY HEATMAP
# ============================================================================

p_jaccard = plot_jaccard_heatmap(long_df_increased)
p_jaccard
ggsave(paste0("output/jaccard_",treatment_i,".pdf"), p_jaccard,width = 5, height=4)
