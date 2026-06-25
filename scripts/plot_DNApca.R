# =============================================================================
# PCA of Pooled Allele Frequency Data
# =============================================================================
# Performs PCA on genome-wide pooled sequencing allele frequencies across all
# samples. Plots PC1 vs PC3 by treatment with generation ellipses, PC2 vs
# generation by cage origin, and a scree plot.
#
# Author: Leah Darwin
#
# Output files:
#   output/pca_allelefreqs.pdf : Combined figure with PC1 vs PC3 by treatment
#                               (left), PC2 vs generation by cage origin and
#                               scree plot (right)
# =============================================================================

# ============================================================================
# LOAD PACKAGES
# ============================================================================

library(ggplot2)      # For creating plots
library(dplyr)        # For data manipulation
library(RColorBrewer) # For color palettes (not actively used)
library(ggrepel)      # For non-overlapping text labels
library(patchwork)    # For combining multiple plots
library(plotly)       # For interactive plots (not actively used)

# ============================================================================
# DATA IMPORT AND PREPROCESSING
# ============================================================================

# Read allele frequency data from pooled sequencing (MAF > 0.01)
# Filter to keep only rows with finite values across all sample columns (4:123)
afs = read.csv("data/pca/joined.sync.MAF01.frq", sep = '\t')  %>%
  filter(if_all(4:123, is.finite)) %>%
  na.omit()

# Read sample labels and create additional grouping variables
labels = read.csv("data/pca/sync_labs.csv") %>%
  mutate(orig4 = case_when(startsWith(cage, "1A") ~ "1A",
                           startsWith(cage, "1B") ~ "1B",
                           startsWith(cage, "2A") ~ "2A",
                           startsWith(cage, "2B") ~ "2B")) %>%
  mutate(orig2 = case_when(startsWith(cage, "1") ~ "1",
                           startsWith(cage, "2") ~ "2")) %>%
  mutate(gen = factor(gen))

# ============================================================================
# PRINCIPAL COMPONENT ANALYSIS
# ============================================================================

# Perform PCA on transposed allele frequency matrix
# - Each row represents a sample, each column represents a SNP
# - center = T: center variables to have mean zero
# - scale. = F: do not scale to unit variance
pc = prcomp(t(as.matrix(afs[,4:123])), center = T, scale. = F)

# Clean up memory
rm(afs)
gc()

# ============================================================================
# PCA RESULTS PREPARATION
# ============================================================================

# Extract PC scores and combine with sample labels
pca_df = as.data.frame(pc$x)
pca_df = cbind(pca_df, labels)

# Calculate variance explained by each PC
pca_var <- pc$sdev^2
pca_var_per <- round(pca_var / sum(pca_var) * 100, 1)

# Define custom colors and labels for plotting
custom_colors = c("C" = "black", "R" = "red")
custom_labels = c("C" = "Control", "R" = "Rotenone")
custom_colors_cage = c("1" = "#9ecae1", "2" = "#3182bd")

# Create interaction group for treatment x generation
pca_df$group = interaction(pca_df$treatment, pca_df$gen)

# ============================================================================
# CALCULATE CENTROIDS FOR ELLIPSES
# ============================================================================

# Calculate mean PC scores for each treatment-generation combination
# Exclude generations 22 and 40 from centroid calculation
centroids <- pca_df %>%
  filter(!gen %in% c(22, 40)) %>%
  group_by(treatment, gen, group) %>%
  summarize(PC1 = mean(PC1), PC2 = mean(PC2), PC3 = mean(PC3), .groups = 'drop')

# ============================================================================
# PREPARE SCREE PLOT DATA
# ============================================================================

# Calculate percentage variance explained by each PC
var_values <- pc$sdev^2
var_pct <- var_values / sum(var_values) * 100

# Create data frame for scree plot (first 10 PCs)
scree_data <- data.frame(
  PC = 1:10,
  Variance = var_pct[1:10]
)

# Ensure PCs stay in order 1-10 on the x-axis
scree_data$PC <- factor(scree_data$PC, levels = scree_data$PC)

# ============================================================================
# PLOT 1: PC1 vs PC3 with Treatment Groups
# ============================================================================

p1 = ggplot(pca_df, aes(x = PC1, y = PC3, color = treatment)) +
  geom_point(size=2) +
  theme_classic()+
  # Add ellipses for treatment-generation groups (excluding gens 22, 40)
  stat_ellipse(data = filter(pca_df, !gen %in% c(22, 40)),
               aes(group = group),
               linetype = 1) +
  # Add labels at centroid of each group
  geom_label_repel(data = centroids, aes(label = paste("F",gen,sep=""), color=treatment),
                   box.padding = 0.5, segment.color = 'grey50') +
  labs(
    x = paste0("PC1 (", pca_var_per[1], "%)"),
    y = paste0("PC3 (", pca_var_per[3], "%)")
  )  +
  scale_color_manual(values=custom_colors)+
  theme(legend.position = "none")

p1

# ============================================================================
# PLOT 2: PC2 vs Generation (by Origin Population)
# ============================================================================

p2 = ggplot(pca_df, aes(x = paste("F",gen,sep=""), y = PC2, color = orig2)) +
  geom_jitter(width = 0.2,  size=2) +
  theme_classic()+
  labs(
    x = paste0("Generation"),
    y = paste0("PC2 (", pca_var_per[2], "%)")
  )  +
  scale_color_manual(values=custom_colors_cage) +
  theme(legend.position = "none")
p2

# ============================================================================
# PLOT 3: Scree Plot (Variance Explained)
# ============================================================================

p3 = ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_col(fill = "grey15") +
  labs( y = "Variance \nExplained (%)",
       x = "Principal Component") +
  theme_classic()

p3

# ============================================================================
# COMBINE PLOTS INTO FINAL FIGURE
# ============================================================================

# Layout: p1 on left, (p2 stacked on p3) on right
# p2 takes up 2/3 of right column height, p3 takes up 1/3
final = (p1 | (p2 / p3)+plot_layout(heights = c(2,1))) +
  plot_layout(widths = c(2, 1))+
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')')
ggsave("output/pca_allelefreqs.pdf", final, width = 8, height = 5)


