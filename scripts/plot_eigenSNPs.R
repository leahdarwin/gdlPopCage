# =============================================================================
# Eigen Analysis of Allele Frequency Change Vectors
# =============================================================================
# Uses the afvaper package to compute eigendecomposition of allele frequency
# change vectors (rotenone vs control) per genomic window along chromosomes 2L,
# 2R, and 3R, with permutation-based significance thresholds.
#
# Author: Leah Darwin
#
# Output files:
#   output/eigenplots_final.pdf : EV1 and EV2 proportion of variance plots for
#                                 chromosomes 2L, 2R, and 3R
# =============================================================================

library(afvaper)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(patchwork)

chrs=c("2L","2R","3R")

clust_bounds = read.csv("data/clust_bounds.csv")
clust_col = "#3182bd"

labels = read.csv("data/popid.csv") %>%
  filter(gen %in% c(20,60), treatment == "R")

labels_all = read.csv("data/popid.csv") %>%
  filter(gen %in% c(20,60))

freq_mat = read.csv("data/joined.sync.MAF01.frq", sep="\t") %>%
  filter(chrom%in%chrs) %>%
  dplyr::select(chrom, pos, all_of(labels_all$pop_id))

# Set window size
window_snps = 500

# Build per-cage, per-treatment lookup of pop_ids at each generation
cage_df = labels_all %>%
  arrange(cage, treatment, gen) %>%
  group_by(cage, treatment) %>%
  summarise(pop_20 = pop_id[gen == min(gen)],
            pop_60 = pop_id[gen == max(gen)],
            .groups = "drop")

# Add difference columns to freq_mat: (freq_60 - freq_20) per cage x treatment
# Prefix with "diff_" so names are valid R identifiers without any transformation
for (i in seq_len(nrow(cage_df))) {
  col_name = paste0("diff_", cage_df$cage[i], "_", cage_df$treatment[i])
  freq_mat[[col_name]] = freq_mat[[cage_df$pop_60[i]]] - freq_mat[[cage_df$pop_20[i]]]
}

# Build vector_list: each element is c(C_diff_col, R_diff_col) for the same cage
# calc_AF_vectors computes vec[2] - vec[1] = R_diff - C_diff
#   = (R_60 - R_20) - (C_60 - C_20)
cages = unique(cage_df$cage)
vector_list = setNames(
  lapply(cages, function(cage_i) {
    c(paste0("diff_", cage_i, "_C"), paste0("diff_", cage_i, "_R"))
  }),
  cages
)

make_plot = function(chr_i){
  
  freq_chr_i = freq_mat %>%
    filter(chrom==chr_i)
  
  af_obj = calc_AF_vectors(vcf=freq_chr_i, 
                           window_size = window_snps,
                           vectors = vector_list,
                           n_cores = 4,
                           data_type = "freq")
  
  # How many permutations to run
  null_perm_N = 10000
  
  # Calculate Allele Frequency Change Vector Matrices
  null_input = calc_AF_vectors(vcf = freq_chr_i,
                                window_size = window_snps,
                                vectors = vector_list,
                                n_cores = 4,
                                null_perms = null_perm_N,
                                data_type = "freq")
  
  # Perform eigen analysis
  eigen_res = lapply(af_obj,eigen_analyse_vectors)
  null_eigen_res = lapply(null_input, eigen_analyse_vectors)
  
  # Proportion of variance explained by EV1 and EV2 per window
  eigen_prop = do.call(rbind, lapply(names(eigen_res), function(w) {
    vals  = eigen_res[[w]]$eigenvals
    total = sum(vals)
    data.frame(
      window = w,
      chr    = sub(":.*", "", w),
      pos    = as.integer(sub(".*:(\\d+)-.*", "\\1", w)),
      EV1    = vals[1] / total,
      EV2    = vals[2] / total
    )
  }))

  
  null_props = do.call(rbind, lapply(null_eigen_res, function(x) {
    vals  = x$eigenvals
    total = sum(vals)
    data.frame(EV1 = vals[1] / total, EV2 = vals[2] / total)
  }))
  
  prop_cutoffs = apply(null_props, 2, quantile, probs = c(0.95, 0.99, 0.999))
  
  chr_bounds = clust_bounds %>% filter(CHR == chr_i)

  prop_p1 = ggplot(eigen_prop, aes(x = pos, y = EV1)) +
    geom_rect(data = chr_bounds,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = clust_col, alpha = 0.2) +
    geom_step() +
    geom_hline(yintercept = prop_cutoffs["99.9%",  "EV1"], linetype = "solid", colour = "red") +
    scale_x_continuous(labels = function(x) x / 1e6) +
    labs(x = paste0(chr_i," Position (Mbp)"), y = expression(paste("EV1/",Sum,"EV"))) +
    theme_classic()

  prop_p2 = ggplot(eigen_prop, aes(x = pos, y = EV2)) +
    geom_rect(data = chr_bounds,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = clust_col, alpha = 0.2) +
    geom_step() +
    geom_hline(yintercept = prop_cutoffs["99%",  "EV2"], linetype = "solid", colour = "red") +
    scale_x_continuous(labels = function(x) x / 1e6) +
    labs(x = paste0(chr_i," Position (Mbp)"), y = expression(paste("EV2/",Sum,"EV"))) +
    theme_classic()
  
  props = prop_p1/prop_p2
  return(props)
  
}

plots = lapply(chrs, make_plot)

final = wrap_elements(plots[[1]]) + wrap_elements(plots[[2]]) + wrap_elements(plots[[3]])

ggsave("output/eigenplots_final.pdf", final, width=8, height=3)


