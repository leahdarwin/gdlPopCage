# =============================================================================
# Contingency Test: Initial Haplotype Composition vs. F60 Frequency
# =============================================================================
# Tests whether haplotype composition at an early generation predicts founder
# frequency at F60, using ILR-transformed PCA of haplotype frequencies and
# linear models with treatment interaction terms.
#
# Author: Leah Darwin
#
# Output files:
#   output/haploPCA_loadings_F[gen].pdf : PC1 CLR loadings per genomic window
#   output/haploPCA_windows_F[gen].pdf  : PCA scores (PC1 vs PC2) per window,
#                                         colored by treatment
#   output/haploPCA1_pred_F[gen].pdf    : PC1 score vs F60 founder frequency
#   output/haploPCA2_pred_F[gen].pdf    : PC2 score vs F60 founder frequency
#   output/haplofreq_pred_F[gen].pdf    : Direct founder frequency at F[gen]
#                                         vs F60 with bootstrap error bars
# =============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(car)
library(ggrepel)
library(compositions)
library(knitr)
library(tibble)

eps = 1e-6

# ============================================================================
# PARAMETERS
# ============================================================================

start_gen = 20

#=============================================================
# DATA IMPORT
# ============================================================================

founders = read.csv("data/haplo_frqs/foundergt.names", header = FALSE, col.names = "name")

labels = read.csv("data/pca/sync_labs.csv") %>%
  mutate(pop_file   = paste0("data/haplo_frqs/f", row_number(), ".tsv"),
         population = paste0("f", row_number())) %>%
  filter(gen == start_gen | gen == 60)

dfs = lapply(labels$pop_file, read.csv, sep = "\t")


# ============================================================================
# FILTER TO WINDOW AND CONVERT TO LONG FORMAT
# ============================================================================

make_long_window = function(window_id, df) {
  
  window_chr   = sub(":.*", "", window_id)
  window_start = as.integer(sub(".*:(\\d+)-.*", "\\1", window_id))
  window_end   = as.integer(sub(".*-(\\d+)$",   "\\1", window_id))
  
  has_sd   <- "boot_sd"   %in% names(df)
  has_bias <- "boot_bias" %in% names(df)
  sep_cols <- c("frequencies", if (has_sd) "boot_sd", if (has_bias) "boot_bias")

  df %>%
    filter(chr == window_chr, pos >= window_start, pos <= window_end) %>%
    separate_rows(all_of(sep_cols), sep = ";") %>%
    mutate(freq          = as.numeric(frequencies) + eps,
           boot_sd_val   = if (has_sd)   as.numeric(boot_sd)   else NA_real_,
           boot_bias_val = if (has_bias) as.numeric(boot_bias) else NA_real_) %>%
    group_by(chr, pos) %>%
    mutate(founder = founders$name[seq_len(n())]) %>%
    ungroup()
}

get_long_df_window = function(window_id){
  
  
  # ============================================================================
  # PARSE WINDOW
  # ============================================================================
  
  window_chr   = sub(":.*", "", window_id)
  window_start = as.integer(sub(".*:(\\d+)-.*", "\\1", window_id))
  window_end   = as.integer(sub(".*-(\\d+)$",   "\\1", window_id))
  
  dfs_long = lapply(dfs, function(df) make_long_window(window_id, df))
  
  # Remove any populations with no data in this window
  has_data = sapply(dfs_long, nrow) > 0
  dfs_long  = dfs_long[has_data]
  labels_w  = labels[has_data, ]
  
  long_df = bind_rows(mapply(function(df, pop) mutate(df, population = pop),
                             dfs_long, labels_w$population, SIMPLIFY = FALSE))
  
}

get_avgs = function(df){
  df %>%
    summarize(freq      = mean(freq),
              pooled_sd = sqrt(mean(boot_sd_val^2, na.rm = TRUE)),
              .by = c(population, founder)) %>%
    left_join(labels, by = join_by(population), relationship="many-to-one")
}

windows_list = c("2L:18000000-23000000",
                 "2R:11000000-12500000",
                 "2R:14000000-16000000",
                 "3R:19400000-23300000")
long_dfs = lapply(windows_list, get_long_df_window)

founders_list = list(c("Ore.Ore","DGRP.375","Bei38"),
                     c("Ore.Ore"),
                     c("Bei42","w1118"),
                     c("DGRP.375","Bei54")
)

avg_dfs = lapply(long_dfs, get_avgs)

format_window_mb = function(window_id) {
  chr   = sub(":.*", "", window_id)
  start = as.numeric(sub(".*:(\\d+)-.*", "\\1", window_id)) / 1e6
  end   = as.numeric(sub(".*-(\\d+)$",   "\\1", window_id)) / 1e6
  paste0(chr, ":", start, "-", end, " Mb")
}

# Back-transform ILR PCA loadings through the basis matrix V to get
# per-founder contributions in CLR space (one value per founder per PC).
plot_pca_loadings = function(pca_res, V, founder_names, highlight_founders = character(0), title = "") {
  founder_load = V %*% pca_res$rotation
  rownames(founder_load) = founder_names

  load_df = as.data.frame(founder_load) %>%
    rownames_to_column("founder") %>%
    pivot_longer(-founder, names_to = "PC", values_to = "loading") %>%
    filter(PC == "PC1") %>%
    mutate(highlighted = founder %in% highlight_founders,
           founder     = reorder(founder, -abs(loading)))

  ggplot(load_df, aes(x = founder, y = loading, fill = highlighted)) +
    geom_col() +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "grey"), guide = "none") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Founder", y = paste0("PC1 CLR loading (F", start_gen, ")"), title = title)
}

# Plot frequency of founder_x at start_gen against frequency of founder_y at F60,
# coloured by treatment. window_id must match an entry in windows_list.
plot_cross_founder = function(window_id, founder_x, founder_y) {
  idx    = match(window_id, windows_list)
  avg_df = avg_dfs[[idx]]

  x_df = avg_df %>%
    filter(gen == start_gen, founder == founder_x) %>%
    select(cage, treatment, freq_x = freq)

  y_df = avg_df %>%
    filter(gen == 60, founder == founder_y) %>%
    select(cage, treatment, freq_y = freq)

  plot_df = inner_join(x_df, y_df, by = c("cage", "treatment"))

  ggplot(plot_df, aes(x = freq_x, y = freq_y, color = treatment)) +
    geom_point() +
    theme_classic() +
    scale_color_manual(values = c("R" = "red", "C" = "black")) +
    labs(x     = paste0(founder_x, " frequency at F", start_gen),
         y     = paste0(founder_y, " frequency at F60"),
         title = paste0(format_window_mb(window_id), " | ", founder_x, " vs ", founder_y))
}

# ============================================================================
# LOOP OVER WINDOWS AND FOUNDERS
# ============================================================================

scree_plots    = list()
pca_plots      = list()
loading_plots  = list()
pc1_plots      = list()
pc2_plots      = list()
direct_plots   = list()
anova_rows     = list()
hyp_rows       = list()

for (idx in seq_along(windows_list)) {

  target_window = format_window_mb(windows_list[idx])
  target_avgs   = avg_dfs[[idx]]

  # --- ILR / PCA (once per window) ------------------------------------------

  f_start_wide = target_avgs %>%
    filter(gen == start_gen) %>%
    pivot_wider(id_cols = c(cage, treatment),
                names_from  = founder,
                values_from = freq)

  row_meta = f_start_wide %>% select(cage, treatment)
  comp_mat = f_start_wide %>% select(-cage, -treatment) %>% as.matrix()
  print(paste("Row sums of comp_mat for", target_window))                                              
  print(rowSums(comp_mat))  
  ilr_mat  = ilr(acomp(comp_mat))
  pca_res  = prcomp(ilr_mat, center = TRUE, scale. = FALSE)
  ilr_V    = ilrBase(acomp(comp_mat))

  scree_df = data.frame(PC  = seq_len(length(pca_res$sdev)),
                        var = pca_res$sdev^2 / sum(pca_res$sdev^2))

  pc_df = row_meta %>%
    mutate(PC1 = pca_res$x[, 1],
           PC2 = pca_res$x[, 2])

  loading_plots[[target_window]] = plot_pca_loadings(pca_res, ilr_V,
                                                     colnames(comp_mat),
                                                     highlight_founders = founders_list[[idx]],
                                                     target_window)

  scree_plots[[target_window]] = ggplot(scree_df, aes(x = PC, y = var)) +
    geom_col() + geom_line() + geom_point() +
    scale_x_continuous(breaks = scree_df$PC) +
    labs(x = "PC", y = "Proportion of variance", title = target_window) +
    theme_classic()

  pca_plots[[target_window]] = ggplot(pc_df, aes(x = PC1, y = PC2, color = treatment, label = cage)) +
    geom_point() +
    geom_text_repel(size = 3, show.legend = FALSE) +
    theme_classic() +
    scale_color_manual(values = c("R" = "red", "C" = "black")) +
    labs(x     = paste0("PC1 (", round(scree_df$var[1] * 100, 1), "%)"),
         y     = paste0("PC2 (", round(scree_df$var[2] * 100, 1), "%)"),
         title = paste0(target_window, " | F", start_gen))

  # --- Per-founder tests -----------------------------------------------------

  for (target_founder in founders_list[[idx]]) {

    combo_key = paste0(target_window, " | ", target_founder)
    print(combo_key)

    f60_freq = target_avgs %>%
      filter(gen == 60, founder == target_founder) %>%
      select(cage, treatment, freq_60 = freq, sd_60 = pooled_sd)

    model_df = left_join(pc_df, f60_freq, by = c("cage", "treatment"))

    pc2_plots[[combo_key]] = ggplot(model_df, aes(x = PC2, y = freq_60, color = treatment)) +
      geom_errorbar(aes(ymin = freq_60 - 1.96*sd_60, ymax = freq_60 + 1.96*sd_60), width = 0.05) +
      geom_point() +
      theme_classic() +
      scale_color_manual(values = c("R" = "red", "C" = "black")) +
      labs(x     = paste0("PC2 of Haplotype Frequencies at F", start_gen),
           y     = "Estimated Average Frequency at F60",
           title = combo_key)

    lm_pca     = lm(freq_60 ~ 0 + treatment + treatment:PC1 + treatment:PC2,
                    data = model_df)
    coefs_pca  = as.data.frame(summary(lm_pca)$coefficients) %>%
      rownames_to_column("term") %>%
      mutate(window = target_window, founder = target_founder, model = "PCA")
    anova_rows[[paste0(combo_key, "_pca")]] = coefs_pca

    lh_pca = if (df.residual(lm_pca) > 0)
      linearHypothesis(lm_pca, "treatmentC:PC1 = treatmentR:PC1")
    else NULL
    hyp_rows[[paste0(combo_key, "_pca")]] = data.frame(
      window  = target_window, founder = target_founder, model = "PCA",
      F       = if (!is.null(lh_pca)) lh_pca[2, "F"]       else NA_real_,
      df_res  = if (!is.null(lh_pca)) lh_pca[2, "Res.Df"]  else NA_real_,
      p       = if (!is.null(lh_pca)) lh_pca[2, "Pr(>F)"]  else NA_real_
    )

    fe      = coef(lm_pca)
    pred_df = bind_rows(lapply(c("C", "R"), function(trt) {
      pc1_term  = paste0("treatment", trt, ":PC1")
      p_val     = coefs_pca$`Pr(>|t|)`[coefs_pca$term == pc1_term]
      if (length(p_val) > 0 && p_val < 0.05) {
        trt_df    = filter(model_df, treatment == trt)
        pc1_range = range(trt_df$PC1)
        intercept = fe[[paste0("treatment", trt)]]
        slope     = fe[[pc1_term]]
        data.frame(treatment = trt,
                   PC1       = pc1_range,
                   freq_60   = intercept + slope * pc1_range)
      }
    }))

    pc1_plots[[combo_key]] = ggplot(model_df, aes(x = PC1, y = freq_60, color = treatment)) +
      geom_errorbar(aes(ymin = freq_60 - 1.96*sd_60, ymax = freq_60 + 1.96*sd_60), width = 0.05) +
      geom_point() +
      { if (nrow(pred_df) > 0) geom_line(data = pred_df, aes(x = PC1, y = freq_60, color = treatment)) } +
      theme_classic() +
      scale_color_manual(values = c("R" = "red", "C" = "black")) +
      labs(x     = paste0("PC1 of Haplotype Frequencies at F", start_gen),
           y     = "Estimated Average Frequency at F60",
           title = combo_key)

    pivot_df = target_avgs %>%
      filter(founder == target_founder) %>%
      pivot_wider(id_cols = c(cage, treatment, founder),
                  names_from = gen, values_from = c(freq, pooled_sd)) %>%
      rename(freq_start = paste0("freq_", start_gen),
             sd_start   = paste0("pooled_sd_", start_gen),
             sd_60      = "pooled_sd_60")

    lm_direct    = lm(freq_60 ~ 0 + treatment + freq_start:treatment, data = pivot_df)
    coefs_direct = as.data.frame(summary(lm_direct)$coefficients) %>%
      rownames_to_column("term") %>%
      mutate(window = target_window, founder = target_founder, model = "direct")
    anova_rows[[paste0(combo_key, "_direct")]] = coefs_direct

    lh_direct = if (df.residual(lm_direct) > 0)
      linearHypothesis(lm_direct, "treatmentC:freq_start = treatmentR:freq_start")
    else NULL
    hyp_rows[[paste0(combo_key, "_direct")]] = data.frame(
      window  = target_window, founder = target_founder, model = "direct",
      F       = if (!is.null(lh_direct)) lh_direct[2, "F"]       else NA_real_,
      df_res  = if (!is.null(lh_direct)) lh_direct[2, "Res.Df"]  else NA_real_,
      p       = if (!is.null(lh_direct)) lh_direct[2, "Pr(>F)"]  else NA_real_
    )

    fe_direct    = coef(lm_direct)
    pred_direct  = bind_rows(lapply(c("C", "R"), function(trt) {
      slope_term = paste0("treatment", trt, ":freq_start")
      p_val      = coefs_direct$`Pr(>|t|)`[coefs_direct$term == slope_term]
      if (length(p_val) > 0 && p_val < 0.05) {
        trt_df     = filter(pivot_df, treatment == trt)
        fs_range   = range(trt_df$freq_start)
        intercept  = fe_direct[[paste0("treatment", trt)]]
        slope      = fe_direct[[slope_term]]
        data.frame(treatment = trt,
                   freq_start = fs_range,
                   freq_60    = intercept + slope * fs_range)
      }
    }))

    direct_plots[[combo_key]] = ggplot(pivot_df, aes(x = freq_start, y = freq_60, color = treatment)) +
      geom_errorbarh(aes(xmin = freq_start - sd_start*1.96, xmax = freq_start + sd_start*1.96), height = 0.005) +
      geom_errorbar(aes(ymin = freq_60 - sd_60*1.96, ymax = freq_60 + sd_60*1.96), width = 0.005) +
      geom_point() +
      { if (nrow(pred_direct) > 0) geom_line(data = pred_direct, aes(x = freq_start, y = freq_60, color = treatment)) } +
      theme_classic() +
      scale_color_manual(values = c("R" = "red", "C" = "black")) +
      labs(x     = paste0("Estimated Average Frequency at F", start_gen),
           y     = "Estimated Average Frequency at F60",
           title = combo_key)
  }
}

# ============================================================================
# PATCHWORKS
# ============================================================================

wrap_plots(scree_plots)
p0 = wrap_plots(loading_plots, axis_titles="collect", guides="collect")
p1 = wrap_plots(pca_plots, axis_titles="collect",guides="collect") 
p2 = wrap_plots(pc1_plots, axis_titles="collect",guides="collect") 
p3 = wrap_plots(pc2_plots, axis_titles="collect",guides="collect") 
p4 = wrap_plots(direct_plots, axis_titles="collect",guides="collect")


ggsave(paste0("output/haploPCA_loadings_F",start_gen,".pdf"),p0, width = 9, height=7)
ggsave(paste0("output/haploPCA_windows_F",start_gen,".pdf"), p1,width = 9, height=7)
ggsave(paste0("output/haploPCA1_pred_F",start_gen,".pdf"), p2, width = 9, height=7)
ggsave(paste0("output/haploPCA2_pred_F",start_gen,".pdf"), p3, width = 9, height=7)
ggsave(paste0("output/haplofreq_pred_F",start_gen,".pdf"), p4, width = 9, height=7)

# ============================================================================
# ANOVA TABLES
# ============================================================================

format_coef_table = function(df) {
  p_col = "Pr(>|t|)"
  df %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)),
           !!p_col := ifelse(.data[[p_col]] < 0.001, "<0.001",
                             as.character(.data[[p_col]])))
}

all_anova = bind_rows(anova_rows) %>%
  select(window, founder, model, term, everything())

kable(all_anova %>% filter(model == "PCA") %>% select(-model) %>% format_coef_table(),
      format = "latex", booktabs = TRUE,
      caption = paste0("lm PCA coefficients: F", start_gen, " haplotype PC1/PC2 predicting F60 frequency"))

kable(all_anova %>% filter(model == "direct") %>% select(-model) %>% format_coef_table(),
      format = "latex", booktabs = TRUE,
      caption = paste0("lm direct coefficients: F", start_gen, " frequency predicting F60 frequency"))

all_hyp = bind_rows(hyp_rows) %>%
  mutate(across(c(F, df_res), ~ round(.x, 3)),
         p = ifelse(p < 0.001, "<0.001", as.character(round(p, 3))))

kable(all_hyp %>% filter(model == "PCA") %>% select(-model),
      format = "latex", booktabs = TRUE,
      caption = paste0("Treatment slope comparison (C vs R) — PCA model: F", start_gen, " PC1"))

kable(all_hyp %>% filter(model == "direct") %>% select(-model),
      format = "latex", booktabs = TRUE,
      caption = paste0("Treatment slope comparison (C vs R) — direct model: F", start_gen, " frequency"))
