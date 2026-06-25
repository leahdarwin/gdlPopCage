# =============================================================================
# Final Figure: Haplotype Frequency Bar Plots
# =============================================================================
# Plots stacked bar charts of mean founder haplotype frequencies across cages
# for specific genomic windows of interest, separately for control and rotenone
# treatments. Generation to plot is set by gen_i in the PARAMETERS section.
#
# Author: Leah Darwin
#
# Output files:
#   output/barfreqs_F[gen_i].pdf : Stacked bar plots of founder frequencies at
#                                  the specified generation for control (top)
#                                  and rotenone (bottom) treatments
# =============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

eps = 1e-6

# ============================================================================
# PARAMETERS
# ============================================================================

# Generation to plot (must exist in sync_labs.csv, e.g. 20, 40, 60)
gen_i = 20

# Set to TRUE to show all founders, FALSE to use founders_to_show
show_all_founders = FALSE

founder_colors = c(
  "w1118"    = "#ec7014",
  "Bei42"    = "#4292c6",
  "Bei23"    = "#084594",
  "Bei54"    = "#08519c",
  "Bei59"    = "#08619e",
  "Bei52"    = "#0570b0",
  "Bei12"    = "#2171b5",
  "Bei38"    = "#6baed6", 
  "ZH26"     = "#9ebcda",##
  "ZH23"     = "#8c96c6", ##
  "ZW144"    = "#8c6bb1",##
  "ZW142"    = "#88419d",##
  "DGRP.375" = "#fec44f", 
  "Ore.Ore"  = "#fe9929",
  "Other"    = "grey90"
)


# ============================================================================
# DATA IMPORT
# ============================================================================

founders = read.csv("data/haplo_frqs/foundergt.names", header = FALSE, col.names = "name")

labels = read.csv("data/pca/sync_labs.csv") %>%
  mutate(pop_file   = paste0("data/haplo_frqs/f", row_number(), ".tsv"),
         population = paste0("f", row_number())) %>%
  filter(gen == gen_i)

dfs = lapply(labels$pop_file, read.csv, sep = "\t")

clust_bounds = read.csv("data/clust_bounds_founders.csv") %>%
  filter()
  

# ============================================================================
# FILTER TO WINDOW AND CONVERT TO LONG FORMAT
# ============================================================================

make_long_window = function(window_id, df) {
  
  window_chr   = sub(":.*", "", window_id)
  window_start = as.integer(sub(".*:(\\d+)-.*", "\\1", window_id))
  window_end   = as.integer(sub(".*-(\\d+)$",   "\\1", window_id))
  
  df %>%
    filter(chr == window_chr, pos >= window_start, pos <= window_end) %>%
    separate_rows(frequencies, sep = ";") %>%
    mutate(freq = as.numeric(frequencies) + eps) %>%
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

# ============================================================================
# PLOT
# ============================================================================

plot_bars = function(long_df, founders_to_show, window_id, gen_i, treatment_i) {
  active_founders = if (show_all_founders) unique(long_df$founder) else founders_to_show
  
  subset = long_df %>%
    left_join(labels, by = "population") %>%
    filter(gen == gen_i, treatment == treatment_i,
           founder %in% active_founders) %>%
    group_by(cage, founder) %>%
    summarise(freq = mean(as.numeric(frequencies), na.rm = TRUE), .groups = "drop")
  
  remainder = subset %>%
    group_by(cage) %>%
    summarise(freq = pmax(0, 1 - sum(freq)), founder = "Other", .groups = "drop")
  
  bind_rows(subset, remainder) %>%
    mutate(founder = factor(founder, levels = c("Other", active_founders))) %>%
    mutate(cage_hist = case_when(grepl("^1", cage) ~ "1",
                                 grepl("^2", cage) ~ "2")) %>%
    ggplot(aes(x = cage, y = freq, fill = founder)) +
    geom_bar(position = "stack", stat = "identity") +
    facet_wrap(~paste("Hist.", cage_hist), strip.position="bottom", scales="free_x") +
    scale_fill_manual(values = founder_colors) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(title = {
           chr   = sub(":.*", "", window_id)
           start = round(as.integer(sub(".*:(\\d+)-.*", "\\1", window_id)) / 1e6,1)
           end   = round(as.integer(sub(".*-(\\d+)$",   "\\1", window_id)) / 1e6,1)
           paste0(chr, ":", start, "-", end, "Mbp")
         }, x = "", y = "Mean estimated\nfounder frequency") +
    theme_classic() +
    theme(axis.text.x  = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 10),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          axis.line    = element_blank())
}


# windows_list =  paste0(clust_bounds$CHR,":",clust_bounds$start,"-",clust_bounds$end)                  
# 
# founders_list  = stringr::str_split(clust_bounds$increased_founders, ";")

windows_list = c("2L:18000000-23000000",
                 "2R:11000000-12500000",
                 "2R:14000000-16000000",
                 "3R:19400000-23300000")

founders_list = list(c("Ore.Ore","DGRP.375","Bei38"),
                     c("Ore.Ore"),
                     c("Bei42","w1118"),
                     c("DGRP.375","Bei54")
                    )

long_dfs = lapply(windows_list, get_long_df_window)

plots_C = mapply(plot_bars,
       long_df = long_dfs,
       founders_to_show = founders_list,
       window_id = windows_list,
       MoreArgs = list(gen_i = gen_i, treatment_i = "C")
)

plots_R = mapply(plot_bars,
                 long_df = long_dfs,
                 founders_to_show = founders_list,
                 window_id = windows_list,
                 MoreArgs = list(gen_i = gen_i, treatment_i = "R")
)

patch_C = wrap_plots(plots_C, nrow=1) +plot_layout(axes = "collect_y") & theme(legend.position = "none")
patch_R = wrap_plots(plots_R, nrow=1) +plot_layout(axes = "collect_y") & theme(legend.position = "none")

final = wrap_elements(patch_C)/wrap_elements(patch_R)+plot_annotation(tag_levels = "a")
final
ggsave(paste0("output/barfreqs_F", gen_i, ".pdf"), final, width=10, height=6)
