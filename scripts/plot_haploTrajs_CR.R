# =============================================================================
# Haplotype Frequency Trajectories: Control vs Rotenone
# =============================================================================
# Plots founder haplotype frequency trajectories across generations for each
# significant SNP cluster region, separately for control and rotenone
# treatments, using the increased founders identified in clust_bounds_founders.csv.
#
# Author: Leah Darwin
#
# Output files:
#   output/haplotraj_2L.pdf : Frequency trajectories for the 2L cluster
#   output/haplotraj_2R.pdf : Frequency trajectories for 2R clusters (stacked)
#   output/haplotraj_3R.pdf : Frequency trajectories for the 3R cluster
# =============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(patchwork)

eps = 1e-6

clust_bounds = read.csv("data/clust_bounds_founders.csv")

founders = read.csv("data/haplo_frqs/foundergt.names", header = FALSE, col.names = "name")

labels = read.csv("data/pca/sync_labs.csv") %>%
  mutate(pop_file   = paste0("data/haplo_frqs/f", row_number(), ".tsv"),
         population = paste0("f", row_number())) 

dfs = lapply(labels$pop_file, read.csv, sep = "\t")

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


windows_list =  paste0(clust_bounds$CHR,":",clust_bounds$start,"-",clust_bounds$end)                  

founders_list  = stringr::str_split(clust_bounds$increased_founders, ";")

# windows_list = c("2L:18000000-23000000",
#                  "2R:11000000-12500000",
#                  "2R:14000000-16000000",
#                  "3R:19400000-23300000")
# founders_list = list(c("Ore.Ore","DGRP.375","Bei38"),
#                      c("Ore.Ore"),
#                      c("Bei42","w1118"),
#                      c("DGRP.375","Bei54")
# )

long_dfs = lapply(windows_list, get_long_df_window)
rm(dfs)
gc()

plot_traj = function(df, window, founders){
  
  plots = list()
  
  for(founder_i in 1:length(founders)){
    
    temp = df %>%
      filter(founder == founders[[founder_i]]) %>%
      summarise(freq = mean(freq), .by = c(population,founder))%>%
      left_join(labels, by = join_by(population)) %>%
      select(-pop_file)
    
    plots[[founder_i]] = ggplot(temp, aes(x = gen, y = freq, color = treatment, group=interaction(treatment,cage))) +
      geom_point() +
      geom_line() +
      facet_wrap(~founder) +
      labs(x = "Generation", y = "Estimated Frequency") +
      theme_bw()+
      scale_color_manual(values = c("C"="black", "R"="red")) +
      theme(strip.background = element_rect(fill = "white", color = "black"),
            strip.text = element_text(color = "black"),
            )
  }
  
  wrap_plots(plots) + plot_layout(guides = "collect", axes="collect") + plot_annotation(title=windows_list[[idx]]) & theme(legend.position = "bottom")
  
}

patches = list()

for(idx in 1:length(clust_bounds)){
  
  patches[[idx]] = plot_traj(long_dfs[[idx]], windows_list[[idx]], founders_list[[idx]])
  
}

ggsave("output/haplotraj_2L.pdf",patches[[1]], width=6.5, height=3)

final_2R = wrap_elements(patches[[4]]) / wrap_elements(patches[[2]]) / wrap_elements(patches[[3]])
ggsave("output/haplotraj_2R.pdf", final_2R, width=6.5, height=9)

ggsave("output/haplotraj_3R.pdf",patches[[5]],width = 6.5, height=3)
