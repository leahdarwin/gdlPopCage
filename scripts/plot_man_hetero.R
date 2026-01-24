library(dplyr)
library(ggplot2)
library(patchwork)

chr_order = c("2L","2R","3L","3R","4","X","Y")

bf = 0.05/3433901

sig = read.csv("../data/treatment_time_repl_F60_100000.csv") %>%
  dplyr::filter(CHR %in% c("2L","2R","3L","3R","4","X","Y"))%>%
  dplyr::mutate(CHR = factor(CHR, levels = chr_order))

labels = read.csv("../data/pca/sync_labs.csv") %>%
  mutate(pop_file = paste("../data/hetero/f", row_number(), "_H.tsv",sep = "")) %>%
  filter(gen == 60)

labels$idx = row_number(labels)

dfs = lapply(labels$pop_file,read.csv, sep="\t")

chr_sizes = sig %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP)) %>%
  arrange(CHR) %>%
  mutate(chr_start = lag(cumsum(chr_len), default = 0))

sig = sig %>%
  left_join(chr_sizes, by = "CHR") %>%
  mutate(pos_cum = BP + chr_start) 

axis_df = chr_sizes %>%
  mutate(center = chr_start + chr_len / 2)

sig = sig %>%
  mutate(chr_index = as.numeric(CHR),
         chr_color = if_else(P<bf, "#3182bd",
           if_else(chr_index %% 2 == 0, "grey60", "grey30")))

p1 = ggplot(sig, aes(x = pos_cum, y = -log10(P))) +
  geom_point(aes(color =chr_color), size = 0.7) +
  scale_color_identity() +
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$CHR
  ) +
  labs(x = "Chromosome", y = expression(-log[10](P))) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  geom_hline(yintercept = -log10(bf), linetype = "dashed", color="black")

plot_H = function(df_C, df_R, chr, title){
  
  df_C = df_C %>%
    filter(chrom == chr)
  df_R = df_R %>%
    filter(chrom == chr)
  
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
      limits = c(0, 0.3),   # optional but recommended
      expand = c(0, 0)
    ) #+ggtitle(title)
  
  
  if(chr=="3R"){
    p = p+geom_vline(xintercept = 21111506/1000000, color = "#3182bd", alpha=0.6)
  }else{
    p = p+geom_vline(xintercept = 14854642/1000000, color = "#3182bd", alpha=0.6)
  }
  
}

style_plots = function(plots){
  
  col1 <- 1:6
  col2 <- 7:12
  new_order <- as.vector(rbind(col1, col2))  # interleave
  plots_list <- plots[new_order]
  
  p = wrap_plots(plots_list ,ncol=2)&
    theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank()
    )
  
  left_idx <- seq(1, 12, by = 2)
  
  for (i in left_idx) {
    p[[i]] <- p[[i]] +
      theme(
        axis.text.y  = element_text(size=6),
        axis.ticks.y = element_line()
      )
  }
  
  bottom_idx <- c(11,12)
  
  for (i in bottom_idx) {
    p[[i]] <- p[[i]] +
      theme(
        axis.text.x  = element_text(size=8),
        axis.ticks.x = element_line()
      )
  }
  
  return(p)
  
}

cages = unique(labels$cage)


p_2R = by(labels, labels$cage, function(x) {
  plot_H(dfs[[x$idx[x$treatment == "C"]]],
         dfs[[x$idx[x$treatment == "R"]]], "2R", x$cage)
})


p2 = style_plots(p_2R)

p_3R = by(labels, labels$cage, function(x) {
  plot_H(dfs[[x$idx[x$treatment == "C"]]],
         dfs[[x$idx[x$treatment == "R"]]], "3R")
})

p3 = style_plots(p_3R)

patch = wrap_elements(p2)|wrap_elements(p3)
patcfinal <- p1  / (patch) &
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')', 
                  theme = theme(plot.margin = margin(1, 1, 1, 1)))

final

ggsave("../output/man_hetero.pdf",plot=final, width=6, height = 5)



