library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)

eps = 1e-6
gen_i = 60
treatment_i = "C"

color_palette <- c(
  "Beijing"   = "#18458e",
  "Zimbabwe"  = "#4fc2ba",
  "OreR"  = "#ec2328",
  "DGRP375"= "#fbbd67",
  "w1118"    = "#c3c1e0"
)

labels = read.csv("../data/pca/sync_labs.csv") %>%
  mutate(pop_file = paste("../data/haplo_frqs/f", row_number(), ".tsv",sep = ""))%>%
  filter(gen == gen_i, treatment == treatment_i)

founders = read.csv("../data/haplo_frqs/foundergt.names",header = FALSE, col.names = "name")

labels$idx = row_number(labels)

dfs = lapply(labels$pop_file,read.csv, sep="\t")

make_long = function(df, chr_i){
  long = df %>%
    filter(chr == chr_i) %>%
    separate_rows(frequencies, sep = ";") %>%
    mutate(
      freq = as.numeric(frequencies) + eps,
    ) %>%
    group_by(chr, pos) %>%                           # each site
    mutate(founder = founders$name[seq_len(n())]) %>% # assign names in order
    ungroup() %>%
    mutate(founder_group = case_when(grepl("^B",founder) ~ "Beijing",
                                     grepl("^Z", founder) ~ "Zimbabwe",
                                     grepl("^O", founder) ~ "OreR",
                                     grepl("^D", founder) ~ "DGRP375",
                                     grepl("^w", founder) ~ "w1118"))
}

plot_frqs = function(df){
  
  df = df %>%
    mutate(
      founder_group = factor(
        founder_group,
        levels = c("OreR","DGRP375",  "w1118", "Beijing", "Zimbabwe")  # your desired order
      )
    )
  
  ggplot(df, aes(x = pos/1000000, y = freq, fill = founder_group)) +
    geom_col() +
    scale_y_continuous(
      breaks = c(0, 0.5, 1),
      expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(
      x = "Genomic position",
      y = "Founder frequency",
      fill = "Founder group"
    ) +
    theme_classic() +
    theme(
      #panel.border = element_rect(color = "black", fill = NA, linewidth=2)
    ) +
    scale_fill_manual(values=color_palette)
  
}


arrange_plots = function(chr){
  
  ##make data long
  df_chr = lapply(dfs, make_long, chr)
  
  ##make plots
  plots = lapply(df_chr, plot_frqs)
  
  ##reorder plots 
  col1 <- 1:6
  col2 <- 7:12
  new_order <- as.vector(rbind(col1, col2))  # interleave
  plots_list <- plots[new_order]
  
  p = wrap_plots(plots_list ,ncol=2)&
    #plot_layout(guides = 'collect')&
    theme(legend.position = "none") &
    theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(), 
      axis.text.x = element_blank(),
      plot.margin = margin(4, 4, 4, 4)
    ) 
  
  left_idx <- seq(1, 12, by = 2)
  
  ##set left side axis
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

p1 = arrange_plots("2R") 
p2 = arrange_plots("3R") 

legend = cowplot::get_legend(p1[[1]]+theme(legend.position = "right"))

final = (wrap_elements(p1) | wrap_elements(p2) )&
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')', 
                  theme = theme(plot.margin = margin(1, 1, 1, 1)))


ggsave(paste("../output/f",gen_i,"_",treatment_i,"_haplofrqs.png",sep=""), final, dpi = 300, width = 5, height = 3)
ggsave("../output/haplocolorsleg.pdf", ggdraw(legend))

