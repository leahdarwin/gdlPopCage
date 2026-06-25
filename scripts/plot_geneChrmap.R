# =============================================================================
# Gene Family Chromosome Map
# =============================================================================
# Maps the genomic positions of CYP450, Complex I, and UGT gene families onto
# a schematic chromosome diagram, overlaying the significant SNP cluster regions
# identified from the population-genetics analysis.
#
# Author: Leah Darwin
#
# Output files:
#   output/genemap.pdf : Chromosome map with gene family positions and cluster
#                        regions highlighted
# =============================================================================

library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)

##load gtf
gtf = read.csv(file="data/dmel-all-r6.57.gtf", sep="\t", header = FALSE) %>% 
  dplyr::filter(V3=="gene") %>%
  dplyr::select(V1, V4, V5, V9)%>%
  mutate(gene_id = str_extract(V9, "FBgn[0-9]+")) %>%
  mutate(gene_name = sub(".*gene_symbol ([^;]+);.*", "\\1", V9)) %>%
  dplyr::select(-V9) 

colnames(gtf) = c("chr","start","end", "gene_id", "gene_name")

gtf = gtf %>% filter(chr%in%c("X","2L","2R","3L","3R","4"))

cyps = read.csv("data/cyp450_genelist.txt",header=F,col.names = c("gene"))
cIs  = read.csv("data/complexI_genelist.txt",header=F,col.names = c("gene"))
ugts = read.csv("data/ugt_genelist.txt",,header=F,col.names = c("gene"))

clust_bounds = read.csv("data/clust_bounds.csv")

# ============================================================================
# CHROMOSOME MAP
# ============================================================================

chr_order = c("2L", "2R", "3L", "3R", "4", "X")
chr_x     = setNames(seq(1, by = 3, length.out = length(chr_order)), chr_order)
chr_w     = 0.3

chr_bounds = gtf %>%
  group_by(chr) %>%
  summarize(chr_start = 0, chr_end = max(end), .groups = "drop") %>%
  mutate(x = chr_x[chr])

gene_df = bind_rows(
  gtf %>% filter(gene_id %in% cyps$gene)  %>% mutate(gene_set = "CYP450"),
  gtf %>% filter(gene_id %in% cIs$gene)   %>% mutate(gene_set = "Complex I"),
  gtf %>% filter(gene_id %in% ugts$gene)  %>% mutate(gene_set = "UGT")
) %>%
  mutate(mid = (start + end) / 2,
         x   = chr_x[chr])

clust_df = left_join(clust_bounds, dplyr::select(chr_bounds, chr, x), by = join_by("CHR"=="chr"))

p1 = ggplot() +
  geom_rect(data = chr_bounds,
            aes(xmin = x - chr_w, xmax = x + chr_w,
                ymin = chr_start,  ymax = chr_end),
            fill = "grey85", colour = "grey50", linewidth = 0.3) +
  geom_rect(data = clust_df,
            aes(xmin = x - chr_w, xmax = x + chr_w,
                ymin = start, ymax = end),
            fill = "#3288bd", colour = NA, alpha=0.5) +
  geom_segment(data = gene_df,
               aes(x = x - chr_w, xend = x + chr_w,
                   y = mid, yend = mid,
                   colour = gene_set),
               linewidth = 0.6) +
  geom_text_repel(data = filter(gene_df, gene_set == "CYP450"),
                  aes(x = x + chr_w, y = mid, label = gene_name, colour = gene_set),
                  size         = 2.5,
                  nudge_x      = 0.5,
                  direction    = "y",
                  hjust        = 0,
                  segment.size = 0.3,
                  force        = 2,
                  show.legend  = FALSE) +
  geom_text_repel(data = filter(gene_df, gene_set %in% c("Complex I", "UGT")),
                  aes(x = x - chr_w, y = mid, label = gene_name, colour = gene_set),
                  size         = 2.5,
                  nudge_x      = -0.2,
                  direction    = "y",
                  hjust        = 1,
                  segment.size = 0.3,
                  force        = 2,
                  show.legend  = FALSE) +
  scale_x_continuous(breaks = chr_x, labels = names(chr_x),
                     expand = expansion(add = 1)) +
  scale_y_continuous(labels = function(x) paste0(x / 1e6, " Mb")) +
  scale_colour_manual(values = c("CYP450" = "#5e3c99", "Complex I" = "#008837", "UGT" = "#d7191c")) +
  theme_classic() +
  theme(axis.line        = element_blank(),
        axis.ticks.x     = element_blank(),
        axis.title       = element_blank(),
        axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        axis.line.y      = element_blank(),
        legend.position = "top") +
  labs(colour = "Gene set") 


ggsave("output/genemap.pdf", p1, width=10, height=8)
