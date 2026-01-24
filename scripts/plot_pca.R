library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggrepel)
library(patchwork)
library(plotly)

afs = read.csv("../data/pca/joined.sync.MAF01.frq",sep = '\t')  %>%
  filter(if_all(4:123, is.finite)) %>%
  na.omit()

labels = read.csv("../data/pca/sync_labs.csv") %>%
  mutate(orig4 = case_when(startsWith(cage, "1A") ~ "1A",
                           startsWith(cage, "1B") ~ "1B", 
                           startsWith(cage, "2A") ~ "2A", 
                           startsWith(cage, "2B") ~ "2B")) %>%
  mutate(orig2 = case_when(startsWith(cage, "1") ~ "1",
                           startsWith(cage, "2") ~ "2")) %>%
  mutate(gen = factor(gen))

pc = prcomp(t(as.matrix(afs[,4:123])), center = T, scale. = F)

rm(afs)
gc()

pca_df = as.data.frame(pc$x)
pca_df = cbind(pca_df, labels)

pca_var <- pc$sdev^2
pca_var_per <- round(pca_var / sum(pca_var) * 100, 1)

custom_colors = c("C" = "black", "R" = "red")
custom_labels = c("C" = "Control", "R" = "Rotenone")
custom_colors_cage = c("1" = "#9ecae1", "2" = "#3182bd")

pca_df$group = interaction(pca_df$treatment, pca_df$gen)

centroids <- pca_df %>%
  filter(!gen %in% c(22, 40)) %>%
  group_by(treatment, gen, group) %>%
  summarize(PC1 = mean(PC1), PC2 = mean(PC2), PC3 = mean(PC3), .groups = 'drop')

var_values <- pc$sdev^2
var_pct <- var_values / sum(var_values) * 100

# 2. Create a data frame for plotting (first 10 PCs is usually enough)
scree_data <- data.frame(
  PC = 1:10,
  Variance = var_pct[1:10]
)

# Ensure the PCs stay in order 1-10 on the x-axis
scree_data$PC <- factor(scree_data$PC, levels = scree_data$PC)


p1 = ggplot(pca_df, aes(x = PC1, y = PC3, color = treatment)) +
  geom_point(size=2) +
  theme_classic()+
  stat_ellipse(data = filter(pca_df, !gen %in% c(22, 40)), 
               aes(group = group), 
               linetype = 1) +
   geom_label_repel(data = centroids, aes(label = paste("F",gen,sep=""), color=treatment),
                    box.padding = 0.5, segment.color = 'grey50') +
  labs(
    x = paste0("PC1 (", pca_var_per[1], "%)"),
    y = paste0("PC3 (", pca_var_per[3], "%)")
  )  +
  scale_color_manual(values=custom_colors)+ 
  theme(legend.position = "none")

p1

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

p3 = ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_col(fill = "grey15") +           # Bar plot
  labs( y = "Variance \nExplained (%)", 
       x = "Principal Component") +
  theme_classic()

p3

(p1 | (p2 / p3)+plot_layout(heights = c(2,1))) + 
  plot_layout(widths = c(2, 1))+ 
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')')

