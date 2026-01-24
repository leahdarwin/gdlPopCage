library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

df_selec = read.csv("../data/mitofreqs.csv")
df_preselec = read.csv("../data/mitofreqs_preselec.csv")

color_palette <- c(
  "Beijing"   = "#18458e",
  "Zimbabwe"  = "#4fc2ba",
  "D. yakuba"  = "#ff7f00",
  "D. simulans"= "#e31a1c"
)

df_long = df_selec %>%
  pivot_longer(
    cols = -c(cage, treatment, gen),   
    names_to = "mt",
    values_to = "freq"
  )%>% # assign names in order
  mutate(mt_group = case_when(grepl("^B",mt) ~ "Beijing",
                                   grepl("^Z", mt) ~ "Zimbabwe",
                                   grepl("^Y", mt) ~ "D. yakuba",
                                   grepl("^S", mt) ~ "D. simulans"))%>%
  mutate(
    mt_group = factor(
      mt_group,
      levels = c( "D. yakuba","D. simulans","Beijing", "Zimbabwe" )  # your desired order
    )
  )

p1 = ggplot(df_long %>% filter(treatment=="C"), aes(x = paste("F",gen,sep=""), y = freq, fill = mt_group)) +
  geom_col(width = 1, linewidth = 0.15, color="white") +
  facet_wrap(~cage, nrow=2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = color_palette) +
  labs(x = "Generation", y = "Estimated haplotype frequency") + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+scale_y_continuous(
    breaks = c(0, 0.5, 1),
    expand = c(0, 0)) + theme(panel.spacing.y = unit(1, "lines"))

p2 = ggplot(df_long %>% filter(treatment=="R"), aes(x = paste("F",gen,sep=""), y = freq, fill = mt_group)) +
  geom_col(width = 1, linewidth = 0.15, color="white") +
  facet_wrap(~cage, nrow=2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = color_palette) +
  labs(x = "Generation", y = "Estimated haplotype frequency") + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+scale_y_continuous(
    breaks = c(0, 0.5, 1),
    expand = c(0, 0)) + theme(panel.spacing.y = unit(1, "lines"))



p1/p2+ plot_layout(guides = "collect") & theme(legend.position = 'bottom')&
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')', 
                  theme = theme(plot.margin = margin(1, 1, 1, 1))) 
