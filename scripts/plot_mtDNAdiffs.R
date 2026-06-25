# =============================================================================
# Mitochondrial Haplotype Ordination and Treatment Difference Heatmaps
# =============================================================================
# Ordination (PCo) of mitochondrial haplotype frequencies using CLR-transformed
# composition distances, with PERMANOVA tests for treatment and generation
# effects. Also plots per-cage heatmaps of rotenone vs control haplotype
# proportion differences at each generation.
#
# Author: Leah Darwin
#
# Output files:
#   output/pco_mitos.pdf : PCo1 and PCo2 vs generation, colored by treatment
#   output/R-C_mts.pdf   : Heatmaps of (Rotenone - Control) haplotype frequency
#                          differences per cage at F20, F40, F50, and F60
# =============================================================================

library(DirichletReg)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

set.seed(99)

df = read.csv("data/mitofreqs.csv") %>%
  rename(yak=Yak.Ore.A, siI=Sm21.Ore.A)
df_preselec = read.csv("data/mitofreqs_preselec.csv")

df$Y = DR_data(df[,4:25])

clr_dat = as.matrix(compositions::clr(df$Y))

df$treatcage = paste(df$cage,df$treatment)

m1 = adonis2(clr_dat ~ treatment * gen,
        data = df,
        permutations = how(nperm = 10000,
                           blocks = df$cage,
                           plots = Plots(strata = df$treatcage, type="free"),
                           within = Within(type="series")),
        method = "euclidean",
        by = "margin")
m1


dist_mat = dist(clr_dat, method = "euclidean")

bd = betadisper(dist_mat , interaction(df$treatment, df$gen))
permutest(bd, pairwise = TRUE)

pco = capscale(clr_dat ~ 1, distance = "euclidean")
scores_df = data.frame(
  scores(pco, display = "sites"),
  treatment = df$treatment,
  gen = df$gen,
  cage = df$cage
)

pco_eig = eigenvals(pco)
pco_var = round(100 * pco_eig / sum(pco_eig), 1)

pco1 = ggplot(scores_df, aes(x=gen, y=MDS1, colour = treatment)) +
  geom_jitter(width=0.9) +
  theme_bw() +
  scale_x_continuous(breaks = c(20, 22, 40, 50, 60)) +
  scale_color_manual(values=c("C"="black","R"="red"))+ theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  labs(x="Generation", y=paste0("PCo1 (", pco_var[1], "%)"))

pco2 = ggplot(scores_df, aes(x=gen, y=MDS2, colour = treatment)) +
  geom_jitter(width=0.9) +
  theme_bw() +
  scale_x_continuous(breaks = c(20, 22, 40, 50, 60)) +
  scale_color_manual(values=c("C"="black","R"="red"))+ theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x="Generation", y=paste0("PCo2 (", pco_var[2], "%)"))

pc = pco1/pco2 + plot_layout(guides="collect",axes="collect")
ggsave("output/pco_mitos.pdf", width=5, height=3)

# ============================================================================
# HEATMAP: R - C PROPORTION DIFFERENCE PER CAGE AT A GIVEN GENERATION
# ============================================================================

haplo_names = colnames(df)[4:25]

plot_trt_diff_heat = function(gen_i) {
  diff_df = df %>%
    filter(gen == gen_i) %>%
    select(cage, treatment, all_of(haplo_names)) %>%
    pivot_longer(cols = all_of(haplo_names), names_to = "haplotype", values_to = "proportion") %>%
    pivot_wider(names_from = treatment, values_from = proportion) %>%
    mutate(diff = R - C)

  ggplot(diff_df, aes(x = cage, y = haplotype, fill = diff)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = round(diff, 2)), size = 2.5) +
    scale_fill_gradient2(low = "grey60", mid = "white", high = "red",
                         midpoint = 0, limits = c(-0.5, 0.5), name = "R \u2212 C") +
    labs(x = "Cage", y = "Haplotype",
         title = paste0("Treatment difference at generation ", gen_i)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

f60 = plot_trt_diff_heat(60)
f20 = plot_trt_diff_heat(20)
f40 = plot_trt_diff_heat(40)
f50 = plot_trt_diff_heat(50)

gens = f20 + f40 + f50 + f60 + plot_layout(nrow=2, guides="collect")  
ggsave("output/R-C_mts.pdf", gens, width=10, height=8)

ggplot(df, aes(x=gen, y=B52, group=interaction(cage,treatment), color=treatment)) +
  geom_line()+
  theme_bw() +
  scale_x_continuous(breaks = c(20, 22, 40, 50, 60)) +
  scale_color_manual(values=c("C"="black","R"="red"))+ theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x="Generation")
