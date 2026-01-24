library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

##rotenone was F73 control was F75
df_F70 = read.csv("../data/pupalcount_July2025.csv") %>%
  group_by(Cage, Treatment, Replicate) %>%
  summarise(TSE = weighted.mean(TSE, w = Count)) %>%
  separate_wider_delim(Cage, ".", names = c("cage", "cage_treatment"))%>%
  mutate(gen = case_when(cage_treatment == "R" ~ "F73",
                         cage_treatment == "C" ~ "F75"))

##F48 rotenone F50 control
df_F50 = read.csv("../data/experimentalResults/eggLay_May2024.csv")%>%
  group_by(cage, cage_treatment, vial, treatment) %>%
  summarise(TSE = weighted.mean(tse, w = count))  %>%
  mutate(gen = case_when(cage_treatment == "R" ~ "F48",
                         cage_treatment == "C" ~ "F50")) %>%
  rename(Treatment=treatment, Replicate = vial)

df = rbind(df_F50, df_F70) %>%
  mutate(CTT = paste(cage_treatment, Treatment)) %>%
  group_by(cage, cage_treatment, Treatment, gen, CTT) %>%
  summarise(TSE = mean(TSE))



ggplot(df, aes(x = Treatment, y = TSE, color = Treatment, fill = Treatment, group = interaction(gen, Treatment))) +
  # Use a slight fill with alpha so it's not overwhelming, and fatten the lines
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5, alpha = 0.1, lwd = 0.7, staplewidth = 0.4) +
  # Manual colors for both border (color) and interior (fill)
  scale_color_manual(values = c("C" = "grey14", "R" = "red")) +
  scale_fill_manual(values = c("C" = "grey14", "R" = "red")) +
  facet_wrap(~cage_treatment, scales = "free_x") +
  # The "Publication" look
  theme_linedraw() +
  theme(
    strip.background = element_rect(fill = "white"), # Light facet headers
    strip.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),             # Remove distracting minor lines
    panel.grid.major.x = element_blank(),           # Remove vertical lines for cleaner grouping
    axis.title = element_text(),
    legend.position = "top",                        # Legend at top saves horizontal space
    legend.title = element_text()
  ) +
  labs(
    x = "Generation",
    y = "Development time (hours)",
    color = "Treatment Group:",
    fill = "Treatment Group:"
  )


