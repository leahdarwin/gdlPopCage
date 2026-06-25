# =============================================================================
# Development Time Analysis
# =============================================================================
# Analyzes egg-to-adult development time (TSE) for rotenone and control
# treatments across two experimental generations (F48/F50 and F73/F75), fitting
# a linear mixed model and computing estimated marginal means and pairwise
# contrasts.
#
# Author: Leah Darwin
#
# Output files:
#   output/devtime.pdf            : Boxplot of development time by treatment
#   output/anova_devtime.tex      : ANOVA table for the linear mixed model
#   output/contrasts_devtime.tex  : Pairwise contrasts by cage treatment and origin
#   output/emm_devtime.tex        : Estimated marginal means table
#   output/contrasts_devtime2.tex : Additional contrasts pooled across generations
# =============================================================================

# ============================================================================
# LOAD PACKAGES
# ============================================================================

library(ggplot2)  # For creating plots
library(dplyr)    # For data manipulation
library(tidyr)    # For data reshaping (separate_wider_delim)
library(patchwork) # For combining plots
library(lme4)
library(lmerTest)
library(emmeans)
library(knitr)
library(kableExtra)

# ============================================================================
# DATA IMPORT AND PROCESSING
# ============================================================================

# F70 Generation Data
# Rotenone treatment was at F73, control was at F75
df_F70 = read.csv("data/pupalcount_July2025.csv") %>%
  group_by(Cage, Treatment, Replicate) %>%
  # Calculate weighted mean TSE (Time Since Egg) based on Count
  summarise(TSE = weighted.mean(TSE, w = Count)) %>%
  # Split Cage column into cage and cage_treatment
  separate_wider_delim(Cage, ".", names = c("cage", "cage_treatment"))%>%
  # Assign generation labels based on treatment
  mutate(gen = case_when(cage_treatment == "R" ~ "F73",
                         cage_treatment == "C" ~ "F75"))

# F50 Generation Data
# Rotenone treatment was at F48, control was at F50
df_F50 = read.csv("data/experimentalResults/eggLay_May2024.csv")%>%
  group_by(cage, cage_treatment, vial, treatment) %>%
  # Calculate weighted mean tse (lowercase) based on count (lowercase)
  summarise(TSE = weighted.mean(tse, w = count))  %>%
  # Assign generation labels based on treatment
  mutate(gen = case_when(cage_treatment == "R" ~ "F48",
                         cage_treatment == "C" ~ "F50")) %>%
  # Rename columns to match F70 data structure
  rename(Treatment=treatment, Replicate = vial)

# ============================================================================
# COMBINE AND AGGREGATE DATA
# ============================================================================

# Combine both generation datasets
df = rbind(df_F50, df_F70) %>%
  # Create combined treatment identifier
  mutate(CTT = paste(cage_treatment, Treatment)) %>%
  # Calculate mean TSE for each cage-treatment-generation combination
  group_by(cage, cage_treatment, Treatment, gen, CTT) %>%
  summarise(TSE = mean(TSE)) %>%
  mutate(cage_orig = case_when(grepl("^1", cage)~"1",
                               grepl("^2", cage)~"2"))

# ============================================================================
# PLOT DEVELOPMENT TIME
# ============================================================================

# Create boxplot of development time by treatment and cage treatment
# Faceted by cage_treatment (R or C from original cage population)
p1 = ggplot(df, aes(x = Treatment, y = TSE, color = Treatment, fill = Treatment, group = interaction(gen, Treatment))) +
  # Boxplot with transparent fill to see individual points
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5, alpha = 0.3, lwd = 0.7, staplewidth = 0.4) +
  # Define colors: Control = grey/black, Rotenone = red
  scale_color_manual(values = c("C" = "grey14", "R" = "red")) +
  scale_fill_manual(values = c("C" = "grey14", "R" = "red")) +
  # 4 side-by-side panels: cage_orig as outer grouping, cage_treatment nested inside
  facet_grid(. ~  cage_treatment + cage_orig) +
  # Publication-style theme
  theme_linedraw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(),
    legend.position = "top",
    legend.title = element_text()
  ) +
  labs(
    x = "Generation",
    y = "Development time (hours)",
    color = "Treatment Exposure:",
    fill = "Treatment Exposure:"
  )

 df_lm = rbind(df_F50, df_F70)%>%
   mutate(cage_orig = case_when(grepl("^1", cage)~"1",
                                grepl("^2", cage)~"2")) %>%
   mutate(gen = case_when(gen%in%c("F73","F75")~"F7*",
                          gen%in%c("F48","F50")~"F5*")) 

linm = lmer(TSE~cage_orig+gen+Treatment+cage_treatment+cage_orig:gen+Treatment:gen+Treatment:gen:cage_treatment+cage_treatment:gen+cage_orig:Treatment + Treatment:cage_treatment+ (1|cage_treatment:cage),
            data=df_lm)
anova(linm)

emm = emmeans(linm, specs = ~ cage_orig*Treatment*cage_treatment*gen)
contrasts_res = contrast(emm, method = "pairwise", by = c("cage_treatment","cage_orig"), adjust="bonferroni")

emm2 = emmeans(linm, specs = ~ Treatment*cage_treatment)
emm2
contrasts_res2 = contrast(emm2, method = "pairwise", by = c("cage_treatment"), adjust="bonferroni")
contrasts_res

emm3 = emmeans(linm, specs = ~ Treatment*cage_treatment*cage_orig)
emm3
contrasts_res3 = contrast(emm3, method = "pairwise", by = c("cage_treatment","cage_orig"), adjust="bonferroni")
contrasts_res


# ============================================================================
# KABLE TABLES: ANOVA AND CONTRASTS
# ============================================================================

# Variable name mapping for display
var_labels = c(
  cage_treatment = "environment",
  cage_orig      = "history",
  Treatment      = "treatment",
  gen            = "generation"
)

rename_terms = function(x) {
  for (old in names(var_labels)) {
    x = gsub(old, var_labels[old], x, fixed = TRUE)
  }
  x
}

fmt_pval = function(p) ifelse(p < 0.0001, "<0.0001", formatC(p, format = "g", digits = 3))

# ANOVA table
anova_df = as.data.frame(anova(linm)) %>%
  tibble::rownames_to_column("Term") %>%
  mutate(Term = rename_terms(Term),
         across(c(`Sum Sq`, `Mean Sq`, `F value`), ~ round(.x, 3)),
         DenDF = round(DenDF, 1),
         `Pr(>F)` = fmt_pval(`Pr(>F)`))

kable(anova_df,
      format = "latex",
      booktabs = TRUE,
      caption = "ANOVA results for linear mixed model of development time.",
      label = "anova_devtime",
      align = c("l", rep("r", ncol(anova_df) - 1))) %>%
  kable_styling(latex_options = c("hold_position")) %>%
  save_kable("output/anova_devtime.tex")

# Contrasts table
contrasts_df = as.data.frame(contrasts_res) %>%
  rename(environment = cage_treatment, history = cage_orig) %>%
  mutate(contrast = rename_terms(as.character(contrast)),
         across(c(estimate, SE, df, t.ratio), ~ round(.x, 3)),
         p.value = fmt_pval(p.value))

kable(contrasts_df,
      format = "latex",
      booktabs = TRUE,
      caption = "Pairwise contrasts for development time (Bonferroni-adjusted), by cage treatment and cage origin.",
      label = "contrasts_devtime",
      align = c("l", rep("r", ncol(contrasts_df) - 1))) %>%
  kable_styling(latex_options = c("hold_position", "scale_down")) %>%
  save_kable("output/contrasts_devtime.tex")

# Combined emmeans table (emm2 + emm3)
emm2_df = as.data.frame(emm2) %>%
  rename(environment = cage_treatment) %>%
  mutate(across(c(emmean, SE, df, lower.CL, upper.CL), ~ round(.x, 3)))

emm3_df = as.data.frame(emm3) %>%
  rename(environment = cage_treatment, history = cage_orig) %>%
  mutate(across(c(emmean, SE, df, lower.CL, upper.CL), ~ round(.x, 3)))

emm_combined = bind_rows(emm2_df, emm3_df) %>%
  select(Treatment, environment, history, everything())

kable(emm_combined,
      format = "latex",
      booktabs = TRUE,
      caption = "Estimated marginal means for development time, pooled across generations (upper rows: Treatment × environment; lower rows: Treatment × environment × history).",
      label = "emm_devtime",
      align = c("l", rep("r", ncol(emm_combined) - 1))) %>%
  kable_styling(latex_options = c("hold_position", "scale_down"))%>%
  save_kable("output/emm_devtime.tex")

# Combined contrasts table (contrasts_res2 + contrasts_res3)
contrasts2_df = as.data.frame(contrasts_res2) %>%
  rename(environment = cage_treatment) %>%
  mutate(contrast = rename_terms(as.character(contrast)),
         across(c(estimate, SE, df, t.ratio), ~ round(.x, 3)),
         p.value = fmt_pval(p.value))

contrasts3_df = as.data.frame(contrasts_res3) %>%
  rename(environment = cage_treatment, history = cage_orig) %>%
  mutate(contrast = rename_terms(as.character(contrast)),
         across(c(estimate, SE, df, t.ratio), ~ round(.x, 3)),
         p.value = fmt_pval(p.value))

contrasts_combined = bind_rows(contrasts2_df, contrasts3_df) %>%
  select(contrast, environment, history, everything())

kable(contrasts_combined,
      format = "latex",
      booktabs = TRUE,
      caption = "Pairwise contrasts for development time (Bonferroni-adjusted), pooled across generations (upper rows: by environment; lower rows: by environment and history).",
      label = "contrasts_devtime2",
      align = c("l", rep("r", ncol(contrasts_combined) - 1))) %>%
  kable_styling(latex_options = c("hold_position", "scale_down")) %>%
  save_kable("output/contrasts_devtime2.tex")

ggsave("output/devtime.pdf", p1, width = 6, height = 3)
