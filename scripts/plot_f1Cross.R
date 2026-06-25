# =============================================================================
# F1 Cross Development Time Analysis
# =============================================================================
# Analyzes development time in reciprocal F1 crosses between control and
# rotenone cage populations, fitting a linear mixed model and computing pairwise
# contrasts by treatment.
#
# Author: Leah Darwin
#
# Output files:
#   output/F1_cross.pdf          : Boxplot of development time by cross type
#                                  and treatment
#   output/anova_f1cross.tex     : ANOVA table for the F1 cross mixed model
#   output/contrasts_f1cross.tex : Pairwise contrasts by treatment
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(knitr)
library(kableExtra)

df = read.csv("data/rcCross_F1_pupalcounts.csv") %>%
  summarise(TSE = weighted.mean(tse, w = count), .by=c(vial, treatment, cage, f_m)) %>%
  separate_wider_delim(f_m, delim = "_", names = c("f", "m"), cols_remove =F)

colors = c("C_C"="black","R_R"="red", "C_R"="#910028", "R_C"="#580d22")

p1 = ggplot(df, aes(x=cage, y=TSE, fill=f_m)) +
  geom_boxplot(alpha=0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), 
             alpha = 0.5, aes(color=f_m)) +
  facet_wrap(~treatment, labeller = as_labeller(c("C"="Control Treated","R"="Rotenone Treated"))) +
  theme_linedraw() +
  labs(y="Development time (hours)",
       x="Population",
       fill="F1 Cross (F_M)")+
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(),
    legend.position = "top",
    legend.title = element_text()
  ) +
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors) +
  guides(color = "none")

ggsave("output/F1_cross.pdf", p1, width=6, height=4)


lin = lmer(TSE~f_m*treatment+(1|cage), data=df)

emm = emmeans(lin, specs = ~ f_m*treatment)
contrasts_res = contrast(emm, method = "pairwise", by = "treatment", adjust="bonferroni")

# ============================================================================
# KABLE TABLES: ANOVA AND CONTRASTS
# ============================================================================

var_labels = c(f_m = "cross", treatment = "treatment")

rename_terms = function(x) {
  for (old in names(var_labels)) {
    x = gsub(old, var_labels[old], x, fixed = TRUE)
  }
  x
}

fmt_pval = function(p) ifelse(p < 0.0001, "<0.0001", formatC(p, format = "g", digits = 3))

# ANOVA table
anova_df = as.data.frame(anova(lin)) %>%
  tibble::rownames_to_column("Term") %>%
  mutate(Term = rename_terms(Term),
         across(c(`Sum Sq`, `Mean Sq`, `F value`), ~ round(.x, 3)),
         DenDF = round(DenDF, 1),
         `Pr(>F)` = fmt_pval(`Pr(>F)`))

kable(anova_df,
      format = "latex",
      booktabs = TRUE,
      caption = "ANOVA results for linear mixed model of F1 cross development time.",
      label = "anova_f1cross",
      align = c("l", rep("r", ncol(anova_df) - 1))) %>%
  kable_styling(latex_options = c("hold_position")) %>%
  save_kable("output/anova_f1cross.tex")

# Contrasts table
contrasts_df = as.data.frame(contrasts_res) %>%
  mutate(contrast = rename_terms(as.character(contrast)),
         across(c(estimate, SE, df, t.ratio), ~ round(.x, 3)),
         p.value = fmt_pval(p.value))

kable(contrasts_df,
      format = "latex",
      booktabs = TRUE,
      caption = "Pairwise contrasts for F1 cross development time (Bonferroni-adjusted), by treatment.",
      label = "contrasts_f1cross",
      align = c("l", rep("r", ncol(contrasts_df) - 1))) %>%
  kable_styling(latex_options = c("hold_position", "scale_down")) %>%
  save_kable("output/contrasts_f1cross.tex")
