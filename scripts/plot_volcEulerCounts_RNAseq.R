# =============================================================================
# RNA-seq Euler Diagrams, Volcano Plots, and Overlap Statistics
# =============================================================================
# Generates Euler (Venn-like) diagrams of DEG overlaps across cage populations,
# volcano plots highlighting genes in significant genomic cluster regions, and
# Jaccard overlap statistics. Also plots normalized read counts for selected genes.
#
# Author: Leah Darwin, Yevgeniy Raynes
#
# Output files:
#   output/eulers_all.pdf       : Euler diagrams of all-gene DEG overlaps across
#                                 three contrasts
#   output/CC_pop_vols.pdf      : Volcano plots for control-population pairwise
#                                 comparisons
#   output/RNAseq_fig_cyp.pdf   : Combined Euler diagrams, volcano plots, and
#                                 normalized count plots for CYP genes
#   output/jaccard_overlaps.tex : Jaccard similarity table for DEG set overlaps
# =============================================================================

library(dplyr)
library(ggrastr)
library(ggrepel)
library(DescTools)
library(DESeq2)
library(ggplot2)
library(patchwork)
library(eulerr)
library(stringr)
library(org.Dm.eg.db)
library(knitr)
library(kableExtra)

FDR.Thresh = 0.05
LFC.Thresh = 0.0

##load count data table 
cts = as.matrix(read.csv("data/GDL_count_table.csv", row.names = 1, check.names = FALSE)) # check.names = FALSE prevents X from being added to column names (because they start with numbers)

pops=c("1A3","2A2","2A3")

# Set to TRUE to exclude the 2A2C.R.1 sample from all plots
# (this sample clusters with untreated samples instead of treated samples)
EXCLUDE_OUTLIERS <- TRUE

##column metadata
coldata = read.csv("data/GDL_group_data.csv", row.names = 1) %>%
  mutate(across(everything(), as.factor))

if (EXCLUDE_OUTLIERS) {
  cat("Excluding outlier samples from all analyses...\n")
  samples_to_keep <- colnames(cts) != "2A2.R.C.2" & colnames(cts) != "2A2.C.R.1" & colnames(cts) != "2A3.R.C.1"
  cts <- cts[, samples_to_keep]
  coldata <- coldata[colnames(cts), ]
}

gtf = read.csv(file="data/dmel-all-r6.57.gtf", sep="\t", header = FALSE) %>% 
  dplyr::filter(V3=="gene") %>%
  dplyr::select(V1, V4, V5, V9)%>%
  mutate(gene_id = str_extract(V9, "FBgn[0-9]+")) %>%
  mutate(gene_name = sub(".*gene_symbol ([^;]+);.*", "\\1", V9)) %>%
  dplyr::select(-V9)

colnames(gtf) = c("chr","start","end", "gene_id", "gene_name")

##make DEseq object
dds = DESeqDataSetFromMatrix(countData = cts,
                             colData = coldata,
                             design = ~ Group)
keep = rowSums(counts(dds)) >= 1 #minimal pre-filtering
table(keep)
dds = dds[keep,]
dds = DESeq(dds)

clust_bounds = read.csv("data/clust_bounds.csv")


# ============================================================================
# FUNCTIONS
# ============================================================================


assign_cluster = function(df, bounds) {
  df$cluster = NA_character_
  for (i in seq_len(nrow(bounds))) {
    hit = df$chr == bounds$CHR[i] &
      df$start <= bounds$end[i] &
      df$end   >= bounds$start[i]
    df$cluster[hit] = bounds$cluster[i]
  }
  df[!is.na(df$cluster), ]
}

get_res = function(contrast){
  
  res <- results(dds, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                 contrast=contrast)
  
  res_df = data.frame(res@rownames, res@listData[["pvalue"]],  res@listData[["padj"]], res@listData[["log2FoldChange"]]) 
  
  colnames(res_df) = c("gene","pval","padj","logFC")
  
  res_df = res_df %>%
    dplyr::left_join(gtf, by = join_by(gene==gene_id)) %>%
    na.omit()%>%
    dplyr::filter(padj<FDR.Thresh)
  
  return(res_df)
  
}

get_intersections = function(sets) {
  nms = names(sets)
  all_genes = unique(unlist(sets))
  membership = sapply(sets, function(s) all_genes %in% s)
  rownames(membership) = all_genes
  
  regions = list()
  for (k in seq_along(nms)) {
    combos = combn(nms, k, simplify = FALSE)
    for (combo in combos) {
      in_these  = rowAlls(membership[, combo, drop = FALSE])
      out_those = if (length(combo) < length(nms))
        rowAlls(!membership[, setdiff(nms, combo), drop = FALSE])
      else rep(TRUE, length(all_genes))
      genes = all_genes[in_these & out_those]
      if (length(genes) > 0)
        regions[[paste(combo, collapse = "&")]] = genes
    }
  }
  regions
}

#' Fisher's exact test comparing DEG counts between two plot_volcano res_dfs
#'
#' @param df1      res_df from plot_volcano (contrast 1)
#' @param df2      res_df from plot_volcano (contrast 2)
#' @param label1   label for contrast 1
#' @param label2   label for contrast 2
#' @return list with contingency table, fisher test result, and summary data frame
test_deg_enrichment = function(df1, df2, label1 = "contrast1", label2 = "contrast2") {
  n_deg1    <- sum(df1$significance != "NS")
  n_nondeg1 <- sum(df1$significance == "NS")
  n_deg2    <- sum(df2$significance != "NS")
  n_nondeg2 <- sum(df2$significance == "NS")

  mat <- matrix(
    c(n_deg1, n_nondeg1, n_deg2, n_nondeg2),
    nrow = 2,
    dimnames = list(
      c("DEG", "non-DEG"),
      c(label1, label2)
    )
  )

  ft <- fisher.test(mat)
  gt <- DescTools::GTest(mat)

  summary_df <- data.frame(
    contrast      = c(label1, label2),
    n_deg         = c(n_deg1, n_deg2),
    n_total       = c(nrow(df1), nrow(df2)),
    pct_deg       = round(100 * c(n_deg1/nrow(df1), n_deg2/nrow(df2)), 2),
    odds_ratio    = c(ft$estimate, NA),
    fisher_p      = c(ft$p.value, NA),
    g_stat        = c(gt$statistic, NA),
    g_p           = c(gt$p.value, NA),
    stringsAsFactors = FALSE
  )

  list(table = mat, fisher = ft, gtest = gt, summary = summary_df)
}

#' Calculate 3-way Jaccard similarity for a named list of 3 gene sets
#' J = |A ∩ B ∩ C| / |A ∪ B ∪ C|
#'
#' @param gene_sets Named list of 3 character vectors (e.g. make_euler(...)$gene_sets)
#' @return Data frame with set sizes, intersection, union, and jaccard
calc_jaccard = function(gene_sets) {
  nms  = names(gene_sets)
  sets = lapply(gene_sets, unique)
  n_int = length(Reduce(intersect, sets))
  n_uni = length(Reduce(union,     sets))
  data.frame(
    sets         = paste(nms, collapse = " & "),
    n_set1       = length(sets[[1]]),
    n_set2       = length(sets[[2]]),
    n_set3       = length(sets[[3]]),
    intersection = n_int,
    union        = n_uni,
    jaccard      = round(n_int / n_uni, 3),
    stringsAsFactors = FALSE
  )
}

#' Build euler plots and intersection data for a single contrast
#'
#' @param g1          First group suffix (e.g. "C.C")
#' @param g2          Second group suffix (e.g. "R.R")
#' @param populations Character vector of population prefixes (default: pops)
#' @param label       Optional label for plot titles; defaults to "g1 v g2"
#' @return List with euler plot panels (all genes and cluster-filtered) and intersect_df
make_euler = function(g1, g2, populations = pops, label = NULL) {
  if (is.null(label)) label = paste0(g1, " v ", g2)

  dfs = lapply(populations, function(p)
    get_res(c("Group", paste0(p, ".", g1), paste0(p, ".", g2)))
  )
  names(dfs) = populations

  gene_sets = lapply(dfs, function(x) unique(x$gene))
  
  dfs_clust      = lapply(dfs, assign_cluster, bounds = clust_bounds)
  gene_sets_clust = lapply(dfs_clust, function(x) unique(x$gene))
  
  fit       = euler(gene_sets)
  fit_clust = euler(gene_sets_clust)
  
  intersections_all   = get_intersections(gene_sets)
  intersections_clust = get_intersections(gene_sets_clust)
  
  build_intersect_df = function(intersections, meta_dfs, extra_cols) {
    gene_info = bind_rows(meta_dfs) %>%
      distinct(gene, .keep_all = TRUE) %>%
      dplyr::select(gene, gene_name, chr, start, end, any_of(extra_cols))
    bind_rows(
      lapply(names(intersections), function(region) {
        data.frame(intersection = region,
                   gene = intersections[[region]],
                   stringsAsFactors = FALSE)
      })
    ) %>%
      left_join(gene_info, by = "gene")
  }
  
  intersect_df_all   = build_intersect_df(intersections_all,   dfs,       c())
  intersect_df_clust = build_intersect_df(intersections_clust, dfs_clust, "cluster")
  
  # These overlap cleanly at alpha = 0.4-0.5
  pop_cols = c("#386CB0", "#BEAED4", "#7FC97F") # green

  list(
    p_all        = plot(fit,
                        quantities = list(fontsize=8),
                        labels = list(fontsize=10),
                        adjust_labels = TRUE,
                        fill = pop_cols,
                        alpha = 0.4,
                        ),
    p_clust      = plot(fit_clust,
                        quantities = list(fontsize=8),
                        labels = list(fontsize=10),
                        adjust_labels = TRUE,
                        fill = pop_cols,
                        alpha = 0.4
                        ),
    intersect_df_all   = intersect_df_all,
    intersect_df_clust = intersect_df_clust,
    gene_sets          = gene_sets,
    gene_sets_clust    = gene_sets_clust
  )
}

plot_volcano = function(contrast, title, pCutoff = FDR.Thresh, FCcutoff = LFC.Thresh, label_genes = NULL, cluster_only = TRUE){
  
  res <- results(dds, alpha = FDR.Thresh, lfcThreshold = LFC.Thresh,
                 contrast=contrast)
  summary(res)
  
  # Convert to data frame for ggplot
  res_df <- as.data.frame(res) %>%
    mutate(gene = rownames(res)) %>%
    filter(!is.na(padj)) %>%
    dplyr::left_join(gtf, by = join_by(gene==gene_id)) %>%
    na.omit()

  res_df <- res_df %>% mutate(neg_log10_padj = -log10(padj))

  if (cluster_only) {
    # Genes that overlap a cluster region get full Up/Down colour;
    # significant genes outside clusters get a darker-grey "Sig" category
    cluster_genes <- assign_cluster(res_df, clust_bounds)$gene
    res_df <- res_df %>%
      mutate(significance = case_when(
        padj < pCutoff & log2FoldChange > FCcutoff  & gene %in% cluster_genes ~ "Up",
        padj < pCutoff & log2FoldChange < -FCcutoff & gene %in% cluster_genes ~ "Down",
        padj < pCutoff ~ "Sig",
        TRUE ~ "NS"
      ))
  } else {
    res_df <- res_df %>%
      mutate(significance = case_when(
        padj < pCutoff & log2FoldChange > FCcutoff  ~ "Up",
        padj < pCutoff & log2FoldChange < -FCcutoff ~ "Down",
        TRUE ~ "NS"
      ))
  }

  # Define colors  ("Sig" = significant but outside cluster bounds)
  sig_colors <- c("Up" = "#ca0020", "Down" = "#0571b0", "Sig" = "grey40", "NS" = "grey80")

  # Count up and down regulated genes (cluster-overlapping only)
  n_up <- sum(res_df$significance == "Up")
  n_down <- sum(res_df$significance == "Down")
  
  # Get axis limits for annotation placement
  x_range <- range(res_df$log2FoldChange)
  y_max <- max(res_df$neg_log10_padj)
  
  p = ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj, color = significance)) +
    geom_point_rast(alpha = 0.6, size = 1.5, raster.dpi = 600) +
    scale_color_manual(values = sig_colors) +
    geom_vline(xintercept = c(-FCcutoff, FCcutoff), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(pCutoff), linetype = "dashed", color = "grey40") +
    annotate("text", x = x_range[2], y = y_max, label = n_up,
             color = "#e41a1c", hjust = 1, vjust = 1, size = 3.5) +
    annotate("text", x = x_range[1], y = y_max, label = n_down,
             color = "#377eb8", hjust = 0, vjust = 1, size = 3.5) +
    labs(
      title = title,
      x = "log2 Fold Change",
      y = "-log10(adjusted p-value)",
      color = "Significance"
    ) +
    theme_classic() +
    theme(
      #plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )

  if (!is.null(label_genes)) {
    label_df <- res_df %>% filter(gene_name %in% label_genes)
    p <- p + geom_text_repel(
      data = label_df,
      aes(label = gene_name),
      color = "black",
      size = 3,
      segment.color = "grey40",
      segment.linewidth = 0.4,
      box.padding = 0.4,
      max.overlaps = Inf
    )
  }

  list(plot = p, res_df = res_df)
}

## Plot normalized read counts for a given gene across populations, faceted by condition
## gene:           gene ID (e.g. "FBgn0001234") or gene symbol matching gtf$gene_name
## group_suffixes: character vector of group suffixes, e.g. c("R.R", "C.C")
## populations:    character vector of population prefixes (default: 1A3, 2A2, 2A3)
plot_norm_counts = function(gene, group_suffixes = c("R.R", "C.C"), populations = pops) {
  
  # Resolve gene symbol -> gene ID if needed
  if (!gene %in% rownames(dds)) {
    hit <- gtf$gene_id[gtf$gene_name == gene]
    if (length(hit) == 0) stop("Gene '", gene, "' not found in rownames(dds) or gtf gene_name.")
    gene <- hit[1]
  }
  
  # Pull normalized counts for this gene
  norm_counts <- counts(dds, normalized = TRUE)[gene, ]
  
  # Build one data frame across all requested conditions
  df <- lapply(group_suffixes, function(gs) {
    target_groups <- paste0(populations, ".", gs)
    as.data.frame(colData(dds)) %>%
      mutate(
        sample     = rownames(colData(dds)),
        norm_count = norm_counts[rownames(colData(dds))]
      ) %>%
      filter(Group %in% target_groups) %>%
      mutate(
        pop       = sub(paste0("\\.", gs, "$"), "", as.character(Group)),
        condition = gs
      )
  }) %>% bind_rows()
  
  gene_label <- gtf$gene_name[gtf$gene_id == gene]
  gene_label <- if (length(gene_label) > 0) gene_label[1] else gene
  
  ggplot(df, aes(x = pop, y = norm_count, color = pop)) +
    geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.4,
                 color = "black", linewidth = 0.5) +
    scale_color_manual(values = setNames(c("#386CB0", "#BEAED4", "#7FC97F"), populations)) +
    facet_wrap(~ condition, nrow=1) +
    labs(title = gene_label,
         x = "Population",
         y = "Normalized\nread count",
         color = "Population") +
    theme_linedraw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color = "black")
    )
}

res_CCvRR = make_euler("C.C", "R.R")
res_CCvRC = make_euler("C.C", "R.C")
res_CCvCR = make_euler("C.C", "C.R")

eulers = wrap_elements(res_CCvCR$p_all) + wrap_elements(res_CCvRR$p_clust) + wrap_elements(res_CCvRC$p_clust)+         
                         plot_layout(nrow=2)
eulers

eulers_all = wrap_elements(res_CCvCR$p_all) + ggtitle('C.C v C.R') +
  wrap_elements(res_CCvRR$p_all) + ggtitle('C.C v R.R') + 
  wrap_elements(res_CCvRC$p_all)+  ggtitle('C.C v R.C') +        
  plot_layout(nrow=1)
ggsave("output/eulers_all.pdf", eulers_all, width=6, height=3)

eulers_all


fmt_pval = function(p) ifelse(p < 0.0001, "<0.0001", formatC(p, format = "g", digits = 3))

# ============================================================================
# KABLE TABLES: DEG COUNTS AND JACCARD OVERLAPS
# ============================================================================

# DEG counts per contrast and population
euler_list = list(
  "C.C v C.R" = res_CCvCR,
  "C.C v R.R" = res_CCvRR,
  "C.C v R.C" = res_CCvRC
)

# Jaccard overlaps (all genes and cluster-filtered)
make_jaccard_row = function(gs, contrast_name, subset_name) {
  pop_names = names(gs)
  sets      = lapply(gs, unique)
  n_int     = length(Reduce(intersect, sets))
  n_uni     = length(Reduce(union,     sets))
  row       = data.frame(
    Contrast     = contrast_name,
    Subset       = subset_name,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(pop_names)) row[[pop_names[i]]] = length(sets[[i]])
  row$Intersection = n_int
  row$Union        = n_uni
  row$Jaccard      = round(n_int / n_uni, 3)
  row
}

jaccard_df = bind_rows(lapply(names(euler_list), function(contrast_name) {
  bind_rows(
    make_jaccard_row(euler_list[[contrast_name]]$gene_sets,       contrast_name, "all genes"),
    make_jaccard_row(euler_list[[contrast_name]]$gene_sets_clust, contrast_name, "cluster genes")
  )
}))

kable(jaccard_df,
      format = "latex",
      booktabs = TRUE,
      caption = "Three-way Jaccard similarity of DEG sets across populations for each C.C contrast.",
      label = "jaccard_overlaps",
      align = c("l", "l", rep("r", ncol(jaccard_df) - 2))) %>%
  kable_styling(latex_options = c("hold_position"))%>%
  save_kable("output/jaccard_overlaps.tex")

res1 = plot_volcano(c("Group","2A2.R.R","2A3.R.R"), "2A2 RR v 2A3 RR")
res2 = plot_volcano(c("Group","2A2.R.R","1A3.R.R"), "2A2 RR v 1A3 RR")
res3 = plot_volcano(c("Group","2A3.R.R","1A3.R.R"), "2A3 RR v 1A3 RR")
vols = res2$plot + res3$plot + res1$plot +plot_layout(guides="collect")&theme(legend.position = "bottom")
vols

res4 = plot_volcano(c("Group","2A2.C.C","2A3.C.C"), "2A2 CC v 2A3 CC", cluster_only = FALSE)
res5 = plot_volcano(c("Group","2A2.C.C","1A3.C.C"), "2A2 CC v 1A3 CC", cluster_only = FALSE)
res6 = plot_volcano(c("Group","2A3.C.C","1A3.C.C"), "2A3 CC v 1A3 CC", cluster_only = FALSE)
vols_C = res5$plot + res6$plot + res4$plot +plot_layout(guides="collect")&theme(legend.position = "bottom")
vols_C
ggsave("output/CC_pop_vols.pdf",vols_C, width=6, height=4)

test_deg_enrichment(res1$res_df, res4$res_df, "2A2v2A3 RR", "2A2v2A3 CC")
test_deg_enrichment(res2$res_df, res5$res_df, "2A2v1A3 RR", "2A2v1A3 CC")
test_deg_enrichment(res3$res_df, res6$res_df, "2A3v1A3 RR", "2A3v1A3 CC")

counts2 = res2$res_df %>%filter(significance!="NS")
paste0(sum(counts2$significance == "Up") + sum(counts2$significance == "Down"),"/", nrow(counts2))

counts1 = res1$res_df %>%filter(significance!="NS")
paste0(sum(counts1$significance == "Up") + sum(counts1$significance == "Down"),"/", nrow(counts1))

counts3 = res3$res_df %>%filter(significance!="NS")
paste0(sum(counts3$significance == "Up") + sum(counts3$significance == "Down"),"/", nrow(counts3))



p2_nc = plot_norm_counts("Cyp6a20", group_suffixes = c("R.R", "C.C","C.R","R.C"))
p3_nc = plot_norm_counts("Cyp49a1", group_suffixes = c("R.R", "C.C","C.R","R.C"))
p3_nc

ncs = p2_nc + p3_nc
ncs

final = (wrap_plots(eulers)+wrap_plots(vols)+plot_layout(widths=c(1,2)))/wrap_plots(ncs)+ plot_layout(heights=c(2,1)) + plot_annotation("a")

ggsave("output/RNAseq_fig_cyp.pdf", final, width=10, height=6)
