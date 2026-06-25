# =============================================================================
# RNA-seq DEG Overlap Permutation Test
# =============================================================================
# Tests whether differentially expressed genes overlap with significant genomic
# cluster regions more than expected by chance, using a permutation test against
# the background of all expressed genes in the RNA-seq dataset.
#
# Author: Leah Darwin
#
# Output files:
#   (results printed to console; null distribution histogram displayed
#   interactively)
# =============================================================================

library(GenomicRanges)

set.seed(99)

##load gtf
gtf = read.csv(file="data/dmel-all-r6.57.gtf", sep="\t", header = FALSE) %>% 
  dplyr::filter(V3=="gene") %>%
  dplyr::select(V1, V4, V5, V9)%>%
  mutate(gene_id = str_extract(V9, "FBgn[0-9]+")) %>%
  mutate(gene_name = sub(".*gene_symbol ([^;]+);.*", "\\1", V9)) %>%
  dplyr::select(-V9)

colnames(gtf) = c("chr","start","end", "gene_id", "gene_name")

##load count data table 
cts = as.matrix(read.csv("data/GDL_count_table.csv", row.names = 1, check.names = FALSE)) 
gene_universe = as.data.frame(row.names(cts[rowSums(cts)>0, ])) 
colnames(gene_universe) = "gene_id"
rm(cts)
gene_universe = gene_universe %>% 
  left_join(gtf, by=join_by(gene_id)) %>%
  unique() %>%
  na.omit()


clust_bounds = read.csv("data/clust_bounds.csv")

selected_region_gr = GRanges(
  seqnames = clust_bounds$CHR,
  ranges = IRanges(start = clust_bounds$start,
                   end = clust_bounds$end)
)

# make GRanges objects
all_genes_gr <- GRanges(
  seqnames = gene_universe$chr,
  ranges = IRanges(start = gene_universe$start, 
                   end = gene_universe$end),
  gene_id = gene_universe$gene_id,
  gene_name = gene_universe$gene_name
)

observed_numerator = 16
observed_denominator = 42

# permutation test
n_perms <- 10000
null_dist <- replicate(n_perms, {
  random_genes <- all_genes_gr[sample(length(all_genes_gr), 
                                      observed_denominator)]
  sum(countOverlaps(random_genes, selected_region_gr) > 0)
})


# empirical p-value
p_value <- sum(null_dist >= observed_numerator) / n_perms
cat("Observed overlap:", observed_numerator, "\n")
cat("Expected overlap:", mean(null_dist), "\n")
cat("P-value:", p_value, "\n")

hist(null_dist, 
     main="Permutation null distribution",
     xlab="Number of overlapping genes",
     col="lightgrey",
     border="white")
abline(v=observed_numerator, col="red", lwd=2)

