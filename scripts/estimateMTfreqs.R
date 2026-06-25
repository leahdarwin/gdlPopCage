# =============================================================================
# Mitochondrial Haplotype Frequency Estimation
# =============================================================================
# Estimates mitochondrial haplotype frequencies from pooled sequencing sync-
# format data using constrained nonlinear optimization (COBYLA), fitting
# observed SNP allele frequencies to a weighted sum of reference haplotypes.
#
# Author: Leah Darwin
#
# Output files:
#   data/mitofreqs.csv          : Estimated haplotype frequencies for all
#                                 experimental evolution samples
#   data/mitofreqs_preselec.csv : Estimated haplotype frequencies for
#                                 pre-selection samples
# =============================================================================

# ============================================================================
# LOAD PACKAGES
# ============================================================================

library(tidyr)         # For data reshaping
library(nloptr)        # For nonlinear optimization
library(ggplot2)       # For plotting (not actively used)
library(dplyr)         # For data manipulation

# ============================================================================
# DATA IMPORT
# ============================================================================

# Read pooled sequencing data (sync format) for mitochondrial genome
# Experimental evolution samples
sync = read.csv("data/joined_mt.sync", sep = "\t", header = FALSE)

# Pre-selection samples (before treatment started)
preselect_sync = read.csv("data/preselec_joined_mt.sync", sep="\t", header=FALSE)

# Sample labels and metadata
labels = read.csv("data/pca/sync_labs.csv")

# Pre-selection replicate labels and time points
rep = c("1","1A","1B","2","2A","2B")
time = c(2,16,16,2,16,16)
preselec_labels = data.frame(rep,time)

# SRA IDs for sequencing runs
sra_id = read.csv("data/sra_id.csv", header=FALSE)

# ============================================================================
# MITOCHONDRIAL SNP DATA
# ============================================================================

# Read SNP table with genotypes for reference strains
snp = read.csv("data/mitoSNP.table.temp", sep="\t")%>%
  select(contains("POS")|contains("REF")|contains("ALT")|any_of(sra_id$V2)|contains("Yak.Ore.A")|contains("Sm21.Ore.A"))%>%
  # Keep only single nucleotide variants (not indels)
  filter(nchar(ALT) == 1)

# Add Zim53 genotype data
snpzim53 = read.csv("data/zim53_var.table", sep="\t") %>%
  select(-c(CHROM, ID, ALT, REF))

colnames(snpzim53) = c("POS","Zim53")

# Merge Zim53 data, filling missing values with REF allele
snp = snp %>%
  left_join(snpzim53, by=join_by(POS)) %>%
  mutate(across(where(is.character), ~coalesce(.x, REF)))

# ============================================================================
# ENCODE SNP TABLE
# ============================================================================

# Convert SNP genotypes to binary encoding
# REF = 0, ALT = 1, missing = NA
snp =  snp %>% mutate(across(-c(POS,REF,ALT), ~ case_when(
  . == REF ~ 0,
  . == ALT ~ 1,
  TRUE ~ NA_real_
))) %>% na.omit

# ============================================================================
# CREATE p_ij MATRIX
# ============================================================================

# p_ij matrix: SNPs (rows) x Haplotypes (columns)
# Each entry is 0 (REF allele) or 1 (ALT allele) for that haplotype
p_ij = snp

rownames(p_ij) = p_ij$POS

p_ij = p_ij %>%
  select(-c(POS,REF,ALT)) %>%
  na.omit()

p_ij = as.matrix(p_ij)

# ============================================================================
# PREPARE SYNC DATA
# ============================================================================

# Extract and format experimental evolution sync data
syncdf = sync %>%
  separate(V1, c(NA,NA,"pos"))%>%
  # Filter to SNPs that are in the SNP table
  filter(pos %in% snp$POS) %>%
  mutate(ALT = snp$ALT) %>%
  mutate(REF = .data$V2) %>%
  select(-V2)

# Extract and format pre-selection sync data
preselec_syncdf = preselect_sync %>%
  dplyr::rename(pos = "V2") %>%
  filter(pos %in% snp$POS) %>%
  mutate(ALT = snp$ALT) %>%
  mutate(REF = .data$V3) %>%
  select(-c(V1,V3))

# ============================================================================
# FUNCTION: Extract Allele Counts from Sync Format
# ============================================================================

#' Extract REF and ALT allele frequencies from sync format
#'
#' @param counts Character vector of sync format counts (A:T:C:G:N:del)
#' @param ref Character vector of reference alleles
#' @param alt Character vector of alternate alleles
#' @return List with alt_freq (ALT frequency) and coverage (REF + ALT depth)
extract_allele_counts = function(counts, ref, alt) {
  # Split sync format string by ":"
  split_counts = strsplit(counts, ":", fixed = TRUE)
  counts_matrix = do.call(rbind, split_counts) %>% as.data.frame() %>% mutate_all(as.numeric)
  colnames(counts_matrix) = c("A", "T", "C", "G", "N", "del")

  # Extract REF and ALT counts by matching base identity
  ref_counts = counts_matrix[cbind(1:nrow(counts_matrix), match(ref, colnames(counts_matrix)))]
  alt_counts = counts_matrix[cbind(1:nrow(counts_matrix), match(alt, colnames(counts_matrix)))]

  # Compute ALT frequency and coverage
  coverage = ref_counts + alt_counts
  alt_freq = alt_counts / coverage

  return(list(alt_freq = alt_freq, coverage = coverage))
}

# ============================================================================
# CALCULATE OBSERVED ALT FREQUENCIES (f_j)
# ============================================================================

# Identify sample columns (exclude pos, REF, ALT)
sample_cols = setdiff(names(syncdf), c("pos", "REF", "ALT"))
preselec_sample_cols = setdiff(names(preselec_syncdf),c("pos", "REF", "ALT") )

# Initialize data frames to store ALT frequencies
f_j_df = syncdf %>% select(pos)
preselec_f_j_df = preselec_syncdf %>% select(pos)

# Calculate ALT frequencies for experimental evolution samples
for (sample in sample_cols) {
  result = extract_allele_counts(syncdf[[sample]], syncdf$REF, syncdf$ALT)
  f_j_df[[sample]] = result$alt_freq
}

# Calculate ALT frequencies for pre-selection samples
for(sample in preselec_sample_cols){
  result = extract_allele_counts(preselec_syncdf[[sample]], preselec_syncdf$REF, preselec_syncdf$ALT)
  preselec_f_j_df[[sample]] = result$alt_freq
}

# ============================================================================
# OPTIMIZATION FUNCTIONS
# ============================================================================

#' Loss function to minimize
#'
#' Calculates sum of squared errors between observed and expected frequencies
#'
#' @param x Vector of haplotype frequencies (to be optimized)
#' @param f_j Vector of observed ALT frequencies
#' @return Sum of squared errors
loss_function = function(x, f_j) {
  # Expected frequency = weighted sum of haplotype frequencies
  expected_f = p_ij %*% x
  sum((f_j - expected_f)^2) # Sum of squared errors
}

#' Equality constraint function
#'
#' Ensures that haplotype frequencies sum to 1
#'
#' @param x Vector of haplotype frequencies
#' @return Deviation from sum = 1 (should equal 0)
sum_constraint <- function(x) {
  sum(x) - 1
}

#' Optimize haplotype frequencies for a single sample
#'
#' @param f_j Vector of observed ALT frequencies for one sample
#' @return Vector of optimized haplotype frequencies
get_opt_freqs = function(f_j){

  # Initial guess: equal proportions for all haplotypes
  x0 <- rep(1 / ncol(p_ij), ncol(p_ij))

  # Bounds: each frequency must be between 0 and 1
  lower_bounds <- rep(0, ncol(p_ij))
  upper_bounds <- rep(1, ncol(p_ij))

  # Run constrained optimization
  # Algorithm: COBYLA (Constrained Optimization BY Linear Approximations)
  result <- nloptr(
    x0 = x0,
    eval_f = function(x) loss_function(x,  f_j),
    eval_g_eq = function(x) sum_constraint(x),
    lb = lower_bounds,
    ub = upper_bounds,
    opts = list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-6)
  )

  # Enforce exact sum-to-1 constraint
  x = result$solution
  x = x / sum(x)

  return(x)
}

# ============================================================================
# ESTIMATE HAPLOTYPE FREQUENCIES FOR ALL SAMPLES
# ============================================================================

# Optimize for experimental evolution samples (columns 2-121)
mtfreqs = sapply(2:121, function(col_i) get_opt_freqs(f_j_df[,col_i]))
t_mtfreqs = t(mtfreqs)

# Optimize for pre-selection samples (columns 2-7)
preselec_mtfreqs = sapply(2:7, function(col_i) get_opt_freqs(preselec_f_j_df[,col_i]) )
t_preselec_mtfreqs = t(preselec_mtfreqs)

# Add haplotype names as column headers
colnames(t_mtfreqs) = colnames(p_ij)
colnames(t_preselec_mtfreqs) = colnames(p_ij)

# ============================================================================
# COMBINE WITH METADATA AND SAVE
# ============================================================================

# Combine estimated frequencies with sample labels
finaldf = cbind(labels,t_mtfreqs)
preselec_finaldf = cbind(preselec_labels, t_preselec_mtfreqs)

# Save results to CSV files
write.csv(finaldf, "data/mitofreqs.csv", quote=FALSE, row.names=FALSE)
write.csv(preselec_finaldf, "data/mitofreqs_preselec.csv", quote=FALSE, row.names=FALSE)
