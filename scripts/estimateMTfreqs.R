library(tidyr)
library(nloptr)
library(ggplot2)
library(DirichletReg)
library(dplyr)

sync = read.csv("../data/joined_mt.sync", sep = "\t", header = FALSE)

preselect_sync = read.csv("../data/preselec_joined_mt.sync", sep="\t", header=FALSE)

labels = read.csv("../data/pca/sync_labs.csv")

rep = c("1","1A","1B","2","2A","2B")
time = c(2,16,16,2,16,16)

preselec_labels = data.frame(rep,time)

sra_id = read.csv("../data/sra_id.csv", header=FALSE)

snp = read.csv("../data/mitoSNP.table.temp", sep="\t")%>%
  select(contains("POS")|contains("REF")|contains("ALT")|any_of(sra_id$V2)|contains("Yak.Ore.A")|contains("Sm21.Ore.A"))%>%
  filter(nchar(ALT) == 1)

##add in Zim53
snpzim53 = read.csv("../data/zim53_var.table", sep="\t") %>%
  select(-c(CHROM, ID, ALT, REF))

colnames(snpzim53) = c("POS","Zim53")

snp = snp %>%
  left_join(snpzim53, by=join_by(POS)) %>%
  mutate(across(where(is.character), ~coalesce(.x, REF)))

##encoded SNP table
snp =  snp %>% mutate(across(-c(POS,REF,ALT), ~ case_when(
  . == REF ~ 0,   
  . == ALT ~ 1,   
  TRUE ~ NA_real_ 
))) %>% na.omit

p_ij = snp

rownames(p_ij) = p_ij$POS 

p_ij = p_ij %>%
  select(-c(POS,REF,ALT)) %>%
  na.omit()

p_ij = as.matrix(p_ij)

syncdf = sync %>%
  separate(V1, c(NA,NA,"pos"))%>%
  filter(pos %in% snp$POS) %>%
  mutate(ALT = snp$ALT) %>%
  mutate(REF = .data$V2) %>%
  select(-V2)

preselec_syncdf = preselect_sync %>%
  dplyr::rename(pos = "V2") %>%
  filter(pos %in% snp$POS) %>%
  mutate(ALT = snp$ALT) %>%
  mutate(REF = .data$V3) %>%
  select(-c(V1,V3)) 

extract_allele_counts = function(counts, ref, alt) {
  split_counts = strsplit(counts, ":", fixed = TRUE)  # Split by ":"
  counts_matrix = do.call(rbind, split_counts) %>% as.data.frame() %>% mutate_all(as.numeric)
  colnames(counts_matrix) = c("A", "T", "C", "G", "N", "del")
  
  # Extract REF and ALT counts
  ref_counts = counts_matrix[cbind(1:nrow(counts_matrix), match(ref, colnames(counts_matrix)))]
  alt_counts = counts_matrix[cbind(1:nrow(counts_matrix), match(alt, colnames(counts_matrix)))]
  
  # Compute ALT frequency and coverage
  coverage = ref_counts + alt_counts
  alt_freq = alt_counts / coverage
  
  return(list(alt_freq = alt_freq, coverage = coverage))
}

sample_cols = setdiff(names(syncdf), c("pos", "REF", "ALT"))  # Identify sample columns
preselec_sample_cols = setdiff(names(preselec_syncdf),c("pos", "REF", "ALT") )

f_j_df = syncdf %>% select(pos) 
preselec_f_j_df = preselec_syncdf %>% select(pos)


##fill in f_j for selec
for (sample in sample_cols) {
  result = extract_allele_counts(syncdf[[sample]], syncdf$REF, syncdf$ALT)
  f_j_df[[sample]] = result$alt_freq
}

##fill in f_j for preselec 
for(sample in preselec_sample_cols){
  result = extract_allele_counts(preselec_syncdf[[sample]], preselec_syncdf$REF, preselec_syncdf$ALT)
  preselec_f_j_df[[sample]] = result$alt_freq
}


# Objective function to minimize
loss_function = function(x, f_j) {
  expected_f = p_ij %*% x  # Compute expected frequencies
  sum((f_j - expected_f)^2) # Weighted sum of squared errors
}

# Equality constraint function: sum(x) - 1 = 0
sum_constraint <- function(x) {
  sum(x) - 1  # Ensure sum(x) = 1
}

get_opt_freqs = function(f_j){
  
  # Initial guess (equal proportions)
  x0 <- rep(1 / ncol(p_ij), ncol(p_ij))
  
  # Bounds (0 ≤ x ≤ 1)
  lower_bounds <- rep(0, ncol(p_ij))
  upper_bounds <- rep(1, ncol(p_ij))
  
  result <- nloptr(
    x0 = x0,
    eval_f = function(x) loss_function(x,  f_j),
    eval_g_eq = function(x) sum_constraint(x),
    lb = lower_bounds,
    ub = upper_bounds,
    opts = list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-6)
  )
  
  x = result$solution
  x = x / sum(x)   # enforce exact sum-to-1
  
  return(x)
  
 # return(result$solution)
}

mtfreqs = sapply(2:121, function(col_i) get_opt_freqs(f_j_df[,col_i]))
t_mtfreqs = t(mtfreqs)

preselec_mtfreqs = sapply(2:7, function(col_i) get_opt_freqs(preselec_f_j_df[,col_i]) )
t_preselec_mtfreqs = t(preselec_mtfreqs)

colnames(t_mtfreqs) = colnames(p_ij)
colnames(t_preselec_mtfreqs) = colnames(p_ij)

finaldf = cbind(labels,t_mtfreqs) 
preselec_finaldf = cbind(preselec_labels, t_preselec_mtfreqs)

write.csv(finaldf, "../data/mitofreqs.csv", quote=FALSE, row.names=FALSE)
write.csv(preselec_finaldf, "../data/mitofreqs_preselec.csv", quote=FALSE, row.names=FALSE)
