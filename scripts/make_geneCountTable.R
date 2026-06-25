# =============================================================================
# Make RNA-seq Count Table and Sample Metadata
# =============================================================================
# Processes featureCounts output into a gene count matrix and derives sample
# group metadata from sample names for use in downstream RNA-seq analyses.
#
# Author: Yevgeniy Raynes
#
# Output files:
#   data/GDL_count_table.csv  : Gene count matrix (rows = gene IDs,
#                               columns = samples)
#   data/GDL_group_data.csv   : Sample metadata table with name, history,
#                               replicate, type, treatment, and group columns
# =============================================================================

library(tidyverse)

# ============================================================================
# BUILD COUNT TABLE
# ============================================================================

table <- data.frame(read.table("data/count_table_2pass_fC.txt",
                                sep = "\t", skip = 1, row.names = 1, header = TRUE))

table <- table[-c(1:5)] # remove featureCounts annotation columns
colnames(table) <- substring(colnames(table), 58, 66) # extract sample names from paths

write.csv(table, "data/GDL_count_table.csv", row.names = TRUE)

# ============================================================================
# BUILD SAMPLE METADATA
# ============================================================================

columnnames <- colnames(table)

groups <- data.frame(Sample = columnnames) %>%
  mutate(Name      = substring(Sample, 1, 3)) %>%
  mutate(History   = substring(Sample, 1, 1)) %>%
  mutate(Rep       = substring(Sample, 9, 9)) %>%
  mutate(Type      = substring(Sample, 5, 5)) %>%
  mutate(Treatment = substring(Sample, 7, 7)) %>%
  mutate(Group     = substring(Sample, 1, 7))

row.names(groups) <- groups$Sample
groups <- groups %>% select(-Sample)

write.csv(groups, "data/GDL_group_data.csv", row.names = TRUE)
