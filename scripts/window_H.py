# =============================================================================
# Windowed Heterozygosity from Pooled Sequencing Data
# =============================================================================
# Computes mean per-site heterozygosity (2p(1-p)) in sliding windows across
# chromosomes for a single pooled sequencing sample. Designed to be run as a
# SLURM array job (window_H.sh), with the sample index passed as a command-line
# argument.
#
# Author: Leah Darwin, Camille Brown
#
# Input files:
#   joined.sync.MAF01.frq : Per-site allele frequencies (MAF-filtered),
#                           one column per sample (f1..f6)
#
# Output files:
#   data/hetero/f{pop_idx}_H.tsv : Windowed mean heterozygosity
#                             (chrom, window_start, window_end, n_sites, H_mean)
# =============================================================================

import pandas as pd
import numpy as np
import sys

# -----------------------------
# PARAMETERS
# -----------------------------
WINDOW_SIZE = 200000   # bp
STEP_SIZE   = 20000    # bp
MIN_SITES   = 50      # minimum SNPs per window

file = "data/joined.sync.MAF01.frq"

df = pd.read_csv(file, sep="\t")

pop_idx = int(sys.argv[1])
p_col = f"f{pop_idx}"

# Rename selected f# column to p
df = df.rename(columns={p_col: "p"})

# enforce finite values
df = df[(df["p"] >= 0) & (df["p"] <= 1)]

# Compute per-site heterozygosity
df["H"] = 2 * df["p"] * (1 - df["p"])

# -----------------------------
# SLIDING WINDOW FUNCTION
# -----------------------------
def sliding_windows(df_chr, window, step):
    results = []

    start = df_chr["pos"].min()
    end   = df_chr["pos"].max()

    for w_start in range(start, end - window + 1, step):
        w_end = w_start + window

        subset = df_chr[
            (df_chr["pos"] >= w_start) &
            (df_chr["pos"] <  w_end)
        ]

        if len(subset) < MIN_SITES:
            continue

        # Mean heterozygosity
        H_mean = subset["H"].mean()

        results.append({
            "chrom": df_chr["chrom"].iloc[0],
            "window_start": w_start,
            "window_end": w_end,
            "n_sites": len(subset),
            "H_mean": H_mean,
        })

    return results


# -----------------------------
# RUN PER CHROMOSOME
# -----------------------------
all_results = []

for chrom, df_chr in df.groupby("chrom"):
    all_results.extend(sliding_windows(df_chr, WINDOW_SIZE, STEP_SIZE))

out = pd.DataFrame(all_results)

# -----------------------------
# SAVE OUTPUT
# -----------------------------
out.to_csv(f"data/hetero/f{pop_idx}_H.tsv",
           sep="\t", index=False)




