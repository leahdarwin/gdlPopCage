# Scripts

---

## Population Genetics

### `glm.sh`
SLURM batch script that runs a generalized linear model on genome-wide pooled sequencing data. Uses the `treatment_time_repl` model: calls `poolFreqDiff_treatment_time_repl.py` on `joined.sync` to generate an R input file, then runs it with `Rscript` to produce the final GLM output. This script primarily uses scripts that have been modified from poolFreqDiff (https://github.com/RAWWiberg/poolFreqDiff/) which are given in the directory tools/poolFreqDiff. 

**Outputs:**
- `${glm_dir}treatment_time_repl.glm` — GLM results file
- `${glm_dir}treatment_time_repl.glm.rin` — intermediate R input script

---

### `window_H.sh` and `window_H.py`
SLURM array job (6 tasks) that computes windowed heterozygosity across chromosomes for pre-selection pooled sequencing data. `window_H.sh` handles cluster setup and launches `window_H.py` once per array task. `window_H.py` reads a per-sample allele frequency file, computes per-site heterozygosity (2p(1−p)), then summarises mean H in 200 kb sliding windows (20 kb step, minimum 50 SNPs) across all chromosomes.

**Outputs:**
- `data/hetero/f{1..6}_H.tsv` — windowed mean heterozygosity per sample (chrom, window_start, window_end, n_sites, H_mean)

---

### `clusterSNPs.R`
Clusters significant SNPs on each chromosome using k-means, subsets haplotype frequency data to those cluster regions, and plots founder frequency trajectories and inter-cage Jaccard similarity heatmaps.

**Outputs:**
- `data/clust_bounds.csv` — Cluster genomic boundaries and founders that increased in frequency (>0.1 change F20→F60)
- `output/jaccard_[treatment].pdf` — Jaccard similarity heatmap between cage populations based on shared increased founders

---

### `plot_DNApca.R`
Performs PCA on genome-wide pooled sequencing allele frequencies across all samples. Plots PC1 vs PC3 by treatment with generation ellipses, PC2 vs generation by cage origin, and a scree plot.

**Outputs:**
- `output/pca_allelefreqs.pdf` — Combined figure with PC1 vs PC3 by treatment (left), PC2 vs generation by cage origin and scree plot (right)

---

### `plot_eigenSNPs.R`
Uses the afvaper package to compute eigendecomposition of allele frequency change vectors (rotenone vs control) per genomic window along chromosomes 2L, 2R, and 3R, with permutation-based significance thresholds.

**Outputs:**
- `output/eigenplots_final.pdf` — EV1 and EV2 proportion of variance plots for chromosomes 2L, 2R, and 3R

---

### `plot_man_hetero.R`
Creates the main combined figure: a genome-wide Manhattan plot of treatment-time interaction P-values stacked above per-cage heterozygosity tracks for chromosomes 2R and 3R.

**Outputs:**
- `output/man_hetero.pdf` — Manhattan plot (top) and heterozygosity tracks for chromosomes 2R and 3R (bottom)

---

### `plot_hetero_supp.R`
Plots windowed heterozygosity (H) for control vs rotenone samples at generation 60, per cage, for chromosomes 2L, 2R, and 3R. Significant SNP positions are marked with vertical lines.

**Outputs:**
- `output/hetero_2L.pdf` — Heterozygosity tracks for chromosome 2L, all cages
- `output/hetero_2R.pdf` — Heterozygosity tracks for chromosome 2R, all cages
- `output/hetero_3R.pdf` — Heterozygosity tracks for chromosome 3R, all cages

---

## Haplotype Frequencies

### `join_vcfFrq.sh`
SLURM batch script that converts a filtered VCF to a per-sample genotype table and then joins it with pooled allele frequencies on CHROM+POS. Runs in two steps: (1) uses `bcftools` to extract genotypes per sample, recoding homozygous REF as 1, homozygous ALT as 0, and heterozygous or missing calls as NA; (2) uses `awk` to left-join the resulting table with a `.frq` allele frequency file, appending frequency columns for matching sites. This file is required for haplotype calling. 

**Outputs:**
- `data/filtered_var.nohet.table` — genotype table (CHROM, POS, REF, ALT, per-sample 0/1/NA)
- `data/var_frq.nohet.tsv` — genotype table joined with pooled allele frequencies

### `plot_haploTrajs_CR.R`
Plots founder haplotype frequency trajectories across generations for each significant SNP cluster region, separately for control and rotenone treatments, using the increased founders identified in `clust_bounds_founders.csv`.

**Outputs:**
- `output/haplotraj_2L.pdf` — Frequency trajectories for the 2L cluster
- `output/haplotraj_2R.pdf` — Frequency trajectories for 2R clusters (stacked)
- `output/haplotraj_3R.pdf` — Frequency trajectories for the 3R cluster

---

### `plot_barhaplofreqs.R`
Plots stacked bar charts of mean founder haplotype frequencies across cages for specific genomic windows of interest, separately for control and rotenone treatments. Generation to plot is set by `gen_i` in the PARAMETERS section.

**Outputs:**
- `output/barfreqs_F[gen_i].pdf` — Stacked bar plots of founder frequencies at the specified generation for control (top) and rotenone (bottom) treatments

---

### `plot_contingTest.R`
Tests whether haplotype composition at an early generation predicts founder frequency at F60, using ILR-transformed PCA of haplotype frequencies and linear models with treatment interaction terms.

**Outputs:**
- `output/haploPCA_loadings_F[gen].pdf` — PC1 CLR loadings per genomic window
- `output/haploPCA_windows_F[gen].pdf` — PCA scores (PC1 vs PC2) per window, colored by treatment
- `output/haploPCA1_pred_F[gen].pdf` — PC1 score vs F60 founder frequency
- `output/haploPCA2_pred_F[gen].pdf` — PC2 score vs F60 founder frequency
- `output/haplofreq_pred_F[gen].pdf` — Direct founder frequency at F[gen] vs F60 with bootstrap error bars

---

## Mitochondrial DNA

### `estimateMTfreqs.R`
Estimates mitochondrial haplotype frequencies from pooled sequencing sync-format data using constrained nonlinear optimization (COBYLA), fitting observed SNP allele frequencies to a weighted sum of reference haplotypes.

**Outputs:**
- `data/mitofreqs.csv` — Estimated haplotype frequencies for all experimental evolution samples
- `data/mitofreqs_preselec.csv` — Estimated haplotype frequencies for pre-selection samples

---

### `plot_mtDNAdiffs.R`
Ordination (PCo) of mitochondrial haplotype frequencies using CLR-transformed composition distances, with PERMANOVA tests for treatment and generation effects. Also plots per-cage heatmaps of rotenone vs control haplotype proportion differences at each generation.

**Outputs:**
- `output/pco_mitos.pdf` — PCo1 and PCo2 vs generation, colored by treatment
- `output/R-C_mts.pdf` — Heatmaps of (Rotenone − Control) haplotype frequency differences per cage at F20, F40, F50, and F60

---

## Gene Mapping

### `plot_geneChrmap.R`
Maps the genomic positions of CYP450, Complex I, and UGT gene families onto a schematic chromosome diagram, overlaying the significant SNP cluster regions identified from the population-genetics analysis.

**Outputs:**
- `output/genemap.pdf` — Chromosome map with gene family positions and cluster regions highlighted

---

## RNA-seq

### `plot_RNAmds.R`
Performs multidimensional scaling (MDS) on VST-transformed RNA-seq count data using DESeq2, then plots samples by treatment and population type for all samples combined and for individual cage populations.

**Outputs:**
- `output/mdsplot_all.png` — MDS plot for all samples combined
- `output/mdsplot_1A3.png` — MDS plot for cage population 1A3
- `output/mdsplot_2A2.png` — MDS plot for cage population 2A2
- `output/mdsplot_2A3.png` — MDS plot for cage population 2A3
- `output/mdsplot_population_panel.png` — Combined panel of per-population MDS plots (1A3, 2A2, 2A3)

---

### `plot_volcEulerCounts_RNAseq.R`
Generates Euler (Venn-like) diagrams of DEG overlaps across cage populations, volcano plots highlighting genes in significant genomic cluster regions, and Jaccard overlap statistics. Also plots normalized read counts for selected genes.

**Outputs:**
- `output/eulers_all.pdf` — Euler diagrams of all-gene DEG overlaps across three contrasts
- `output/CC_pop_vols.pdf` — Volcano plots for control-population pairwise comparisons
- `output/RNAseq_fig_cyp.pdf` — Combined Euler diagrams, volcano plots, and normalized count plots for CYP genes
- `output/jaccard_overlaps.tex` — Jaccard similarity table for DEG set overlaps

---

### `rnaSeqOverlap_permTest.R`
Tests whether differentially expressed genes overlap with significant genomic cluster regions more than expected by chance, using a permutation test against the background of all expressed genes in the RNA-seq dataset.

**Outputs:** Results printed to console; null distribution histogram displayed interactively.

---

## Phenotype

### `plot_devtime.R`
Analyzes egg-to-adult development time (TSE) for rotenone and control treatments across two experimental generations (F48/F50 and F73/F75), fitting a linear mixed model and computing estimated marginal means and pairwise contrasts.

**Outputs:**
- `output/devtime.pdf` — Boxplot of development time by treatment
- `output/anova_devtime.tex` — ANOVA table for the linear mixed model
- `output/contrasts_devtime.tex` — Pairwise contrasts by cage treatment and origin
- `output/emm_devtime.tex` — Estimated marginal means table
- `output/contrasts_devtime2.tex` — Additional contrasts pooled across generations

---

### `plot_f1Cross.R`
Analyzes development time in reciprocal F1 crosses between control and rotenone cage populations, fitting a linear mixed model and computing pairwise contrasts by treatment.

**Outputs:**
- `output/F1_cross.pdf` — Boxplot of development time by cross type and treatment
- `output/anova_f1cross.tex` — ANOVA table for the F1 cross mixed model
- `output/contrasts_f1cross.tex` — Pairwise contrasts by treatment
