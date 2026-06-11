################################################################################
# Immune Deconvolution — PFF/fHP Training Cohort
#
# Description:
#   Processes raw RNA-seq count data from the PFF/fHP cohort through three
#   sequential steps:
#     1. Low-expression filtering of genes
#     2. TPM normalization using gene lengths (via scuttle)
#     3. Immune cell deconvolution using quanTIseq and MCP-counter
#        (via immunedeconv)
#
# Input files:
#   - counts.cmb3.csv               : Raw count matrix with gene annotation
#                                      columns (cols 1–20) followed by samples;
#                                      must contain 'Symbol' and 'Length' fields
#   - latent.clust.C2.allGenes.csv  : Sample metadata for subsetting to cohort
#
# Output files:
#   - pff_fhp_counts.flt.csv            : Filtered count matrix with gene lengths
#   - pff_fhp_tpm.csv                   : TPM-normalized expression matrix
#   - pff_fhp_tpm_decon_quantiseq.csv   : quanTIseq immune deconvolution results
#   - pff_fhp_tpm_decon_mcp_counter.csv : MCP-counter immune deconvolution results
#
# Dependencies:
#   BiocManager::install("scuttle")
#   install.packages("immunedeconv")
#   install.packages("CIBERSORT")   # required by immunedeconv for some methods
#
# Notes:
#   - TPM values are used as input (not log2-transformed) per immunedeconv
#     recommendations for quanTIseq and MCP-counter
#   - Gene filtering retains genes with > 10 counts in at least 10% of samples
################################################################################


# --- 0. Load Libraries --------------------------------------------------------

library(scuttle)
library(immunedeconv)
library(CIBERSORT)


# --- 1. Load Raw Count Data ---------------------------------------------------

df1 <- read.csv("counts.cmb3.csv", row.names = 1)

# Separate gene annotation columns (cols 1–19, after row names consumed) from
# count data; adjust index if annotation column count differs
annot <- df1[, 1:19]

# Load sample metadata
df2 <- read.csv("latent.clust.C2.allGenes.csv", row.names = 1)


# --- 2. Align Samples ---------------------------------------------------------

idx   <- intersect(rownames(df2), colnames(df1))
dd1   <- df1[, idx]
meta  <- df2[idx, ]


# --- 3. Filter Low-Expression Genes -------------------------------------------

# Retain genes with > 10 counts in at least 10% of samples
min_samples <- 0.1 * ncol(dd1)
keep        <- rowSums(dd1 > 10) > min_samples

counts_flt  <- dd1[keep, ]
annot_flt   <- annot[keep, ]
gene_length <- annot_flt$Length

write.csv(
  cbind(gLen = gene_length, counts_flt),
  file      = "pff_fhp_counts.flt.csv",
  row.names = TRUE
)


# --- 4. TPM Normalization -----------------------------------------------------

counts_mtx <- as.matrix(counts_flt)
tpm        <- calculateTPM(counts_mtx, lengths = gene_length)

write.csv(tpm, file = "pff_fhp_tpm.csv", row.names = TRUE)


# --- 5. Immune Deconvolution --------------------------------------------------

# quanTIseq: absolute immune cell fractions
# Input: raw TPM (not log2-transformed)
results_quantiseq <- deconvolute(tpm, method = "quantiseq")

write.csv(
  results_quantiseq,
  file      = "pff_fhp_tpm_decon_quantiseq.csv",
  row.names = FALSE
)

# MCP-counter: immune and stromal cell scores
# Input: raw TPM (not log2-transformed)
results_mcp <- deconvolute(tpm, method = "mcp_counter")

write.csv(
  results_mcp,
  file      = "pff_fhp_tpm_decon_mcp_counter.csv",
  row.names = FALSE
)


# --- Session Info (for reproducibility) ---------------------------------------

sessionInfo()
