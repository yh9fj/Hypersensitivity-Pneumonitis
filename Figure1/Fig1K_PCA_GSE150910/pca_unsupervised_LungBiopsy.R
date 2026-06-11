################################################################################
# Figure 1H: Unsupervised PCA — Lung Biopsy Cohort (GSE150910)
#
# Description:
#   Performs unsupervised PCA on the top 1,000 most variable genes from
#   voom-normalized RNA-seq data (GSE150910 lung biopsy cohort). Confidence
#   ellipses are plotted per endotype (Major / Minor). PC1 separation between
#   endotypes is assessed by Welch two-sample t-test (reported p = 0.01225).
#   Mirrors the analysis in pca_unsupervised.R (PFF/fHP), pca_unsupervised_UCD.R,
#   and pca_unsupervised_NJH.R. Note: top 1,000 genes used here vs. 2,000 in
#   other cohorts, reflecting the smaller sample size of this dataset.
#
# Input files:
#   - tpm.voom.csv                  : Voom-normalized expression matrix
#                                      (genes x samples)
#   - meta.biopsie_endotype.csv     : Sample metadata including endotype labels
#
# Output:
#   - PCA biplot with 95% confidence ellipses (rendered to active device)
#   - Welch t-test result for PC1 ~ Endotype printed to console
#
# Dependencies:
#   install.packages(c("FactoMineR", "dplyr"))
################################################################################


# --- 0. Load Libraries --------------------------------------------------------

library(FactoMineR)
library(dplyr)

data_dir <- "../../Data/Lung.GSE150910"


# --- 1. Load Data -------------------------------------------------------------

# Expression matrix (genes x samples)
df1 <- read.csv(file.path(data_dir, "tpm.voom.csv"), row.names = 1)

# Sample metadata (endotype labels)
df2 <- read.csv(file.path(data_dir, "Biopsy/meta.biopsie_endotype.csv"), row.names = 1)


# --- 2. Align Samples ---------------------------------------------------------

idx      <- intersect(rownames(df2), colnames(df1))
expr_mtx <- df1[, idx]
meta     <- df2[idx, ]

Endotype <- meta$Endotype     # Note: source column may be "Endotyp" — verify if needed


# --- 3. Select Top 1,000 Most Variable Genes ----------------------------------
# Note: 1,000 genes used (vs. 2,000 in other cohorts) due to smaller sample size

gene_sd   <- apply(expr_mtx, 1, sd)
top_genes <- names(sort(gene_sd, decreasing = TRUE))[1:1000]

pca_input <- as.data.frame(t(expr_mtx[top_genes, ]))


# --- 4. PCA -------------------------------------------------------------------

# Endotype added as qualitative supplementary variable (column 1)
pca_dat <- cbind(Endotype, pca_input)

res.pca <- PCA(pca_dat, quali.sup = 1, graph = FALSE)


# --- 5. Plot: PCA with Confidence Ellipses ------------------------------------

par(
  mar      = c(3.5, 3.5, 0.5, 0.5),
  mgp      = c(2, 0.7, 0),
  cex      = 1.25,
  cex.main = 1
)

plotellipses(
  res.pca,
  keepvar    = 1,
  label      = "none",
  col.hab    = c(Major = "blue", Minor = "brown"),
  title      = "",
  graph.type = "classic"
)


# --- 6. Welch t-test: PC1 ~ Endotype ------------------------------------------

pc1_coords <- res.pca$ind$coord[, 1]

cat("\n--- Welch Two-Sample t-test: PC1 by Endotype ---\n")
print(t.test(pc1_coords ~ Endotype))


# --- Session Info (for reproducibility) ---------------------------------------

sessionInfo()
