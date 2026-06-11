################################################################################
# Fig.1E: Unsupervised PCA of Biopsy Samples by Endotype
#
# Description:
#   Performs unsupervised PCA on the top 2,000 most variable genes from
#   voom-normalized RNA-seq data. Confidence ellipses are plotted per endotype
#   (Major / Minor). PC1 separation between endotypes is assessed by
#   Welch two-sample t-test.
#
# Input files:
#   - pff.tpm.voom.fhp.ipf.csv      : Voom-normalized expression matrix
#                                      (genes x samples)
#   - latend.coord.C2.allGenes.csv  : Sample metadata including endotype labels
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


data_dir <- "../../Data/PFF.fHP.ipf_RNAseq"


# --- 1. Load Data -------------------------------------------------------------

# Expression matrix (genes x samples)
df1 <- read.csv(file.path(data_dir, "pff.tpm.voom.fhp.ipf.csv"), row.names = 1)

# Sample metadata
df2 <- read.csv(file.path(data_dir, "latend.coord.C2.allGenes.csv"), row.names = 1)


# --- 2. Align Samples ---------------------------------------------------------

idx     <- intersect(rownames(df2), colnames(df1))
expr_mtx <- df1[, idx]
meta     <- df2[idx, ]

Endotype <- meta$Endotype


# --- 3. Select Top 2,000 Most Variable Genes ----------------------------------

gene_sd   <- apply(expr_mtx, 1, sd)
top_genes <- names(sort(gene_sd, decreasing = TRUE))[1:2000]

# Transpose to samples x genes for PCA input
pca_input <- as.data.frame(t(expr_mtx[top_genes, ]))


# --- 4. PCA -------------------------------------------------------------------

# Endotype added as qualitative supplementary variable (column 1)
pca_dat <- cbind(Endotype, pca_input)

res.pca <- PCA(pca_dat, quali.sup = 1, graph = FALSE)


# --- 5. Plot: PCA with Confidence Ellipses ------------------------------------

par(
  mar     = c(3.5, 3.5, 0.5, 0.5),
  mgp     = c(2, 0.7, 0),
  cex     = 1.25,
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