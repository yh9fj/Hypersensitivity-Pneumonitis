################################################################################
# Figure 1D: DDRTree Dimensionality Reduction — Lung Biopsy Cohort (GSE150910)
#
# Description:
#   Applies DDRTree to voom-normalized RNA-seq data from the GSE150910 lung
#   biopsy cohort to visualize trajectory structure and separate Major/Minor
#   endotypes. This script mirrors the analysis in DDRTree.R (PFF/fHP),
#   DDRTree_UCD.R, and DDRTree_NJH.R.
#
# Input files:
#   - tpm.voom.csv                  : Voom-normalized TPM expression matrix
#                                      (genes x samples)
#   - meta.biopsie_endotype.csv     : Sample metadata including endotype labels
#                                      and rs35705950 genotype
#
# Output files:
#   - ddrtree.lung.biopsy.fhp.csv   : DDRTree 2D coordinates merged with
#                                      sample metadata
#   - Two base-R plots (rendered to active device or saveable via pdf())
#
# Dependencies:
#   install.packages("DDRTree")
#
# Reference:
#   Qi Mao et al. Dimensionality reduction via graph structure learning.
#   KDD 2015. https://doi.org/10.1145/2783258.2783309
################################################################################


# --- 0. Setup -----------------------------------------------------------------

library(DDRTree)

data_dir <- "../../Data/Lung.GSE150910"


# --- 1. Load Data -------------------------------------------------------------

# Expression matrix (genes x samples)
df1 <- read.csv(file.path(data_dir, "tpm.voom.csv"), row.names = 1)

# Sample metadata (endotype labels and rs35705950 genotype)
df2 <- read.csv(file.path(data_dir, "Biopsy/meta.biopsie_endotype.csv"), row.names = 1)


# --- 2. Align Samples ---------------------------------------------------------

idx      <- intersect(rownames(df2), colnames(df1))
expr_mtx <- as.matrix(df1[, idx])   # genes x samples expression matrix
meta     <- df2[idx, ]               # metadata aligned to same sample order


# --- 3. Define Sample Annotations ---------------------------------------------

Endotype   <- meta$Endotype      # Primary grouping variable (Major / Minor)
rs35705950 <- meta$rs35705950    # MUC5B promoter variant genotype


# --- 4. Run DDRTree -----------------------------------------------------------

set.seed(1234)

ddr_res <- DDRTree(
  expr_mtx,
  dimensions  = 2,
  maxIter     = 5,
  sigma       = 1e-2,
  lambda      = 1,
  ncenter     = 2,      # Number of centers = number of endotypes
  param.gamma = 10,
  tol         = 1e-2,
  verbose     = FALSE
)

# Extract DDRTree outputs
Z      <- ddr_res$Z      # Reduced-dimension sample coordinates (2 x n)
Y      <- ddr_res$Y      # Smooth principal tree node positions
stree  <- ddr_res$stree  # Spanning tree structure

Dim1 <- Z[1, ]
Dim2 <- Z[2, ]


# --- 5. Assign Plot Colors by Endotype ----------------------------------------

col_map <- c(Major = "blue", Minor = "brown")
pt_col  <- col_map[Endotype]


# --- 6. Plot A: DDRTree Reduced Dimensions ------------------------------------

par(mar = c(4, 4, 2, 1))

plot(
  Dim1, Dim2,
  col      = pt_col,
  pch      = 19,
  cex      = 1,
  cex.lab  = 1.5,
  cex.axis = 1.5,
  xlab     = "Dimension 1",
  ylab     = "Dimension 2",
  main     = "DDRTree Reduced Dimensions"
)

legend(
  "bottomright",
  legend  = names(col_map),
  col     = col_map,
  pch     = 19,
  pt.cex  = 1.5,
  cex     = 1.2,
  bty     = "n"
)


# --- 7. Plot B: Smooth Principal Tree Nodes -----------------------------------

plot(
  Y[1, ], Y[2, ],
  col  = "blue",
  pch  = 17,
  xlab = "Dimension 1",
  ylab = "Dimension 2",
  main = "DDRTree Smooth Principal Curves"
)


# --- 8. Export Results --------------------------------------------------------

colnames(Z) <- colnames(expr_mtx)
rownames(Z) <- c("PC1", "PC2")

coords_df <- as.data.frame(t(Z))

write.csv(
  cbind(coords_df, meta[rownames(coords_df), ]),
  file      = "ddrtree.lung.biopsy.fhp.csv",
  row.names = TRUE
)


# --- Session Info (for reproducibility) ---------------------------------------

sessionInfo()
