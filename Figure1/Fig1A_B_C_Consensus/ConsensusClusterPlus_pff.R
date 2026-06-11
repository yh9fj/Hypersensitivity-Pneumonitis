################################################################################
# Fig.1A-C: Consensus Clustering of Biopsy Samples
#
# Description:
#   Unsupervised consensus clustering of biopsy RNA-seq data (TPM/voom-normalized)
#   using the ConsensusClusterPlus package to identify robust patient subgroups.
#
# Input files:
#   - tpm.voom.csv         : Normalized expression matrix (genes x samples)
#   - meta.chp.biopsy2.csv : Sample metadata including latent cluster labels
#
# Output files:
#   - ConsensusClusterResults/ : PDF plots from ConsensusClusterPlus
#   - Biopsie.fHP_ICL/         : ICL diagnostic plots
#   - BiopsieCluster_Assignments_km.csv : Sample cluster assignments for optimal K
#
# Dependencies:
#   install.packages("BiocManager")
#   BiocManager::install("ConsensusClusterPlus")
#
# Reference:
#   Wilkerson MD, Hayes DN. ConsensusClusterPlus: a class discovery tool with
#   confidence assessments and item tracking. Bioinformatics. 2010;26(12):1572-3.
################################################################################


# --- 0. Load Libraries --------------------------------------------------------

library(ConsensusClusterPlus)


# --- 1. Load Data -------------------------------------------------------------

# Expression matrix (genes x samples)
df <- read.csv("tpm.voom.csv", row.names = 1)

# Sample metadata
df2 <- read.csv("meta.chp.biopsy2.csv", row.names = 1)

# Subset expression matrix to samples present in metadata
mtx <- as.matrix(df[, rownames(df2)])

# Retrieve latent cluster labels from metadata (used for downstream comparison)
latClust <- df2$latClust


# --- 2. Consensus Clustering --------------------------------------------------

set.seed(1234)

consensusResult <- ConsensusClusterPlus(
  d          = mtx,
  maxK       = 6,           # Maximum number of clusters to evaluate
  reps       = 100,         # Number of subsampling iterations
  pItem      = 0.8,         # Proportion of samples to subsample per iteration
  pFeature   = 1,           # No feature subsampling (clustering samples)
  clusterAlg = "km",        # k-means clustering
  distance   = "euclidean", # Euclidean distance
  seed       = 1234,
  plot       = "pdf",       # Save consensus plots to PDF
  title      = "ConsensusClusterResults"
)


# --- 3. ICL Diagnostics -------------------------------------------------------

# Item consensus and cluster consensus scores
resICL <- calcICL(consensusResult, title = "Biopsie.fHP")


# --- 4. Extract Cluster Assignments at Optimal K ------------------------------

optimalK <- 2   # Selected based on consensus CDF and delta area plots

clusters <- consensusResult[[optimalK]]$consensusClass

write.csv(
  data.frame(Sample = names(clusters), Cluster = clusters),
  file      = "BiopsieCluster_Assignments_km.csv",
  row.names = FALSE
)


# --- 5. Plot Consensus Dendrogram ---------------------------------------------

plot(
  consensusResult[[optimalK]]$consensusTree,
  main  = paste0("Consensus Dendrogram (K = ", optimalK, ")"),
  xlab  = "",
  sub   = "",
  cex   = 0.7
)


# --- Session Info (for reproducibility) ---------------------------------------

sessionInfo()