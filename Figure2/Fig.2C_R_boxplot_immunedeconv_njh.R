################################################################################
# Figure: Immune Cell Deconvolution Boxplots by Endotype — NJH Cohort
#
# Description:
#   Visualizes quanTIseq immune cell fractions across eight cell types
#   (Myeloid dendritic cell, Monocyte, Macrophage M2, B cell, NK cell,
#   CD8+ T cell, Tregs, CD4+ non-regulatory) stratified by endotype
#   (L-fHP vs N-fHP) in the NJH validation cohort. Group differences assessed
#   by Wilcoxon rank-sum test. Note: 8 cell types displayed vs. 5 in PFF/UCD
#   cohorts, reflecting additional myeloid populations captured in this dataset.
#   Mirrors the analysis in R_boxplot_immunedeconv_pff.R and
#   R_boxplot_immunedeconv_ucd.R.
#
# Input files:
#   - njh_fhp_decon_quantiseq3.csv  : quanTIseq deconvolution results for NJH
#                                      cohort; must contain 'SampleID' and
#                                      'Endotype' columns alongside cell-type
#                                      fractions
#
# Output files:
#   - boxplot_endotype_quantiseq_njh.png  : Raster figure (300 dpi, 20 x 6 in)
#   - boxplot_endotype_quantiseq_njh.pdf  : Vector figure (20 x 6 in)
#
# Dependencies:
#   install.packages(c("ggplot2", "tidyr", "dplyr", "ggpubr"))
################################################################################


# --- 0. Load Libraries --------------------------------------------------------

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)


# --- 1. Load and Prepare Data -------------------------------------------------

df <- read.csv("njh_fhp_decon_quantiseq3.csv", row.names = 1)

# Map numeric endotype codes to descriptive labels
df$Endotype_label <- factor(
  ifelse(df$Endotype == 0, "L-fHP", "N-fHP"),
  levels = c("L-fHP", "N-fHP")
)


# --- 2. Define Cell Types and Sample-Size Labels ------------------------------

# Cell type columns: exclude SampleID (col 1), Endotype (col 2), and col 11
# Adjust indices if column layout of input file differs
cell_types <- colnames(df)[-c(1, 10)]   # excludes Endotype and Endotype_label

n0 <- sum(df$Endotype == 0)
n1 <- sum(df$Endotype == 1)

x_labels <- c(
  "L-fHP" = sprintf("L-fHP\n(n=%d)", n0),
  "N-fHP" = sprintf("N-fHP\n(n=%d)", n1)
)


# --- 3. Reshape to Long Format ------------------------------------------------

cell_type_levels <- c(
  "Myeloid dendritic cell", "Monocyte", "Macrophage M2",
  "B cell", "NK cell", "CD8+ T cell", "Tregs", "CD4+ (non-regulatory)"
)

df_long <- df %>%
  pivot_longer(
    cols      = all_of(cell_types),
    names_to  = "CellType",
    values_to = "Fraction"
  ) %>%
  mutate(
    CellType = recode(CellType,
      "Myeloid_dendritic_cell" = "Myeloid dendritic cell",
      "Monocyte"               = "Monocyte",
      "Macrophage_M2"          = "Macrophage M2",
      "B.cell"                 = "B cell",
      "NK.cell"                = "NK cell",
      "CD8"                    = "CD8+ T cell",
      "Treg"                   = "Tregs",
      "CD4_non_regulatory"     = "CD4+ (non-regulatory)"
    ),
    CellType = factor(CellType, levels = cell_type_levels)
  )


# --- 4. Plot ------------------------------------------------------------------

colors <- c("L-fHP" = "#4878CF", "N-fHP" = "#D65F5F")

p <- ggplot(df_long,
            aes(x = Endotype_label, y = Fraction, fill = Endotype_label)) +

  geom_boxplot(
    width         = 0.5,
    outlier.shape = NA,
    alpha         = 0.7,
    linewidth     = 0.8
  ) +

  geom_jitter(width = 0.15, size = 1.5, alpha = 0.4) +

  # Significance stars
  stat_compare_means(
    method       = "wilcox.test",
    label        = "p.signif",
    comparisons  = list(c("L-fHP", "N-fHP")),
    size         = 6,
    bracket.size = 0.8
  ) +

  # Exact p-value
  stat_compare_means(
    method      = "wilcox.test",
    label       = "p.format",
    label.y.npc = 0.95,
    size        = 4
  ) +

  facet_wrap(~ CellType, nrow = 1, scales = "free_y") +

  scale_fill_manual(values = colors) +
  scale_x_discrete(labels = x_labels) +

  labs(
    title = "Immune Cell Deconvolution by Endotype (quanTIseq)",
    x     = NULL,
    y     = "Cell Fraction"
  ) +

  theme_bw(base_size = 14) +
  theme(
    plot.title         = element_text(size = 16, face = "bold", hjust = 0.5),
    strip.text         = element_text(size = 14, face = "bold"),
    strip.background   = element_rect(fill = "gray95"),
    axis.text.x        = element_text(size = 13),
    axis.text.y        = element_text(size = 12),
    axis.title.y       = element_text(size = 14),
    legend.position    = "bottom",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 13),
    legend.key.size    = unit(1.2, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )


# --- 5. Save Outputs ----------------------------------------------------------

# Width increased to 20 in to accommodate 8 facets
ggsave("boxplot_endotype_quantiseq_njh.png", p,
       width = 20, height = 6, dpi = 300)

ggsave("boxplot_endotype_quantiseq_njh.pdf", p,
       width = 20, height = 6)

cat("Saved: boxplot_endotype_quantiseq_njh.png and .pdf\n")


# --- Session Info (for reproducibility) ---------------------------------------

sessionInfo()
