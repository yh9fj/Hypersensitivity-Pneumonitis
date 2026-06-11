################################################################################
# Figure: Immune Cell Deconvolution Boxplots by Endotype — UCD Cohort
#
# Description:
#   Visualizes quanTIseq immune cell fractions across five cell types
#   (Neutrophil, B cell, NK cell, CD8+ T cell, Tregs) stratified by endotype
#   (L-fHP vs N-fHP) in the UCD validation cohort. Sample UCD2017106 is
#   excluded as an outlier. Group differences assessed by Wilcoxon rank-sum test.
#   Mirrors the analysis in R_boxplot_immunedeconv_pff.R (PFF/fHP cohort).
#
# Input files:
#   - ucd_fhp_decon_quantiseq3.csv  : quanTIseq deconvolution results for UCD
#                                      cohort; must contain 'SampleID' and
#                                      'Endotype' columns alongside cell-type
#                                      fractions
#
# Output files:
#   - boxplot_endotype_quantiseq_ucd.png  : Raster figure (300 dpi, 16 x 6 in)
#   - boxplot_endotype_quantiseq_ucd.pdf  : Vector figure (16 x 6 in)
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

df <- read.csv("ucd_fhp_decon_quantiseq3.csv", row.names = 1)

# Remove outlier sample
df <- df[rownames(df) != "UCD2017106", ]

# Map numeric endotype codes to descriptive labels
df$Endotype_label <- factor(
  ifelse(df$Endotype == 0, "L-fHP", "N-fHP"),
  levels = c("L-fHP", "N-fHP")
)


# --- 2. Define Cell Types and Sample-Size Labels ------------------------------

# Cell type columns: exclude SampleID (col 1), Endotype (col 2), and col 8
# Adjust indices if column layout of input file differs
cell_types <- colnames(df)[-c(1, 7)]   # excludes Endotype and Endotype_label

n0 <- sum(df$Endotype == 0)
n1 <- sum(df$Endotype == 1)

x_labels <- c(
  "L-fHP" = sprintf("L-fHP\n(n=%d)", n0),
  "N-fHP" = sprintf("N-fHP\n(n=%d)", n1)
)


# --- 3. Reshape to Long Format ------------------------------------------------

df_long <- df %>%
  pivot_longer(
    cols      = all_of(cell_types),
    names_to  = "CellType",
    values_to = "Fraction"
  ) %>%
  mutate(
    CellType = recode(CellType,
      "Neutrophil" = "Neutrophil",
      "B.cell"     = "B cell",
      "NK.cell"    = "NK cell",
      "CD8"        = "CD8+ T cell",
      "Treg"       = "Tregs"
    ),
    CellType = factor(CellType,
      levels = c("Neutrophil", "B cell", "NK cell", "CD8+ T cell", "Tregs")
    )
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

ggsave("boxplot_endotype_quantiseq_ucd.png", p,
       width = 16, height = 6, dpi = 300)

ggsave("boxplot_endotype_quantiseq_ucd.pdf", p,
       width = 16, height = 6)

cat("Saved: boxplot_endotype_quantiseq_ucd.png and .pdf\n")


# --- Session Info (for reproducibility) ---------------------------------------

sessionInfo()
