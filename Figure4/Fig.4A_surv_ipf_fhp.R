################################################################################
# Figure 4A: Kaplan-Meier Survival — fHP/IPF Three-Group Comparison
#            L-fHP vs N-fHP vs IPF
#
# Description:
#   Plots unadjusted Kaplan-Meier transplant-free survival (TFS) curves for
#   three groups: L-fHP, N-fHP, and IPF, with a combined risk table.
#
# Input files:
#   - PFF.survival.IPF.fHP.workingFile.csv  : Survival data with columns:
#       ssid, Surv_time, TFS, Age_draw, Sex, Endotype,
#       FVCpct_base, DLCOpct_Base
#
# Output files:
#   - km_fHP_IPF.png / .pdf   : KM plot with risk table (7 x 5 in)
#
# Dependencies:
#   install.packages(c("survival", "survminer", "ggsurvfit", "ggfortify",
#                      "cowplot", "ggplot2", "dplyr", "tibble"))
################################################################################


# --- 0. Load Libraries --------------------------------------------------------

library(survival)
library(survminer)
library(ggsurvfit)
library(ggfortify)
library(cowplot)
library(ggplot2)
library(dplyr)
library(tibble)


# --- 1. Load and Prepare Data -------------------------------------------------

df <- read.csv("PFF.survival.IPF.fHP.workingFile.csv", row.names = 1)

data <- data.frame(
  time   = as.numeric(df$Surv_time),
  status = as.numeric(df$TFS),
  Dx     = df$Endotype,
  Age    = df$Age_draw,
  Sex    = df$Sex,
  FVCpp  = df$FVCpct_base,
  DLCOpp = df$DLCOpct_Base
)

labs <- levels(factor(data$Dx))


# --- 2. Kaplan-Meier Plot -----------------------------------------------------

km_fit <- survfit(Surv(time, status) ~ Dx, data = data)

km_plot <- ggsurvplot(
  km_fit,
  data             = data,
  fun              = "pct",
  linetype         = c("solid", "solid", "solid"),
  pval             = TRUE,
  conf.int         = FALSE,
  break.x.by      = 6,
  risk.table       = FALSE,
  fontsize         = 6,
  surv.median.line = "hv",
  ggtheme          = theme_light(),
  palette          = c("#8B7355", "#0072B2", "#D55E00"),
  legend.title     = "",
  legend.labs      = labs,
  font.main        = c(16, "bold", "black"),
  font.x           = c(14, "bold.italic", "black"),
  font.y           = c(14, "bold.italic", "black"),
  font.tickslab    = c(12, "plain", "black")
)

risk_table <- ggrisktable(
  km_fit,
  data          = data,
  break.time.by = 6,
  xlim          = c(0, 36)
)

km_combined <- plot_grid(
  km_plot$plot, risk_table,
  ncol        = 1,
  rel_heights = c(3, 1),
  align       = "v"
)

print(km_combined)


# --- 3. Save Outputs ----------------------------------------------------------

ggsave("km_fHP_IPF.png", km_combined, width = 7, height = 5, dpi = 300)
ggsave("km_fHP_IPF.pdf", km_combined, width = 7, height = 5)
cat("Saved: km_fHP_IPF.png and .pdf\n")


# --- Session Info (for reproducibility) ---------------------------------------

sessionInfo()
