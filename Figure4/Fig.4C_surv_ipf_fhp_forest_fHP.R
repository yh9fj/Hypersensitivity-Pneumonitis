################################################################################
# Figure 4C: Kaplan-Meier Survival and Multivariable Cox PH Forest Plot
#            Endotype-Adjusted Analysis — fHP Cohort
#
# Description:
#   1. Kaplan-Meier curves for transplant-free survival (TFS) by endotype
#      (L-fHP vs N-fHP), unadjusted
#   2. Univariable Cox PH: TFS ~ Endotype
#   3. Multivariable Cox PH: TFS ~ Endotype + Age + Sex + FVC% + DLCO% + Steroid
#   4. Custom ggplot2 forest plot of multivariable Cox HR estimates
#   5. survminer::ggforest() shaded forest plot (alternative)
#
# Input files:
#   - PFF.survival.fHP.workingFile.csv  : Survival data with columns:
#       ssid, Surv_time, TFS, Age_draw, Sex, FVCpct_base, DLCOpct_Base,
#       Endotype, Steroid_before_bloodDraw
#
# Output files:
#   - km_fHP_endotype.png / .pdf           : KM plot with risk table (7 x 5 in)
#   - forest_fHP_adj.png / .pdf            : Custom forest plot (7 x 5 in)
#   - forest_fHP_adj_shade.png / .pdf      : ggforest shaded plot (7 x 5 in)
#
# Dependencies:
#   install.packages(c("survival", "survminer", "ggsurvfit", "ggfortify",
#                      "cowplot", "forestplot", "ggplot2", "dplyr", "tibble"))
#
# Reported results (for verification):
#   Univariable Cox: N-fHP HR 3.36 (1.77–6.38), p=0.0002, Concordance=0.639
#   Multivariable Cox (n=116, events=37):
#     N-fHP HR 2.16 (1.10–4.22), p=0.025; FVC% HR 0.955, p=0.003
################################################################################


# --- 0. Load Libraries --------------------------------------------------------

library(survival)
library(survminer)
library(ggsurvfit)
library(ggfortify)
library(cowplot)
library(forestplot)
library(ggplot2)
library(dplyr)
library(tibble)


# --- 1. Load and Prepare Data -------------------------------------------------

df <- read.csv("PFF.survival.fHP.workingFile.csv", row.names = 1)

Time    <- as.numeric(df$Surv_time)
TFS     <- as.numeric(df$TFS)
Age     <- as.numeric(df$Age_draw)
Sex     <- ifelse(tolower(df$Sex) == "male", 1, 0)
FVCpp   <- df$FVCpct_base
DLCOpp  <- df$DLCOpct_Base
Endotype <- df$Endotype
Steroid  <- as.integer(as.logical(df$Steroid_before_bloodDraw))

labs <- levels(factor(Endotype))


# --- 2. Assemble Analysis Data Frame ------------------------------------------

data <- data.frame(
  time    = Time,
  status  = TFS,
  Dx      = Endotype,
  Age     = Age,
  Sex     = Sex,
  FVCpp   = FVCpp,
  DLCOpp  = DLCOpp,
  Steroid = Steroid
)


# --- 3. Kaplan-Meier Plot -----------------------------------------------------

km_fit <- survfit(Surv(time, status) ~ Dx, data = data)

km_plot <- ggsurvplot(
  km_fit,
  data             = data,
  fun              = "pct",
  linetype         = c("solid", "solid"),
  pval             = TRUE,
  conf.int         = FALSE,
  break.x.by      = 6,
  risk.table       = FALSE,
  fontsize         = 6,
  surv.median.line = "hv",
  ggtheme          = theme_light(),
  palette          = c("#0072B2", "#D55E00"),
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

ggsave("km_fHP_endotype.png", km_combined, width = 7, height = 5, dpi = 300)
ggsave("km_fHP_endotype.pdf", km_combined, width = 7, height = 5)
cat("Saved: km_fHP_endotype.png and .pdf\n")


# --- 4. Univariable Cox PH: TFS ~ Endotype ------------------------------------

coxfit_uni <- coxph(Surv(time, status) ~ Dx, data = data, ties = "exact")
cat("\n--- Univariable Cox PH ---\n")
print(summary(coxfit_uni))


# --- 5. Multivariable Cox PH --------------------------------------------------

coxfit_multi <- coxph(
  Surv(time, status) ~ Dx + Age + Sex + FVCpp + DLCOpp + Steroid,
  data  = data,
  ties  = "exact"
)
cat("\n--- Multivariable Cox PH ---\n")
print(summary(coxfit_multi))


# --- 6. Custom Forest Plot ----------------------------------------------------

summ <- summary(coxfit_multi)

forest_data <- data.frame(
  variable = c("N-fHP vs L-fHP", "Age (years)", "Sex (M vs F)",
               "FVC (% predicted)", "DLCO (% predicted)", "Steroid (Yes vs No)"),
  HR       = summ$conf.int[, "exp(coef)"],
  CI_lower = summ$conf.int[, "lower .95"],
  CI_upper = summ$conf.int[, "upper .95"],
  p_value  = summ$coefficients[, "Pr(>|z|)"]
)

# Reverse factor order so first variable plots at top
forest_data$variable <- factor(forest_data$variable,
                                levels = rev(forest_data$variable))

p_forest <- ggplot(forest_data, aes(y = variable, x = HR)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red",
             linewidth = 0.8) +
  geom_point(size = 4, shape = 18) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper),
                 height = 0.3, linewidth = 0.8) +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10)) +
  labs(
    x     = "Hazard Ratio (95% CI)",
    y     = NULL,
    title = "Multivariable Cox Regression"
  ) +
  theme_bw() +
  theme(
    axis.text.y      = element_text(size = 12, color = "black"),
    axis.text.x      = element_text(size = 11),
    axis.title.x     = element_text(size = 13, face = "bold"),
    plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  ) +
  annotate(
    "text",
    x     = max(forest_data$CI_upper) * 1.5,
    y     = seq_len(nrow(forest_data)),
    label = paste0(
      "HR: ", sprintf("%.2f", forest_data$HR),
      "\np=", format.pval(forest_data$p_value, digits = 2)
    ),
    size = 3.5
  )

print(p_forest)

ggsave("forest_fHP_adj.png", p_forest, width = 7, height = 5, dpi = 300)
ggsave("forest_fHP_adj.pdf", p_forest, width = 7, height = 5)
cat("Saved: forest_fHP_adj.png and .pdf\n")


# --- 7. ggforest Shaded Forest Plot (Alternative) -----------------------------

p_ggforest <- ggforest(coxfit_multi, data = data)
print(p_ggforest)

ggsave("forest_fHP_adj_shade.png", p_ggforest, width = 7, height = 5, dpi = 300)
ggsave("forest_fHP_adj_shade.pdf", p_ggforest, width = 7, height = 5)
cat("Saved: forest_fHP_adj_shade.png and .pdf\n")


# --- Session Info (for reproducibility) ---------------------------------------

sessionInfo()
