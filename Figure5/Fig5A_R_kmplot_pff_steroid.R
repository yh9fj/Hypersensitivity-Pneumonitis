################################################################################
# Figure 3C: Kaplan-Meier Survival and Cox Proportional Hazards Analysis
#            Endotype × Steroid Use — PFF/fHP Cohort
#
# Description:
#   Examines transplant-free survival stratified by endotype (L-fHP / N-fHP)
#   and steroid use before blood draw (steroid+ / steroid-), forming four
#   groups: LS-, LS+, NS-, NS+. Analysis includes:
#     1. Kaplan-Meier curves with risk table
#     2. Cox proportional hazards model (Group as predictor)
#     3. Pairwise log-rank tests
#     4. Pairwise Cox-PH contrasts with HR, 95% CI, and BH-adjusted p-values
#     5. DLCO t-tests by steroid use within each endotype
#
# Input files:
#   - pff.steroid.endotype.surv3.csv  : Survival data with columns including
#       endotype, Steroid_before_bloodDraw, Surv_time, VS, DLCOpct_Base
#
# Output files:
#   - km_steroid_pff.png / .pdf       : KM plot with risk table (7 x 7 in)
#   - pairwise_cox_results.csv        : Pairwise Cox HR table with adjusted p
#
# Dependencies:
#   install.packages(c("survival", "survminer", "ggsurvfit", "ggfortify",
#                      "cowplot", "ggplot2", "dplyr", "tibble", "knitr"))
#
# Reported results (for verification):
#   Cox PH (n=119, events=38): Concordance=0.69
#   Pairwise log-rank: LS- vs NS- p=1.3e-05; LS- vs LS+ p=0.007
#   DLCO t-test (all): steroid- 45.9% vs steroid+ 39.1%, p=0.019
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
library(knitr)


# --- 1. Load and Prepare Data -------------------------------------------------

df <- read.csv("pff.steroid.endotype.surv3.csv", row.names = 1)

# Combine endotypes 0/1 and normalize steroid flag to 0/1
df$dlco  <- as.numeric(df$DLCOpct_Base)
df$st.bf <- as.integer(as.logical(df$Steroid_before_bloodDraw))
df$time  <- as.numeric(df$Surv_time)

# Recode survival status: alive=0, dead/transplant=1
df$status <- ifelse(df$VS == "alive", 0, 1)


# --- 2. Define Four-Group Variable --------------------------------------------

group_code <- paste(df$endotype, df$st.bf, sep = "_")
df$Group <- factor(
  dplyr::recode(group_code,
    "0_0" = "LS-", "0_1" = "LS+",
    "1_0" = "NS-", "1_1" = "NS+"
  ),
  levels = c("LS-", "LS+", "NS-", "NS+")
)

labs <- levels(df$Group)

# Endotype subsets (used for DLCO t-tests)
dd0 <- subset(df, endotype == 0)
dd1 <- subset(df, endotype == 1)


# --- 3. Kaplan-Meier Plot -----------------------------------------------------

km_fit <- survfit(Surv(time, status) ~ Group, data = df)

km_plot <- ggsurvplot(
  km_fit,
  data            = df,
  fun             = "pct",
  linetype        = c("solid", "solid", "solid", "solid"),
  size            = c(1, 4, 1, 4),      # LS+/NS+ drawn thicker
  pval            = TRUE,
  conf.int        = FALSE,
  break.x.by     = 6,
  risk.table      = FALSE,
  fontsize        = 4,
  surv.median.line = "hv",
  ggtheme         = theme_light(),
  palette         = c("#0072B2", "#0072B2", "#D55E00", "#D55E00"),
  legend.title    = "",
  legend.labs     = labs,
  font.main       = c(10, "bold", "black"),
  font.x          = c(10, "bold.italic", "black"),
  font.y          = c(10, "bold.italic", "black"),
  font.tickslab   = c(10, "plain", "black")
)

risk_table <- ggrisktable(
  km_fit,
  data           = df,
  break.time.by  = 6,
  xlim           = c(0, 36)
)

km_combined <- plot_grid(
  km_plot$plot, risk_table,
  ncol        = 1,
  rel_heights = c(3, 1),
  align       = "v"
)

ggsave("km_steroid_pff.png", km_combined, width = 7, height = 7, dpi = 300)
ggsave("km_steroid_pff.pdf", km_combined, width = 7, height = 7)
cat("Saved: km_steroid_pff.png and .pdf\n")


# --- 4. Cox Proportional Hazards Model ----------------------------------------

coxfit <- coxph(Surv(time, status) ~ Group, data = df, ties = "exact")
print(summary(coxfit))


# --- 5. Pairwise Log-Rank Tests -----------------------------------------------

cat("\n--- Pairwise Log-Rank Tests (unadjusted) ---\n")
print(pairwise_survdiff(
  Surv(time, status) ~ Group,
  data            = df,
  p.adjust.method = "none"
))


# --- 6. Pairwise Cox-PH Contrasts ---------------------------------------------

coefs        <- coef(coxfit)
vcov_matrix  <- vcov(coxfit)
group_levels <- levels(df$Group)
n_levels     <- length(group_levels)

pairwise_results <- data.frame()

# Reference (LS-) vs all others — directly from model summary
summ <- summary(coxfit)
for (i in seq_along(coefs)) {
  pairwise_results <- rbind(pairwise_results, data.frame(
    Group1   = group_levels[1],
    Group2   = gsub("Group", "", names(coefs)[i]),
    HR       = summ$conf.int[i, "exp(coef)"],
    CI_lower = summ$conf.int[i, "lower .95"],
    CI_upper = summ$conf.int[i, "upper .95"],
    p_value  = summ$coefficients[i, "Pr(>|z|)"]
  ))
}

# All other pairwise contrasts via linear combination
for (i in 2:(n_levels - 1)) {
  for (j in (i + 1):n_levels) {
    contrast    <- rep(0, length(coefs))
    contrast[i - 1] <-  1
    contrast[j - 1] <- -1

    log_hr   <- sum(contrast * coefs)
    se       <- as.numeric(sqrt(t(contrast) %*% vcov_matrix %*% contrast))
    hr       <- exp(log_hr)
    ci_lower <- exp(log_hr - 1.96 * se)
    ci_upper <- exp(log_hr + 1.96 * se)
    z_score  <- log_hr / se
    p_value  <- 2 * pnorm(abs(z_score), lower.tail = FALSE)

    pairwise_results <- rbind(pairwise_results, data.frame(
      Group1   = group_levels[i],
      Group2   = group_levels[j],
      HR       = hr,
      CI_lower = ci_lower,
      CI_upper = ci_upper,
      p_value  = p_value
    ))
  }
}

# BH-adjusted p-values
pairwise_results$p_adjusted <- p.adjust(pairwise_results$p_value, method = "BH")

cat("\n--- Pairwise Cox-PH Contrasts ---\n")
print(pairwise_results)

write.csv(pairwise_results, file = "pairwise_cox_results.csv", row.names = FALSE)
cat("Saved: pairwise_cox_results.csv\n")


# --- 7. Flipped HR: NS+ vs NS- ------------------------------------------------
# The model computes NS- vs NS+; invert to report NS+ vs NS- direction

ns_row       <- pairwise_results[pairwise_results$Group1 == "NS-" &
                                   pairwise_results$Group2 == "NS+", ]
hr_flipped   <- 1 / ns_row$HR
ci_l_flipped <- 1 / ns_row$CI_upper   # bounds swap on inversion
ci_u_flipped <- 1 / ns_row$CI_lower

cat(sprintf(
  "\nNS+ vs NS-: HR = %.3f, 95%% CI [%.3f, %.3f]\n",
  hr_flipped, ci_l_flipped, ci_u_flipped
))


# --- 8. DLCO t-tests by Steroid Use -------------------------------------------

cat("\n--- DLCO by Steroid Use ---\n")

cat("\nL-fHP (endotype 0):\n")
print(t.test(dd0$dlco ~ dd0$st.bf))

cat("\nN-fHP (endotype 1):\n")
print(t.test(dd1$dlco ~ dd1$st.bf))

cat("\nAll samples:\n")
print(t.test(df$dlco ~ df$st.bf))


# --- Session Info (for reproducibility) ---------------------------------------

sessionInfo()
