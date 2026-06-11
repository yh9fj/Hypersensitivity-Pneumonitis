################################################################################
# Figure 5C: Cox PH — Endotype × Steroid Interaction Plot
#            Adjusted for Age, Sex, FVC%, DLCO%
#
# Description:
#   Fits a multivariable Cox proportional hazards model with a Steroid ×
#   Endotype interaction term to assess whether steroid benefit on
#   transplant-free survival (TFS) differs by endotype (L-fHP vs N-fHP).
#   Steroid HRs and 95% CIs are computed separately for each endotype via
#   linear contrasts, tested with glht (multcomp), and displayed as an
#   interaction plot.
#
# Input files:
#   - ../pff.steroid.endotype.surv3.csv  : Survival data with columns:
#       endotype, Steroid_before_bloodDraw, Surv_time, VS,
#       Age_draw, Sex, FVCpct_base, DLCOpct_Base
#
# Output files:
#   - interaction_plot_steroid_endotype.png / .pdf  : Interaction HR plot
#   - steroid_effect_by_endotype.csv               : HR table with p-values
#
# Dependencies:
#   install.packages(c("survival", "survminer", "ggsurvfit", "ggfortify",
#                      "cowplot", "ggplot2", "dplyr", "tibble", "multcomp"))
#
# Reported results (for verification):
#   Interaction p = 0.048
#   L-fHP steroid HR: 2.83 (95% CI 1.07–7.49), p = 0.036
#   N-fHP steroid HR: 0.72 (95% CI 0.28–1.85), p = 0.491
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
library(multcomp)


# --- 1. Load and Prepare Data -------------------------------------------------

df <- read.csv("../pff.steroid.endotype.surv3.csv", row.names = 1)

df$dlco   <- as.numeric(df$DLCOpct_Base)
df$st.bf  <- as.integer(as.logical(df$Steroid_before_bloodDraw))
df$time   <- as.numeric(df$Surv_time)
df$status <- ifelse(df$VS == "alive", 0, 1)

# Numeric sex coding: female = 0, male = 1
df$Sex <- ifelse(tolower(df$Sex) == "male", 1, 0)

df$FVCpp  <- df$FVCpct_base
df$DLCOpp <- df$DLCOpct_Base
df$Age    <- df$Age_draw
df$Steroid  <- df$st.bf


# --- 2. Define Four-Group Variable (for reference / downstream use) -----------

group_code <- paste(df$endotype, df$st.bf, sep = "_")
df$Group <- factor(
  dplyr::recode(group_code,
    "0_0" = "LS-", "0_1" = "LS+",
    "1_0" = "NS-", "1_1" = "NS+"
  ),
  levels = c("LS-", "LS+", "NS-", "NS+")
)


# --- 3. Multivariable Cox PH with Interaction ---------------------------------

coxfit <- coxph(
  Surv(time, status) ~ Steroid * Endotype + Age + Sex + FVCpp + DLCOpp,
  data = df,
  ties = "exact"
)

print(summary(coxfit))


# --- 4. Steroid Effect by Endotype (Linear Contrasts) -------------------------

coefs    <- coef(coxfit)
vcov_mat <- vcov(coxfit)

# L-fHP (Endotype = 0, reference): steroid effect = coef["Steroid"]
log_hr_lfhp <- coefs["Steroid"]
se_lfhp     <- sqrt(vcov_mat["Steroid", "Steroid"])

# N-fHP (Endotype = 1): steroid effect = coef["Steroid"] + coef["Steroid:Endotype"]
log_hr_nfhp <- coefs["Steroid"] + coefs["Steroid:Endotype"]
var_nfhp    <- vcov_mat["Steroid", "Steroid"] +
               vcov_mat["Steroid:Endotype", "Steroid:Endotype"] +
               2 * vcov_mat["Steroid", "Steroid:Endotype"]
se_nfhp     <- sqrt(var_nfhp)


# --- 5. Hypothesis Tests via glht ---------------------------------------------

contrast_lfhp <- glht(coxfit, linfct = c("Steroid = 0"))
contrast_nfhp <- glht(coxfit, linfct = c("Steroid + `Steroid:Endotype` = 0"))

p_lfhp <- summary(contrast_lfhp)$test$pvalues
p_nfhp <- summary(contrast_nfhp)$test$pvalues


# --- 6. Assemble Results Table ------------------------------------------------

steroid_effect <- data.frame(
  Endotype = c("L-fHP", "N-fHP"),
  HR       = c(exp(log_hr_lfhp), exp(log_hr_nfhp)),
  SE       = c(se_lfhp, se_nfhp),
  Lower    = c(exp(log_hr_lfhp - 1.96 * se_lfhp),
               exp(log_hr_nfhp - 1.96 * se_nfhp)),
  Upper    = c(exp(log_hr_lfhp + 1.96 * se_lfhp),
               exp(log_hr_nfhp + 1.96 * se_nfhp)),
  p_value  = c(p_lfhp, p_nfhp)
)

print(steroid_effect)

write.csv(steroid_effect,
          file      = "steroid_effect_by_endotype.csv",
          row.names = FALSE)
cat("Saved: steroid_effect_by_endotype.csv\n")


# --- 7. Interaction Plot ------------------------------------------------------

p <- ggplot(steroid_effect, aes(x = Endotype, y = HR, group = 1)) +
  geom_line(linewidth = 1,   color = "darkblue") +
  geom_point(size = 4,       color = "darkblue") +
  geom_errorbar(
    aes(ymin = Lower, ymax = Upper),
    width = 0.1, color = "darkblue"
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    x        = "Endotype",
    y        = "TFS Hazard Ratio of Steroid",
    title    = "Adjusted TFS HR of Steroid–Endotype Interaction",
    subtitle = "Significant interaction (p = 0.048)"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_text(size = 18, face = "plain")
  )

print(p)


# --- 8. Save Outputs ----------------------------------------------------------

ggsave("interaction_plot_steroid_endotype.png", p,
       width = 7, height = 5, dpi = 300)

ggsave("interaction_plot_steroid_endotype.pdf", p,
       width = 7, height = 5)

cat("Saved: interaction_plot_steroid_endotype.png and .pdf\n")


# --- Session Info (for reproducibility) ---------------------------------------

sessionInfo()
