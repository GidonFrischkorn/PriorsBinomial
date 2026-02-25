# ==============================================================================
# Script 04b: Post-process BF Calibration RE Results
# ==============================================================================
#
# PURPOSE:
#   Load per-condition results from SimDesign (script 03b), combine into a
#   clean data frame, and save for use in the BF calibration RE Quarto report.
#
# INPUT:  output/res_bf_calibration_re.rds
# OUTPUT: output/bf_calibration_re_summaries.rds
#         (tidy data frame, one row per condition)
#
# ==============================================================================

library(SimDesign)
library(dplyr)
library(here)

load(here("output", "res_bf_calibration_re.rds"))

# Extract summary statistics (one row per condition)
# Stats are stored directly as columns on the SimDesign results object
design_cols <- c("sd_prior_re", "true_b1", "n_subjects", "n_trials")
stat_cols   <- c("mean_log_BF10", "P_BF10_gt3", "P_BF10_gt10",
                 "P_BF01_gt3", "P_BF01_gt10", "n_failed", "n_valid",
                 "mean_true_p0")

summaries <- select(res, all_of(c(design_cols, stat_cols)))

saveRDS(summaries, here("output", "bf_calibration_re_summaries.rds"))
message(
  "Saved: output/bf_calibration_re_summaries.rds  (",
  nrow(summaries), " conditions)"
)
