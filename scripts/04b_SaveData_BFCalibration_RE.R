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
summaries <- SimExtract(res, what = "results")

# Add design condition columns
design_cols <- c("sd_prior_re", "true_b1", "n_subjects", "n_trials")
summaries <- bind_cols(
  select(res, all_of(design_cols)),
  summaries
)

saveRDS(summaries, here("output", "bf_calibration_re_summaries.rds"))
message(
  "Saved: output/bf_calibration_re_summaries.rds  (",
  nrow(summaries), " conditions)"
)
