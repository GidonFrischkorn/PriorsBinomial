# ==============================================================================
# Script 04: Post-process BF Calibration Results
# ==============================================================================
#
# PURPOSE:
#   Load per-condition results from SimDesign, combine into a clean data frame,
#   and save for use in the BF calibration Quarto report.
#
# INPUT:  output/res_bf_calibration.rds
# OUTPUT: output/bf_calibration_summaries.rds  (tidy data frame, one row per condition)
#
# ==============================================================================

library(SimDesign)
library(dplyr)
library(here)

load(here("output", "res_bf_calibration.rds"))

# Extract summary statistics (one row per condition)
summaries <- SimExtract(res, what = "results")

# Add design condition columns
design_cols <- c("dist_b0", "dist_b1", "sd_b1",
                 "true_b1", "n_subjects", "n_trials")
summaries <- bind_cols(
  select(res, all_of(design_cols)),
  summaries
)

saveRDS(summaries, here("output", "bf_calibration_summaries.rds"))
message("Saved: output/bf_calibration_summaries.rds  (", nrow(summaries), " conditions)")
