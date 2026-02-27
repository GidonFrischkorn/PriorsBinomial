# ==============================================================================
# Script 04: Post-process BF Calibration Results
# ==============================================================================
#
# PURPOSE:
#   Load per-condition results from SimDesign, combine into a clean data frame,
#   and save for use in the BF calibration Quarto report.
#   Additionally, compute Bayesian posteriors (Beta-Binomial and
#   Dirichlet-Multinomial) for detection-rate uncertainty.
#
# INPUT:  output/res_bf_calibration.rds
#         output/Simulation_BFCalibration/BFCalib_Cond-{i}.rds
# OUTPUT: output/bf_calibration_summaries.rds  (tidy data frame, one row per condition)
#         output/bf_calibration_posteriors.rds  (one row per condition x threshold)
#
# ==============================================================================

library(SimDesign)
library(dplyr)
library(here)

source(here("R", "bf_helpers.R"))

load(here("output", "res_bf_calibration.rds"))

# --- Original summaries (backward compatible) --------------------------------

design_cols <- c("dist_b0", "dist_b1", "sd_b1",
                 "true_b1", "n_subjects", "n_trials")
stat_cols   <- c("mean_log_BF10", "P_BF10_gt3", "P_BF10_gt10",
                 "P_BF01_gt3", "P_BF01_gt10", "n_failed", "n_valid",
                 "mean_true_p0")

summaries <- select(res, all_of(c(design_cols, stat_cols)))

saveRDS(summaries, here("output", "bf_calibration_summaries.rds"))
message("Saved: output/bf_calibration_summaries.rds  (", nrow(summaries), " conditions)")

# --- Continuous BF summaries (median, percentiles on log10 scale) ------------

results_dir_cont <- here("output", "Simulation_BFCalibration")

if (!dir.exists(results_dir_cont)) {
  message("Per-condition files not found in ", results_dir_cont,
          ". Skipping continuous BF computation.")
} else {
  n_cond <- nrow(summaries)
  cont_rows <- vector("list", n_cond)

  for (i in seq_len(n_cond)) {
    cond_file <- file.path(results_dir_cont, paste0("BFCalib_Cond-", i, ".rds"))
    if (!file.exists(cond_file)) {
      warning("Missing: ", cond_file)
      next
    }
    cond_data   <- readRDS(cond_file)
    bf10_vec    <- cond_data$results[, "BF10"]
    bf10_valid  <- bf10_vec[!is.na(bf10_vec)]
    log10_bf    <- log10(bf10_valid)

    cont_rows[[i]] <- c(
      as.list(cond_data$condition),
      list(
        median_log10_BF10 = median(log10_bf),
        mean_log10_BF10   = mean(log10_bf),
        q10_log10_BF10    = unname(quantile(log10_bf, 0.10)),
        q25_log10_BF10    = unname(quantile(log10_bf, 0.25)),
        q75_log10_BF10    = unname(quantile(log10_bf, 0.75)),
        q90_log10_BF10    = unname(quantile(log10_bf, 0.90)),
        n_valid           = length(bf10_valid)
      )
    )
  }

  continuous_summaries <- bind_rows(cont_rows)
  saveRDS(continuous_summaries, here("output", "bf_calibration_continuous.rds"))
  message("Saved: output/bf_calibration_continuous.rds  (",
          nrow(continuous_summaries), " conditions)")
}

# --- Bayesian posteriors for detection-rate uncertainty -----------------------

results_dir  <- here("output", "Simulation_BFCalibration")

if (!dir.exists(results_dir)) {
  message("Per-condition files not found in ", results_dir,
          ". Skipping posterior computation.")
} else {
  n_conditions <- nrow(summaries)
  thresholds   <- c(3, 10)

  posterior_rows <- vector("list", n_conditions * length(thresholds))
  idx <- 0L

  for (i in seq_len(n_conditions)) {
    cond_file <- file.path(results_dir, paste0("BFCalib_Cond-", i, ".rds"))
    if (!file.exists(cond_file)) {
      warning("Missing: ", cond_file)
      next
    }
    cond_data <- readRDS(cond_file)
    bf10_vec  <- cond_data$results[, "BF10"]
    design_i  <- cond_data$condition

    for (thr in thresholds) {
      idx <- idx + 1L
      counts <- compute_bf_counts(bf10_vec, threshold = thr)

      power_post <- beta_posterior_summary(counts["n_h1"], counts["n_valid"])
      spec_post  <- beta_posterior_summary(counts["n_h0"], counts["n_valid"])
      dir_post   <- dirichlet_posterior_summary(
        counts["n_h1"], counts["n_h0"], counts["n_inconclusive"]
      )

      posterior_rows[[idx]] <- c(
        as.list(design_i),
        list(threshold = thr),
        as.list(counts),
        setNames(as.list(power_post), paste0("power_", names(power_post))),
        setNames(as.list(spec_post),  paste0("spec_",  names(spec_post))),
        as.list(dir_post)
      )
    }
  }

  posteriors <- bind_rows(posterior_rows)

  saveRDS(posteriors, here("output", "bf_calibration_posteriors.rds"))
  message("Saved: output/bf_calibration_posteriors.rds  (",
          nrow(posteriors), " rows = ",
          n_conditions, " conditions x ", length(thresholds), " thresholds)")
}
