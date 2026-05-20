# ==============================================================================
# Script 04c: Post-process BF Calibration — Consistency Results
# ==============================================================================
# Builds summaries directly from per-condition files (conditions 1-30,
# n_subjects <= 200). Does not require the main res object, which may be
# incomplete if the simulation crashed on the n_subjects = 400 conditions.
# ==============================================================================

library(dplyr)
library(here)

results_dir <- here("output", "Simulation_BFCalibration_Consistency")

if (!dir.exists(results_dir)) {
  stop("Per-condition directory not found: ", results_dir)
}

# Conditions 1-30 cover n_subjects <= 200 (6 effects x 5 sample sizes).
# Increase n_cond to 36 once the n_subjects = 400 conditions are available.
n_cond <- 30

cont_rows <- vector("list", n_cond)
raw_rows  <- vector("list", n_cond)

for (i in seq_len(n_cond)) {
  cond_file_qs  <- file.path(results_dir,
                             paste0("BFCalib_Consistency_Cond-", i))
  cond_file_rds <- file.path(results_dir,
                             paste0("BFCalib_Consistency_Cond-", i, ".rds"))

  if (file.exists(cond_file_qs)) {
    cond_data <- qs2::qd_read(cond_file_qs)
  } else if (file.exists(cond_file_rds)) {
    cond_data <- readRDS(cond_file_rds)
  } else {
    warning("Missing: condition ", i, " — skipping")
    next
  }

  bf10_vec   <- cond_data$results[, "BF10"]
  bf10_valid <- bf10_vec[!is.na(bf10_vec)]
  log10_bf   <- log10(bf10_valid)

  cont_rows[[i]] <- c(
    as.list(cond_data$condition),
    list(
      median_log10_BF10 = if (length(log10_bf) > 0) median(log10_bf)               else NA_real_,
      q25_log10_BF10    = if (length(log10_bf) > 0) unname(quantile(log10_bf, 0.25)) else NA_real_,
      q75_log10_BF10    = if (length(log10_bf) > 0) unname(quantile(log10_bf, 0.75)) else NA_real_,
      n_valid           = length(bf10_valid),
      n_failed          = sum(is.na(bf10_vec))
    )
  )

  raw_rows[[i]] <- cbind(
    as.data.frame(as.list(cond_data$condition)),
    data.frame(BF10 = bf10_vec, log10_BF10 = log10(bf10_vec))
  )
}

continuous_summaries <- bind_rows(cont_rows)

saveRDS(continuous_summaries,
        here("output", "bf_calibration_consistency_continuous.rds"))
message("Saved: bf_calibration_consistency_continuous.rds  (",
        nrow(continuous_summaries), " conditions)")

# Also save as summaries so consistency_available is TRUE in the QMD
saveRDS(continuous_summaries,
        here("output", "bf_calibration_consistency_summaries.rds"))
message("Saved: bf_calibration_consistency_summaries.rds  (",
        nrow(continuous_summaries), " conditions)")

raw_df <- bind_rows(raw_rows)
saveRDS(raw_df, here("output", "bf_calibration_consistency_raw.rds"))
message("Saved: bf_calibration_consistency_raw.rds  (",
        nrow(raw_df), " replications)")
