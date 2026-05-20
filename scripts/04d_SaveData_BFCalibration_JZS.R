# ==============================================================================
# Script 04d: Post-process BF Calibration — JZS Comparison Results
# ==============================================================================

library(SimDesign)
library(dplyr)
library(here)

load(here("output", "res_bf_calibration_jzs.rds"))

# --- Summary statistics -------------------------------------------------------

design_cols <- c("dist_b1", "sd_b1", "true_b1", "n_subjects", "n_trials")
stat_cols   <- c("median_log10_BF10", "mean_log_BF10",
                 "P_BF10_gt3", "P_BF10_gt10",
                 "P_BF01_gt3", "P_BF01_gt10",
                 "n_failed", "n_valid", "mean_true_p0")

summaries <- select(res, all_of(c(design_cols, stat_cols)))

saveRDS(summaries, here("output", "bf_calibration_jzs_summaries.rds"))
message("Saved: output/bf_calibration_jzs_summaries.rds  (",
        nrow(summaries), " conditions)")

# --- Power / Type I error table -----------------------------------------------
# For each prior condition, summarise P(BF > 3 | H0) and P(BF > 3 | H1)

power_table <- summaries |>
  mutate(
    prior_label = case_when(
      dist_b1 == "logistic" ~ "Matched Logistic(0, 0.25)",
      dist_b1 == "cauchy"   ~ "JZS Cauchy(0, √2/2)"
    )
  ) |>
  select(prior_label, true_b1, n_subjects, n_trials,
         P_BF10_gt3, P_BF01_gt3, n_valid)

saveRDS(power_table, here("output", "bf_calibration_jzs_power_table.rds"))
message("Saved: output/bf_calibration_jzs_power_table.rds")

# --- Load per-condition raw BF10 values ---------------------------------------

results_dir <- here("output", "Simulation_BFCalibration_JZS")

if (!dir.exists(results_dir)) {
  message("Per-condition files not found. Skipping raw BF computation.")
  quit(save = "no")
}

n_cond    <- nrow(summaries)
cont_rows <- vector("list", n_cond)
raw_rows  <- vector("list", n_cond)

for (i in seq_len(n_cond)) {
  # SimDesign >= 2.25 saves via qs2::qd_write (no .rds extension)
  cond_file_qs  <- file.path(results_dir,
                             paste0("BFCalib_JZS_Cond-", i))
  cond_file_rds <- file.path(results_dir,
                             paste0("BFCalib_JZS_Cond-", i, ".rds"))

  if (file.exists(cond_file_qs)) {
    cond_data <- qs2::qd_read(cond_file_qs)
  } else if (file.exists(cond_file_rds)) {
    cond_data <- readRDS(cond_file_rds)
  } else {
    warning("Missing: ", cond_file_qs)
    next
  }
  bf10_vec   <- cond_data$results[, "BF10"]
  bf10_valid <- bf10_vec[!is.na(bf10_vec)]
  log10_bf   <- log10(bf10_valid)

  cont_rows[[i]] <- c(
    as.list(cond_data$condition),
    list(
      median_log10_BF10 = median(log10_bf),
      q25_log10_BF10    = unname(quantile(log10_bf, 0.25)),
      q75_log10_BF10    = unname(quantile(log10_bf, 0.75)),
      n_valid           = length(bf10_valid)
    )
  )

  raw_rows[[i]] <- cbind(
    as.data.frame(as.list(cond_data$condition)),
    data.frame(BF10 = bf10_vec, log10_BF10 = log10(bf10_vec))
  )
}

continuous_summaries <- bind_rows(cont_rows)
saveRDS(continuous_summaries,
        here("output", "bf_calibration_jzs_continuous.rds"))
message("Saved: output/bf_calibration_jzs_continuous.rds  (",
        nrow(continuous_summaries), " conditions)")

raw_df <- bind_rows(raw_rows)
saveRDS(raw_df, here("output", "bf_calibration_jzs_raw.rds"))
message("Saved: output/bf_calibration_jzs_raw.rds  (",
        nrow(raw_df), " replications)")
