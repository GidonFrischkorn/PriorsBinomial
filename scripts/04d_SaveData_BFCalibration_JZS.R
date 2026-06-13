# ==============================================================================
# Script 04d: Post-process BF Calibration — JZS Comparison Results
# ==============================================================================

library(SimDesign)
library(dplyr)
library(tidyr)
library(here)

load(here("output", "res_bf_calibration_jzs_v3.rds"))

# --- Summary statistics -------------------------------------------------------

design_cols <- c("sd_b1", "true_b1", "n_subjects", "n_trials")

stat_cols <- c(
  paste0(c("median_log10_BF10", "mean_log_BF10",
           "P_BF10_gt3", "P_BF10_gt10",
           "P_BF01_gt3", "P_BF01_gt10",
           "n_failed", "n_valid"), "_logistic"),
  paste0(c("median_log10_BF10", "mean_log_BF10",
           "P_BF10_gt3", "P_BF10_gt10",
           "P_BF01_gt3", "P_BF01_gt10",
           "n_failed", "n_valid"), "_cauchy"),
  "mean_true_p0"
)

summaries <- select(res, all_of(c(design_cols, stat_cols)))

saveRDS(summaries, here("output", "bf_calibration_jzs_summaries.rds"))
message("Saved: output/bf_calibration_jzs_summaries.rds  (",
        nrow(summaries), " conditions)")

# --- Operating-characteristic rates table (long format, one row per prior × condition) --

rates_table <- summaries |>
  pivot_longer(
    cols      = matches("^(P_BF10_gt3|P_BF01_gt3|n_valid)_(logistic|cauchy)$"),
    names_to  = c(".value", "prior"),
    names_sep = "_(?=logistic|cauchy)"
  ) |>
  mutate(
    prior_label = case_when(
      prior == "logistic" ~ paste0("Logistic(0, ", sd_b1, ")"),
      prior == "cauchy"   ~ paste0("Cauchy(0, ",   sd_b1, ")")
    )
  ) |>
  select(prior_label, sd_b1, prior, true_b1, n_subjects, n_trials,
         P_BF10_gt3, P_BF01_gt3, n_valid)

saveRDS(rates_table, here("output", "bf_calibration_jzs_rates_table.rds"))
message("Saved: output/bf_calibration_jzs_rates_table.rds")

# --- Load per-condition raw BF10 values ---------------------------------------

results_dir <- here("output", "Simulation_BFCalibration_JZS_v3")

if (!dir.exists(results_dir)) {
  message("Per-condition files not found. Skipping raw BF computation.")
  quit(save = "no")
}

n_cond    <- nrow(summaries)
cont_rows <- list()
raw_rows  <- vector("list", n_cond)

for (i in seq_len(n_cond)) {
  cond_file_qs  <- file.path(results_dir, paste0("BFCalib_JZS_v3_Cond-", i))
  cond_file_rds <- file.path(results_dir,
                             paste0("BFCalib_JZS_v3_Cond-", i, ".rds"))

  if (file.exists(cond_file_qs)) {
    cond_data <- qs2::qd_read(cond_file_qs)
  } else if (file.exists(cond_file_rds)) {
    cond_data <- readRDS(cond_file_rds)
  } else {
    warning("Missing: ", cond_file_qs, " (or .rds)")
    next
  }

  cond_df <- as.data.frame(as.list(cond_data$condition))

  # Continuous quantile summaries — one row per condition per prior
  for (prior in c("logistic", "cauchy")) {
    col      <- paste0("BF10_", prior)
    bf10_vec <- cond_data$results[, col]
    bf10_val <- bf10_vec[!is.na(bf10_vec)]
    log10_bf <- log10(bf10_val)

    cont_rows[[length(cont_rows) + 1L]] <- cbind(
      cond_df,
      data.frame(
        prior             = prior,
        median_log10_BF10 = median(log10_bf),
        q25_log10_BF10    = unname(quantile(log10_bf, 0.25)),
        q75_log10_BF10    = unname(quantile(log10_bf, 0.75)),
        n_valid           = length(bf10_val)
      )
    )
  }

  # Raw replication-level data — long format (one row per replicate per prior)
  raw_rows[[i]] <- cond_df |>
    slice(rep(1L, nrow(cond_data$results))) |>
    bind_cols(as.data.frame(
      cond_data$results[, c("BF10_logistic", "BF10_cauchy", "true_p0")]
    )) |>
    pivot_longer(
      cols      = c(BF10_logistic, BF10_cauchy),
      names_to  = "prior",
      names_prefix = "BF10_",
      values_to = "BF10"
    ) |>
    mutate(log10_BF10 = log10(BF10))
}

continuous_summaries <- bind_rows(cont_rows)
saveRDS(continuous_summaries,
        here("output", "bf_calibration_jzs_continuous.rds"))
message("Saved: output/bf_calibration_jzs_continuous.rds  (",
        nrow(continuous_summaries), " rows)")

raw_df <- bind_rows(raw_rows)
saveRDS(raw_df, here("output", "bf_calibration_jzs_raw.rds"))
message("Saved: output/bf_calibration_jzs_raw.rds  (",
        nrow(raw_df), " replications)")
