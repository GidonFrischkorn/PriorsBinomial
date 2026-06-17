# ==============================================================================
# Script 03e: BF Calibration — Probit Confirmation
# ==============================================================================
# Purpose: Confirm that the recommended matched prior scale transfers from the
# logit link to the probit link. The probability integral transform makes the
# matched prior an anchor on the probability scale, so the recommended scale
# (0.25) should yield comparable BF calibration under either link.
#
# This is a focused CONFIRMATION run, not a full factorial: only the recommended
# matched probit configuration (Normal(0, 0.75) intercept + Normal(0, 0.25)
# effect) is fitted, across the focal effect sizes and two design cells that
# also exist in the logit calibration (so rates compare cell-for-cell).
# Mirrors scripts/03_BFCalibration.R (between-subjects, random intercept).
# ==============================================================================

library(SimDesign)
library(brms)
library(dplyr)
library(here)

source(here("R", "link_functions.R"))
source(here("R", "bf_helpers.R"))

# ------------------------------------------------------------------------------
# Design grid (reduced; matched-probit configuration is fixed, not varied)
# ------------------------------------------------------------------------------
# Cells chosen to overlap the logit calibration grid (n_subjects in {30, 60},
# n_trials in {20, 50}) so detection / null-support rates can be compared
# directly against the matched logit results in bf_calibration_summaries.rds.

Design <- createDesign(
  true_b1    = c(0.00, 0.10, 0.20),
  n_subjects = c(30, 60),
  n_trials   = c(20, 50)
)

# ------------------------------------------------------------------------------
# SimDesign functions
# ------------------------------------------------------------------------------

Generate <- function(condition, fixed_objects = NULL) {
  Attach(condition)

  link <- fixed_objects$link

  # Sample intercept uniformly on the probability scale, convert to link scale.
  # Keeping b0_range on the probability scale (0.4-0.9) makes the data-generating
  # baseline identical to the logit run; only the link transform differs.
  true_p0    <- runif(1, fixed_objects$b0_range[1], fixed_objects$b0_range[2])
  true_b0    <- apply_link(true_p0, link)
  true_sd_re <- fixed_objects$true_sd_re

  g_inv <- function(x) apply_inverse_link(x, link = link)

  # Subject-level random intercepts (between-subject variability)
  n_total <- 2L * n_subjects
  cond    <- rep(c(1L, -1L), each = n_subjects)
  u       <- rnorm(n_total, mean = 0, sd = true_sd_re)

  p_subj <- g_inv(true_b0 + true_b1 * cond + u)

  data.frame(
    subject_id = seq_len(n_total),
    y          = rbinom(n_total, n_trials, p_subj),
    n          = n_trials,
    condition  = cond,
    true_p0    = true_p0
  )
}

Analyse <- function(condition, dat, fixed_objects = NULL) {
  Attach(condition)
  true_p0 <- dat$true_p0[1]

  bf10 <- tryCatch(
    fit_and_get_bf_sd_re(
      dat         = dat,
      link        = fixed_objects$link,        # "probit"
      dist_b0     = fixed_objects$dist_b0,     # "normal"  (matched for probit)
      sd_b0       = fixed_objects$sd_b0,       # 0.75
      dist_b1     = fixed_objects$dist_b1,     # "normal"  (matched for probit)
      sd_b1       = fixed_objects$sd_b1,       # 0.25
      sd_prior_re = fixed_objects$sd_prior_re  # "exponential"
    ),
    error = function(e) NA_real_
  )

  c(BF10 = bf10, BF01 = 1 / bf10, true_p0 = true_p0)
}

Summarise <- function(condition, results, fixed_objects = NULL) {
  bf10 <- results[, "BF10"]
  bf01 <- results[, "BF01"]

  c(
    mean_log_BF10 = mean(log(bf10), na.rm = TRUE),
    P_BF10_gt3    = mean(bf10 > 3,  na.rm = TRUE),
    P_BF10_gt10   = mean(bf10 > 10, na.rm = TRUE),
    P_BF01_gt3    = mean(bf01 > 3,  na.rm = TRUE),
    P_BF01_gt10   = mean(bf01 > 10, na.rm = TRUE),
    n_failed      = sum(is.na(bf10)),
    n_valid       = sum(!is.na(bf10)),
    mean_true_p0  = mean(results[, "true_p0"], na.rm = TRUE)
  )
}

# ------------------------------------------------------------------------------
# Run simulation
# ------------------------------------------------------------------------------

smoke_test  <- FALSE
force_rerun <- FALSE

# Matched probit configuration held fixed across all conditions
fixed_cfg <- list(
  link        = "probit",
  dist_b0     = "normal",
  dist_b1     = "normal",
  sd_b0       = 0.75,
  sd_b1       = 0.25,
  sd_prior_re = "exponential",
  b0_range    = c(0.4, 0.9),
  true_sd_re  = 0.25
)

out_file     <- here("output", "res_bf_calibration_probit.rds")
summary_file <- here("output", "bf_calibration_probit_summaries.rds")

if (smoke_test) {
  res_test <- runSimulation(
    design        = Design[1, , drop = FALSE],
    replications  = 2,
    generate      = Generate,
    analyse       = Analyse,
    summarise     = Summarise,
    fixed_objects = fixed_cfg,
    parallel      = FALSE,
    packages      = c("brms", "posterior")
  )
  message("Smoke test complete.")
  print(res_test)

} else if (!file.exists(out_file) || force_rerun) {
  res <- runSimulation(
    design        = Design,
    replications  = 100,
    generate      = Generate,
    analyse       = Analyse,
    summarise     = Summarise,
    fixed_objects = fixed_cfg,
    parallel      = TRUE,
    ncores        = parallel::detectCores() - 2,
    packages      = c("brms", "posterior")
  )
  save(res, file = out_file)
  message("Simulation complete. Results saved to: ", out_file)

  # --- Flat summary table for the supplement chunk --------------------------
  design_cols <- c("true_b1", "n_subjects", "n_trials")
  stat_cols   <- c("mean_log_BF10", "P_BF10_gt3", "P_BF10_gt10",
                   "P_BF01_gt3", "P_BF01_gt10", "n_failed", "n_valid",
                   "mean_true_p0")
  probit_summaries <- dplyr::select(res, dplyr::all_of(c(design_cols, stat_cols)))
  saveRDS(probit_summaries, summary_file)
  message("Saved: ", summary_file, "  (", nrow(probit_summaries), " conditions)")

} else {
  message("Results file already exists. Load with: load('", out_file, "')")
}
