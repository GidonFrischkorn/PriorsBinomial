# ==============================================================================
# Script 03c: BF Calibration — Information & Lindley Consistency
# ==============================================================================
# Purpose: Demonstrate information consistency (BF grows with effect size at
# fixed N=100) and Lindley consistency (BF grows with N at fixed b1=0.20) for
# the matched Logistic(0, 0.25) prior.
#
# BF estimation via bridge sampling (fit_and_get_bf_bs_re). The Savage-Dickey
# KDE approach used in the main calibration simulation is unreliable here:
# when the posterior is far from zero (large effects, large N), the KDE
# bandwidth drives the density at 0 to near-machine-precision zero, inflating
# BF to ~10^15–10^30 even for moderate effects at moderate N.
# ==============================================================================

library(SimDesign)
library(brms)
library(here)

source(here("R", "link_functions.R"))
source(here("R", "bf_helpers.R"))

# ------------------------------------------------------------------------------
# Design grid
# ------------------------------------------------------------------------------

Design <- createDesign(
  dist_b0    = "logistic",     # matched intercept prior only
  dist_b1    = "logistic",     # matched effect prior only
  sd_b1      = 0.25,           # recommended default
  true_b1    = c(0.00, 0.05, 0.10, 0.20, 0.30, 0.50),  # plausible effect range
  n_subjects = c(20, 30, 50, 100, 200),           # Lindley consistency range
  n_trials   = 20              # fixed at sparse scenario
)

# ------------------------------------------------------------------------------
# SimDesign functions
# ------------------------------------------------------------------------------

Generate <- function(condition, fixed_objects = NULL) {
  Attach(condition)

  true_p0    <- runif(1, fixed_objects$b0_range[1], fixed_objects$b0_range[2])
  true_b0    <- apply_link(true_p0, "logit")
  true_sd_re <- fixed_objects$true_sd_re

  g_inv <- function(x) apply_inverse_link(x, link = "logit")

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
  sd_b0       <- fixed_objects$sd_b0
  sd_prior_re <- fixed_objects$sd_prior_re
  true_p0     <- dat$true_p0[1]

  # Bridge sampling gives reliable BF even when posterior is far from zero,
  # avoiding the KDE underflow that inflates Savage-Dickey estimates at large
  # effects / large N.
  bf10 <- tryCatch(
    fit_and_get_bf_bs_re(
      dat         = dat,
      link        = "logit",
      dist_b0     = dist_b0,
      sd_b0       = sd_b0,
      dist_b1     = dist_b1,
      sd_b1       = sd_b1,
      sd_prior_re = sd_prior_re
    ),
    error = function(e) NA_real_
  )

  c(BF10 = bf10, BF01 = 1 / bf10, true_p0 = true_p0)
}

Summarise <- function(condition, results, fixed_objects = NULL) {
  bf10 <- results[, "BF10"]
  bf01 <- results[, "BF01"]

  c(
    median_log10_BF10 = median(log10(bf10), na.rm = TRUE),
    mean_log_BF10     = mean(log(bf10),     na.rm = TRUE),
    P_BF10_gt3        = mean(bf10 > 3,      na.rm = TRUE),
    P_BF10_gt10       = mean(bf10 > 10,     na.rm = TRUE),
    P_BF01_gt3        = mean(bf01 > 3,      na.rm = TRUE),
    P_BF01_gt10       = mean(bf01 > 10,     na.rm = TRUE),
    n_failed          = sum(is.na(bf10)),
    n_valid           = sum(!is.na(bf10)),
    mean_true_p0      = mean(results[, "true_p0"], na.rm = TRUE)
  )
}

# ------------------------------------------------------------------------------
# Run simulation
# ------------------------------------------------------------------------------

smoke_test   <- FALSE
force_rerun  <- TRUE

out_file <- here("output", "res_bf_calibration_consistency.rds")

if (smoke_test) {
  test_design <- Design[1:4, , drop = FALSE]
  res_test <- runSimulation(
    design        = test_design,
    replications  = 2,
    generate      = Generate,
    analyse       = Analyse,
    summarise     = Summarise,
    fixed_objects = list(b0_range = c(0.4, 0.9), sd_b0 = 0.75,
                         true_sd_re = 0.25, sd_prior_re = "exponential"),
    parallel      = TRUE,
    ncores        = 4,
    packages      = c("brms", "posterior")
  )
  message("Smoke test complete.")
  print(res_test)

} else if (!file.exists(out_file) || force_rerun) {
  res <- runSimulation(
    design        = Design,
    replications  = 50,
    generate      = Generate,
    analyse       = Analyse,
    summarise     = Summarise,
    fixed_objects = list(b0_range = c(0.4, 0.9), sd_b0 = 0.75,
                         true_sd_re = 0.25, sd_prior_re = "exponential"),
    save_results  = TRUE,
    save_details  = list(
      safe                  = TRUE,
      out_rootdir           = here("output"),
      save_results_dirname  = "Simulation_BFCalibration_Consistency",
      save_results_filename = "BFCalib_Consistency_Cond"
    ),
    parallel = TRUE,
    ncores   = parallel::detectCores() - 2,
    packages = c("brms", "posterior")
  )
  save(res, file = out_file)
  message("Simulation complete. Results saved to: ", out_file)

} else {
  message("Results file already exists. Load with: load('", out_file, "')")
}
