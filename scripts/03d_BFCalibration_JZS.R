# ==============================================================================
# Script 03d: BF Calibration — JZS/Cauchy Comparison
# ==============================================================================
# Purpose: Compare the matched Logistic(0, 0.25) prior against the JZS/Cauchy
# default Cauchy(0, sqrt(2)/2) as used in the BayesFactor R package.
# ==============================================================================

library(SimDesign)
library(brms)
library(here)

source(here("R", "link_functions.R"))
source(here("R", "bf_helpers.R"))

# ------------------------------------------------------------------------------
# Design grid
# ------------------------------------------------------------------------------

# Create all combinations, then keep only the intended diagonal:
# matched Logistic at sd = 0.25 vs. JZS Cauchy at sd = sqrt(2)/2
jzs_scale <- round(sqrt(2) / 2, 6)  # ≈ 0.707107

Design_full <- createDesign(
  dist_b1    = c("logistic", "cauchy"),
  sd_b1      = c(0.25, jzs_scale),
  true_b1    = c(0.00, 0.10, 0.20, 0.50),
  n_subjects = c(30, 60, 100),
  n_trials   = c(20, 50)
)

Design <- subset(
  Design_full,
  (dist_b1 == "logistic" & sd_b1 == 0.25) |
    (dist_b1 == "cauchy"   & sd_b1 == jzs_scale)
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

  bf10 <- tryCatch(
    fit_and_get_bf_sd_re(
      dat         = dat,
      link        = "logit",
      dist_b0     = "logistic",
      sd_b0       = sd_b0,
      dist_b1     = dist_b1,
      sd_b1       = sd_b1,
      sd_prior_re = sd_prior_re
    ),
    error = function(e) NA_real_
  )

  # Posterior density at 0 can underflow to machine precision for large effects
  # or wide priors (Cauchy). Both NA (Stan error) and Inf (KDE underflow) encode
  # BF -> Inf; cap at 1e30 so Summarise always receives a finite value.
  if (is.na(bf10) || is.infinite(bf10)) bf10 <- 1e30

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

smoke_test  <- FALSE
force_rerun <- FALSE

out_file <- here("output", "res_bf_calibration_jzs.rds")

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
    replications  = 100,
    generate      = Generate,
    analyse       = Analyse,
    summarise     = Summarise,
    fixed_objects = list(b0_range = c(0.4, 0.9), sd_b0 = 0.75,
                         true_sd_re = 0.25, sd_prior_re = "exponential"),
    save_results  = TRUE,
    save_details  = list(
      safe                  = TRUE,
      out_rootdir           = here("output"),
      save_results_dirname  = "Simulation_BFCalibration_JZS",
      save_results_filename = "BFCalib_JZS_Cond"
    ),
    parallel = TRUE,
    ncores   = 10,
    packages = c("brms", "posterior")
  )
  save(res, file = out_file)
  message("Simulation complete. Results saved to: ", out_file)

} else {
  message("Results file already exists. Load with: load('", out_file, "')")
}
