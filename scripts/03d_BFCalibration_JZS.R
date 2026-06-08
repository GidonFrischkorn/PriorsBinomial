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

# Both logistic and cauchy priors are fit on the same data within each
# condition, so dist_b1 is not a design factor — sd_b1 is.
jzs_scale <- round(sqrt(2) / 2, 3)  # ≈ 0.707

Design <- createDesign(
  sd_b1      = c(0.25, jzs_scale),
  true_b1    = c(0.00, 0.10, 0.20),
  n_subjects = c(20, 30, 50),
  n_trials   = c(20, 50)
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

  # Within-subjects: each subject contributes one block of n_trials per condition.
  # Independent random intercept and random slope drawn once per subject.
  true_sd_slope <- fixed_objects$true_sd_slope
  u_int   <- rnorm(n_subjects, mean = 0, sd = true_sd_re)
  u_slope <- rnorm(n_subjects, mean = 0, sd = true_sd_slope)
  cond    <- rep(c(1L, -1L), each = n_subjects)

  p_subj <- g_inv(true_b0 + true_b1 * cond +
                    rep(u_int,   times = 2) +
                    rep(u_slope, times = 2) * cond)

  data.frame(
    subject_id = rep(seq_len(n_subjects), times = 2),
    y          = rbinom(2L * n_subjects, n_trials, p_subj),
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

  fit_bf <- function(dist_b1) {
    bf10 <- tryCatch(
      fit_and_get_bf_sd_re(
        dat          = dat,
        link         = "logit",
        dist_b0      = "logistic",
        sd_b0        = sd_b0,
        dist_b1      = dist_b1,
        sd_b1        = sd_b1,
        sd_prior_re  = sd_prior_re,
        random_slope = TRUE
      ),
      error = function(e) NA_real_
    )
    # Cap Inf (KDE underflow, extreme evidence for H1) at 1e30
    if (is.na(bf10) || is.infinite(bf10)) bf10 <- 1e30
    bf10
  }

  bf10_logistic <- fit_bf("logistic")
  bf10_cauchy   <- fit_bf("cauchy")

  c(
    BF10_logistic = bf10_logistic,
    BF10_cauchy   = bf10_cauchy,
    BF01_logistic = 1 / bf10_logistic,
    BF01_cauchy   = 1 / bf10_cauchy,
    true_p0       = true_p0
  )
}

Summarise <- function(condition, results, fixed_objects = NULL) {
  summarise_bf <- function(bf10, suffix) {
    bf01 <- 1 / bf10
    setNames(
      c(
        median(log10(bf10), na.rm = TRUE),
        mean(log(bf10),     na.rm = TRUE),
        mean(bf10 > 3,      na.rm = TRUE),
        mean(bf10 > 10,     na.rm = TRUE),
        mean(bf01 > 3,      na.rm = TRUE),
        mean(bf01 > 10,     na.rm = TRUE),
        sum(is.na(bf10)),
        sum(!is.na(bf10))
      ),
      paste0(c("median_log10_BF10", "mean_log_BF10",
               "P_BF10_gt3", "P_BF10_gt10",
               "P_BF01_gt3", "P_BF01_gt10",
               "n_failed", "n_valid"),
             suffix)
    )
  }

  c(
    summarise_bf(results[, "BF10_logistic"], "_logistic"),
    summarise_bf(results[, "BF10_cauchy"],   "_cauchy"),
    mean_true_p0 = mean(results[, "true_p0"], na.rm = TRUE)
  )
}

# ------------------------------------------------------------------------------
# Run simulation
# ------------------------------------------------------------------------------

smoke_test  <- FALSE
force_rerun <- FALSE

out_file <- here("output", "res_bf_calibration_jzs_v3.rds")

if (smoke_test) {
  test_design <- Design[1:4, , drop = FALSE]
  res_test <- runSimulation(
    design        = test_design,
    replications  = 2,
    generate      = Generate,
    analyse       = Analyse,
    summarise     = Summarise,
    fixed_objects = list(b0_range = c(0.4, 0.9), sd_b0 = 0.75,
                         true_sd_re = 0.25, true_sd_slope = 0.10,
                         sd_prior_re = "exponential"),
    parallel      = TRUE,
    ncores        = 4,
    packages      = c("brms", "posterior")
  )
  message("Smoke test complete.")
  print(res_test)

} else if (!file.exists(out_file) || force_rerun) {
  res <- runSimulation(
    design        = Design,
    replications  = 300,
    generate      = Generate,
    analyse       = Analyse,
    summarise     = Summarise,
    fixed_objects = list(b0_range = c(0.4, 0.9), sd_b0 = 0.75,
                         true_sd_re = 0.25, true_sd_slope = 0.10,
                         sd_prior_re = "exponential"),
    save_results  = TRUE,
    save_details  = list(
      safe                  = TRUE,
      out_rootdir           = here("output"),
      save_results_dirname  = "Simulation_BFCalibration_JZS_v3",
      save_results_filename = "BFCalib_JZS_v3_Cond"
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
