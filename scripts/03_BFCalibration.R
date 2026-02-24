# ==============================================================================
# Script 03: BF Calibration Simulation
# ==============================================================================
#
# PURPOSE:
#   Goal 2: For each combination of prior specification, true effect size, and
#   sample size, simulate binomial datasets and compute Bayes factors to
#   assess power (P(BF10 > threshold | H1)) and type-I specificity
#   (P(BF01 > threshold | H0)).
#
# MODEL:
#   Two-condition design with sum-to-zero contrast x in {+1, -1}.
#   k_i ~ Binomial(n_trials, p_i)
#   g(p_i) = b0 + b1 * x_i
#   H1: b0 ~ matched_prior, b1 ~ effect_prior
#   H0: b0 ~ matched_prior, b1 constrained to 0
#
# DESIGN (288 conditions x 200 reps):
#   dist_b0     in {logistic, normal}      [matched vs. misfit for logit]
#   dist_b1     in {normal, logistic, cauchy}
#   sd_b1       in {0.25, 0.50}
#   true_b1     in {0.00, 0.25, 0.50, 0.75}
#   n_subjects  in {30, 60, 100}
#   n_trials    in {20, 50}
#
# FIXED:
#   link     = "logit" (probit is symmetric by the matched prior principle)
#   b0_range = c(0.4, 0.9): true_b0 is sampled per replication from
#              Uniform(0.4, 0.9) on the probability scale, then transformed
#              to the logit scale. Covers realistic baseline accuracies in
#              cognitive experiments (2-AFC to 4-AFC).
#   sd_b0    = 0.75 (informed by Goal 1)
#
# REDUCTIONS vs. original 1,620-condition grid:
#   link:     fixed to logit (probit analogous by matched prior principle)
#   sd_b1:    dropped 1.00 (too diffuse, established in Goal 1)
#   true_b1:  dropped 1.00 (ceiling effects; trend clear from 0.75)
#   n_trials: dropped 100 (monotone improvement; extrapolates from 20/50)
#
# OUTPUT:
#   output/Simulation_BFCalibration/BFCalib_Cond_*.rds
#   output/res_bf_calibration.rds
#
# RUNTIME: Very long (~hours). Run on compute server with many cores.
#          Test with smoke_test = TRUE before full run.
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
  dist_b0    = c("logistic", "normal"),    # matched vs. misfit intercept prior
  dist_b1    = c("normal", "logistic", "cauchy"),
  sd_b1      = c(0.25, 0.50),
  true_b1    = c(0.00, 0.25, 0.50),
  n_subjects = c(30, 60, 100),
  n_trials   = c(20, 50)
)
# Total: 2 x 3 x 2 x 4 x 3 x 2 = 288 conditions


# ------------------------------------------------------------------------------
# SimDesign functions
# ------------------------------------------------------------------------------

#' Generate: simulate a binomial dataset with subject-level random intercepts
#'
#' @param condition  One row of Design.
#' @param fixed_objects List with: b0_range (numeric[2]), sd_b0 (numeric),
#'   sd_re (numeric) between-subject SD on the logit scale.
#' @return data.frame with columns y, n, condition, subject_id, true_p0.
Generate <- function(condition, fixed_objects = NULL) {
  Attach(condition)

  # Sample intercept uniformly on probability scale, convert to link scale
  true_p0 <- runif(1, fixed_objects$b0_range[1], fixed_objects$b0_range[2])
  true_b0 <- apply_link(true_p0, "logit")
  sd_re   <- fixed_objects$sd_re

  g_inv <- function(x) apply_inverse_link(x, link = "logit")

  # Subject-level random intercepts (between-subject variability)
  n_total <- 2L * n_subjects
  cond    <- rep(c(1L, -1L), each = n_subjects)
  u       <- rnorm(n_total, mean = 0, sd = sd_re)

  p_subj <- g_inv(true_b0 + true_b1 * cond + u)

  data.frame(
    y          = rbinom(n_total, n_trials, p_subj),
    n          = n_trials,
    condition  = cond,
    subject_id = seq_len(n_total),
    true_p0    = true_p0
  )
}


# example_dat <- Generate(condition = Design[1,], fixed_objects = list(b0_range = c(0.4, 0.9), sd_b0 = 0.75, sd_re = 0.25))

#' Analyse: fit brms model with random intercept and compute BF10
#'
#' @param condition  One row of Design.
#' @param dat        Output of Generate.
#' @param fixed_objects List with: sd_b0 (numeric), sd_prior_re (character).
#' @return Named numeric vector with BF10, BF01, and true_p0.
Analyse <- function(condition, dat, fixed_objects = NULL) {
  Attach(condition)
  sd_b0       <- fixed_objects$sd_b0
  sd_prior_re <- fixed_objects$sd_prior_re
  true_p0     <- dat$true_p0[1]

  bf10 <- tryCatch(
    fit_and_get_bf_sd_re(
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


# res <- Analyse(condition = Design[1,],
#                dat = example_dat,
#                fixed_objects = list(b0_range = c(0.4, 0.9), sd_b0 = 0.75, sd_re = 0.25))

#' Summarise: compute power and specificity metrics
#'
#' @param condition  One row of Design.
#' @param results    Matrix (replications x 3) with BF10, BF01, true_p0 columns.
#' @param fixed_objects Unused.
#' @return Named numeric vector of summary statistics.
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

# Set smoke_test = TRUE to run a quick check on a single condition with 2 reps
smoke_test <- FALSE

out_file <- here("output", "res_bf_calibration.rds")

if (smoke_test) {
  test_design <- Design[1, , drop = FALSE]
  res_test <- runSimulation(
    design        = test_design,
    replications  = 2,
    generate      = Generate,
    analyse       = Analyse,
    summarise     = Summarise,
    fixed_objects = list(b0_range = c(0.4, 0.9), sd_b0 = 0.75,
                         sd_re = 0.25, sd_prior_re = "exponential"),
    parallel      = FALSE
  )
  message("Smoke test complete.")
  print(res_test)

} else if (!file.exists(out_file)) {
  res <- runSimulation(
    design        = Design,
    replications  = 100,
    generate      = Generate,
    analyse       = Analyse,
    summarise     = Summarise,
    fixed_objects = list(b0_range = c(0.4, 0.9), sd_b0 = 0.75,
                         sd_re = 0.25, sd_prior_re = "exponential"),
    save_results  = TRUE,
    save_details  = list(
      safe                  = TRUE,
      out_rootdir           = here("output"),
      save_results_dirname  = "Simulation_BFCalibration",
      save_results_filename = "BFCalib_Cond"
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
