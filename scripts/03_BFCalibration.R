# ==============================================================================
# Script 03: BF Calibration Simulation
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
  true_b1    = c(0.00, 0.05, 0.10, 0.20, 0.50),
  n_subjects = c(30, 60, 100),
  n_trials   = c(20, 50)
)

# ------------------------------------------------------------------------------
# SimDesign functions
# ------------------------------------------------------------------------------

Generate <- function(condition, fixed_objects = NULL) {
  Attach(condition)

  # Sample intercept uniformly on probability scale, convert to link scale
  true_p0    <- runif(1, fixed_objects$b0_range[1], fixed_objects$b0_range[2])
  true_b0    <- apply_link(true_p0, "logit")
  true_sd_re <- fixed_objects$true_sd_re

  g_inv <- function(x) apply_inverse_link(x, link = "logit")

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
                         true_sd_re = 0.25, sd_prior_re = "exponential"),
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
                         true_sd_re = 0.25, sd_prior_re = "exponential"),
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
