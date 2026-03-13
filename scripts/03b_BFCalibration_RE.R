# ==============================================================================
# Script 03b: BF Calibration — Within-Subjects with Random Intercepts & Slopes
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
  sd_prior_re = c("default", "gamma", "exponential"),
  true_b1     = c(0.00, 0.05, 0.10, 0.20, 0.50),
  n_subjects  = c(30, 60, 100),
  n_trials    = c(20, 50)
)

# ------------------------------------------------------------------------------
# SimDesign functions
# ------------------------------------------------------------------------------

Generate <- function(condition, fixed_objects = NULL) {
  Attach(condition)
  true_p0       <- runif(1, fixed_objects$b0_range[1], fixed_objects$b0_range[2])
  true_b0       <- apply_link(true_p0, "logit")
  true_sd_re    <- fixed_objects$true_sd_re
  true_sd_slope <- fixed_objects$true_sd_slope
  g_inv         <- function(x) apply_inverse_link(x, link = "logit")

  # Per-subject random intercept and random slope
  u_i <- rnorm(n_subjects, mean = 0, sd = true_sd_re)
  v_i <- rnorm(n_subjects, mean = 0, sd = true_sd_slope)

  data.frame(
    subject_id = rep(seq_len(n_subjects), times = 2),
    y          = c(
      rbinom(n_subjects, n_trials, g_inv(true_b0 + u_i + (true_b1 + v_i))),
      rbinom(n_subjects, n_trials, g_inv(true_b0 + u_i - (true_b1 + v_i)))
    ),
    n         = n_trials,
    condition = rep(c(1L, -1L), each = n_subjects),
    true_p0   = true_p0
  )
}

Analyse <- function(condition, dat, fixed_objects = NULL) {
  Attach(condition)
  sd_b0   <- fixed_objects$sd_b0
  true_p0 <- dat$true_p0[1]

  bf10 <- tryCatch(
    fit_and_get_bf_sd_re(
      dat          = dat,
      link         = "logit",
      dist_b0      = "logistic",
      sd_b0        = sd_b0,
      dist_b1      = "normal",
      sd_b1        = 0.25,
      sd_prior_re  = sd_prior_re,
      random_slope = TRUE
    ),
    error = function(e) NA_real_
  )

  c(BF10 = bf10, BF01 = 1 / bf10, true_p0 = true_p0)
}

Summarise <- function(condition, results, fixed_objects = NULL) {
  bf10    <- results[, "BF10"]
  bf01    <- results[, "BF01"]

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

out_file <- here("output", "res_bf_calibration_re.rds")

fixed_objects_re <- list(
  b0_range      = c(0.4, 0.9),
  sd_b0         = 0.75,
  true_sd_re    = 0.25,
  true_sd_slope = 0.15
)

if (smoke_test) {
  test_design <- Design[Design$sd_prior_re == "gamma", , drop = FALSE][1, ]
  res_test <- runSimulation(
    design        = test_design,
    replications  = 2,
    generate      = Generate,
    analyse       = Analyse,
    summarise     = Summarise,
    fixed_objects = fixed_objects_re,
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
    fixed_objects = fixed_objects_re,
    save_results  = TRUE,
    save_details  = list(
      safe                  = TRUE,
      out_rootdir           = here("output"),
      save_results_dirname  = "Simulation_BFCalibration_RE",
      save_results_filename = "BFCalib_RE_Cond"
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
