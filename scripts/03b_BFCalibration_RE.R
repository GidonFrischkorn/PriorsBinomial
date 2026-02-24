# ==============================================================================
# Script 03b: BF Calibration — Within-Subjects with Random Intercepts & Slopes
# ==============================================================================
#
# PURPOSE:
#   Extend Goal 2 to a within-subjects design. Each subject contributes one
#   binomial observation per condition. Adds (1 + condition | subject_id) to
#   the brms model: subjects vary both in baseline probability (random
#   intercept) and in how strongly condition affects them (random slope).
#   Varies the SD prior for both random effects to assess its effect on BF
#   calibration (power and specificity).
#
# MODEL:
#   Within-subjects, sum-to-zero contrast x in {+1, -1}.
#   y_{i,c} ~ Binomial(n_trials, p_{i,c})
#   g(p_{i,c}) = (b0 + u_i) + (b1 + v_i) * x_c
#   u_i ~ Normal(0, sigma_u),  v_i ~ Normal(0, sigma_v)
#   corr(u_i, v_i) ~ LKJ(1)  [uniform, brms default]
#   H1: b0 ~ intercept_prior, b1 ~ effect_prior,
#       sigma_u ~ sd_prior_re, sigma_v ~ sd_prior_re
#   H0: Savage-Dickey at b1 = 0 (identical RE structure as H1)
#
# DESIGN (72 conditions x 200 reps):
#   sd_prior_re in {"default", "gamma", "exponential"}
#   true_b1     in {0.00, 0.25, 0.50, 0.75}
#   n_subjects  in {30, 60, 100}
#   n_trials    in {20, 50}
#
# FIXED:
#   link        = "logit"    (matched prior from Goal 1)
#   dist_b0     = "logistic" (matched intercept prior for logit link)
#   dist_b1     = "normal"   (recommended from script 03)
#   sd_b1       = 0.25       (recommended from script 03)
#   true_b0     = 0          (50% baseline at null)
#   sd_b0       = 0.75       (recommended from Goal 1)
#   true_sd_re  = 0.25       (realistic between-subject SD on logit scale,
#                             ~6% between-subject SD in probability at p=0.5)
#   true_sd_slope = 0.15     (realistic between-subject slope SD on logit scale;
#                             smaller than intercept SD — subjects vary less in
#                             their condition sensitivity than in baseline)
#
# REDUCTIONS vs. initial 432-condition grid:
#   dist_b1, sd_b1 fixed to recommended values from 03 — the RE simulation
#   focuses on the SD prior comparison; effect prior variation is already
#   covered in script 03.
#
# SD PRIOR LEVELS (on logit/probit scale):
#   "default"     — student_t(3, 0, 2.5): brms default, very wide
#   "gamma"       — gamma(2, 4): mean 0.5, P(sd < 0.5) ~= 0.71
#   "exponential" — exponential(4): mean 0.25, P(sd < 0.5) ~= 0.86
#
# OUTPUT:
#   output/Simulation_BFCalibration_RE/BFCalib_RE_Cond_*.rds
#   output/res_bf_calibration_re.rds
#
# RUNTIME: Long. Run on compute server with many cores.
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
  sd_prior_re = c("default", "gamma", "exponential"),
  true_b1     = c(0.00, 0.25, 0.50, 0.75),
  n_subjects  = c(30, 60, 100),
  n_trials    = c(20, 50)
)
# Total: 3 x 4 x 3 x 2 = 72 conditions


# ------------------------------------------------------------------------------
# SimDesign functions
# ------------------------------------------------------------------------------

#' Generate: simulate a within-subjects binomial dataset with random intercepts
#' and random slopes
#'
#' Each subject contributes two observations (one per condition). A subject-
#' level random intercept u_i ~ Normal(0, true_sd_re) models between-subject
#' variability in baseline probability; a random slope v_i ~ Normal(0,
#' true_sd_slope) models between-subject variability in the condition effect.
#'
#' @param condition    One row of Design.
#' @param fixed_objects List with: true_b0 (numeric), sd_b0 (numeric),
#'                    true_sd_re (numeric), true_sd_slope (numeric).
#' @return data.frame with columns subject_id, y, n, condition.
Generate <- function(condition, fixed_objects = NULL) {
  Attach(condition)
  true_b0      <- fixed_objects$true_b0
  true_sd_re   <- fixed_objects$true_sd_re
  true_sd_slope <- fixed_objects$true_sd_slope
  g_inv        <- function(x) apply_inverse_link(x, link = "logit")

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
    condition = rep(c(1L, -1L), each = n_subjects)
  )
}


#' Analyse: fit brms RE model with random slopes and compute BF10 via
#' Savage-Dickey ratio
#'
#' @param condition    One row of Design.
#' @param dat          Output of Generate.
#' @param fixed_objects List with: sd_b0 (numeric).
#' @return Named numeric vector with BF10 and BF01.
Analyse <- function(condition, dat, fixed_objects = NULL) {
  Attach(condition)
  sd_b0 <- fixed_objects$sd_b0

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

  c(BF10 = bf10, BF01 = 1 / bf10)
}


#' Summarise: compute power and specificity metrics
#'
#' @param condition    One row of Design.
#' @param results      Matrix (replications x 2) with BF10 and BF01 columns.
#' @param fixed_objects Unused.
#' @return Named numeric vector of summary statistics.
Summarise <- function(condition, results, fixed_objects = NULL) {
  bf10    <- results[, "BF10"]
  bf01    <- results[, "BF01"]
  n_valid <- sum(!is.na(bf10))

  c(
    mean_log_BF10 = mean(log(bf10), na.rm = TRUE),
    P_BF10_gt3    = mean(bf10 > 3,  na.rm = TRUE),
    P_BF10_gt10   = mean(bf10 > 10, na.rm = TRUE),
    P_BF01_gt3    = mean(bf01 > 3,  na.rm = TRUE),
    P_BF01_gt10   = mean(bf01 > 10, na.rm = TRUE),
    n_failed      = sum(is.na(bf10)),
    n_valid       = n_valid
  )
}


# ------------------------------------------------------------------------------
# Run simulation
# ------------------------------------------------------------------------------

# Set smoke_test = TRUE to run a quick check on a single condition with 2 reps
smoke_test <- FALSE

out_file <- here("output", "res_bf_calibration_re.rds")

fixed_objects_re <- list(
  true_b0       = 0,
  sd_b0         = 0.75,
  true_sd_re    = 0.25,
  true_sd_slope = 0.15
)

if (smoke_test) {
  test_design <- Design[
    Design$sd_prior_re == "gamma" & Design$dist_b1 == "normal",
    , drop = FALSE
  ][1, ]
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
