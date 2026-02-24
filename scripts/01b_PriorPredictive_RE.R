# ==============================================================================
# Script 01b: Prior Predictive Distribution — Random Effects Hyperprior Comparison
# ==============================================================================
#
# PURPOSE:
#   Extend Goal 1 to a hierarchical (random-effects) model. Instead of fixing
#   the between-subject SD, we draw it from each candidate hyperprior. This
#   gives the true hierarchical prior predictive:
#
#     sd_re  ~ hyperprior          (one per "study")
#     u_i    ~ Normal(0, sd_re)    (one per subject)
#     b0, b1 ~ fixed-effects priors
#     p_ij   = g^{-1}(b0 + u_i ± (b1 + v_i))
#
#   The three candidate hyperpriors correspond to the levels in
#   make_brms_sd_prior(): brms default student_t(3, 0, 2.5), gamma(2, 4),
#   and exponential(4). The same hyperprior family is used for both sd_re
#   and sd_slope (when re_structure = "intercept_slope").
#
#   Fixed-effects priors are fixed to the recommended values from script 01:
#   logit link, Logistic(0, 0.75) intercept, Normal(0, 0.25) effect.
#
# DESIGN (6 conditions x 10 reps):
#   sd_prior_re  in {"student_t", "gamma", "exponential"}
#   re_structure in {"intercept_only", "intercept_slope"}
#
# FIXED:
#   link       = "logit"
#   dist_b0    = "logistic",  sd_b0 = 0.75
#   dist_b1    = "normal",    sd_b1 = 0.25
#   n_studies  = 1000     (study-level prior draws per replication)
#   n_subjects = 50       (subjects per study — needed for within-study stats)
#
# HYPERPRIOR SAMPLING (positive half-distributions, matching brms parameterisation):
#   student_t   ~ |t(df=3)| × 2.5   [half-Student-t, brms default]
#   gamma       ~ Gamma(shape=2, rate=4)   [mean=0.5, P(sd<0.5)~.71]
#   exponential ~ Exponential(rate=4)      [mean=0.25, P(sd<0.5)~.86]
#
# OUTPUT:
#   output/Simulation_PriorPredictive_RE/PriorPred_RE_Cond_*.rds
#   output/res_prior_predictive_re.rds
# ==============================================================================

library(SimDesign)
library(here)

source(here("R", "link_functions.R"))
source(here("R", "prior_sampling.R"))


# ------------------------------------------------------------------------------
# Helper: draw sd values from a named hyperprior
# ------------------------------------------------------------------------------

#' Draw random-effect SDs from a candidate hyperprior
#'
#' @param n           Integer. Number of draws.
#' @param sd_prior_re Character. One of "student_t", "gamma", "exponential".
#' @return Numeric vector of length n, all >= 0.
sample_sd_hyperprior <- function(n, sd_prior_re) {
  switch(sd_prior_re,
    student_t   = abs(rt(n, df = 3)) * 2.5,   # half-t(3, 0, 2.5)
    gamma       = rgamma(n, shape = 2, rate = 4),
    exponential = rexp(n, rate = 4)
  )
}


# ------------------------------------------------------------------------------
# Design grid
# ------------------------------------------------------------------------------

Design <- createDesign(
  sd_prior_re  = c("student_t", "gamma", "exponential"),
  re_structure = c("intercept_only", "intercept_slope")
)
# Total: 3 x 2 = 6 conditions


# ------------------------------------------------------------------------------
# SimDesign functions
# ------------------------------------------------------------------------------

#' Generate: draw hierarchical prior predictive dataset
#'
#' For each of n_studies "studies", draws study-level parameters (b0, b1,
#' sd_re) and subject-level random effects (u_i, and v_i if intercept_slope).
#'
#' @param condition    One row of Design.
#' @param fixed_objects List with: n_studies, n_subjects, dist_b0, sd_b0,
#'                     dist_b1, sd_b1.
#' @return data.frame with n_studies * n_subjects rows and columns:
#'   study_id, b0, b1, sd_re, sd_slope, u, v.
Generate <- function(condition, fixed_objects = NULL) {
  Attach(condition)
  n_studies  <- fixed_objects$n_studies
  n_subjects <- fixed_objects$n_subjects
  n_total    <- n_studies * n_subjects

  # Study-level fixed-effects draws
  b0 <- sample_prior(n_studies, dist = fixed_objects$dist_b0,
                     scale = fixed_objects$sd_b0)
  b1 <- sample_prior(n_studies, dist = fixed_objects$dist_b1,
                     scale = fixed_objects$sd_b1)

  # Study-level hyperprior draws for sd_re
  sd_re_study <- sample_sd_hyperprior(n_studies, sd_prior_re)

  # For intercept+slope: also draw sd_slope from the same hyperprior
  if (re_structure == "intercept_slope") {
    sd_slope_study <- sample_sd_hyperprior(n_studies, sd_prior_re)
  } else {
    sd_slope_study <- rep(0, n_studies)
  }

  # Expand all study-level quantities to subject level
  study_id     <- rep(seq_len(n_studies), each = n_subjects)
  b0_exp       <- rep(b0,            each = n_subjects)
  b1_exp       <- rep(b1,            each = n_subjects)
  sd_re_exp    <- rep(sd_re_study,   each = n_subjects)
  sd_slope_exp <- rep(sd_slope_study, each = n_subjects)

  # Subject-level random effects
  u <- rnorm(n_total, mean = 0, sd = sd_re_exp)
  # Guard: when re_structure = "intercept_only", sd_slope = 0 everywhere
  v <- if (re_structure == "intercept_slope") {
    rnorm(n_total, mean = 0, sd = sd_slope_exp)
  } else {
    rep(0, n_total)
  }

  data.frame(
    study_id  = study_id,
    b0        = b0_exp,
    b1        = b1_exp,
    sd_re     = sd_re_exp,
    sd_slope  = sd_slope_exp,
    u         = u,
    v         = v
  )
}


#' Analyse: compute marginal and within-study prior predictive summaries
#'
#' @param condition    One row of Design.
#' @param dat          Output of Generate.
#' @param fixed_objects List with: n_subjects (integer).
#' @return Named numeric vector of summary statistics.
Analyse <- function(condition, dat, fixed_objects = NULL) {
  n_subjects <- fixed_objects$n_subjects
  g_inv      <- function(x) apply_inverse_link(x, link = "logit")

  # Subject-level implied probabilities
  # Baseline (intercept prior only, no effect)
  p0_subj     <- g_inv(dat$b0 + dat$u)

  # Condition-specific (marginal subject-level effect)
  linear_pos  <- dat$b0 + dat$u + dat$b1 + dat$v
  linear_neg  <- dat$b0 + dat$u - dat$b1 - dat$v
  p_pos_subj  <- g_inv(linear_pos)
  p_neg_subj  <- g_inv(linear_neg)
  delta_p_subj <- p_pos_subj - p_neg_subj

  # ── Marginal stats (over all n_studies × n_subjects rows) ──────────────────

  # Floor/ceiling of baseline: P(p0 < 0.05 or p0 > 0.95)
  prob_floor_ceiling <- mean(p0_subj < 0.05 | p0_subj > 0.95)

  # Effect size on probability scale
  adp <- abs(delta_p_subj)
  adp_q50 <- as.numeric(quantile(adp, 0.50))
  adp_q75 <- as.numeric(quantile(adp, 0.75))
  adp_q90 <- as.numeric(quantile(adp, 0.90))
  adp_q95 <- as.numeric(quantile(adp, 0.95))
  adp_q99 <- as.numeric(quantile(adp, 0.99))

  # Mass on large effects
  prob_dp_gt10 <- mean(adp > 0.10)
  prob_dp_gt20 <- mean(adp > 0.20)
  prob_dp_gt30 <- mean(adp > 0.30)
  prob_dp_gt50 <- mean(adp > 0.50)

  # ── Hyperprior-level stats (one value per study, then averaged) ────────────

  # Unique sd_re per study (first row of each study is sufficient)
  study_rows    <- seq(1, nrow(dat), by = n_subjects)
  sd_re_studies <- dat$sd_re[study_rows]
  mean_sd_re    <- mean(sd_re_studies)
  p90_sd_re     <- as.numeric(quantile(sd_re_studies, 0.90))

  # Within-study SD of delta_p (effect heterogeneity across subjects)
  # For intercept_only, this should be 0 because v = 0 but u shifts b0 only,
  # so delta_p = g_inv(b0+u+b1) - g_inv(b0+u-b1) still varies across subjects
  # via u. For intercept_slope, v additionally varies b1 per subject.
  study_id        <- dat$study_id
  sd_dp_by_study  <- tapply(delta_p_subj, study_id, sd)
  mean_sd_dp      <- mean(sd_dp_by_study)

  # Within-study sign reversal: proportion of subjects whose delta_p_i
  # has the opposite sign from b1 (averaged across studies)
  b1_sign_study   <- sign(dat$b1[study_rows])
  b1_sign_exp     <- rep(b1_sign_study, each = n_subjects)
  sign_rev_by_subj <- as.integer(sign(delta_p_subj) != b1_sign_exp)
  prop_sign_rev_by_study <- tapply(sign_rev_by_subj, study_id, mean)
  mean_prop_sign_rev <- mean(prop_sign_rev_by_study)

  c(
    prob_floor_ceiling  = prob_floor_ceiling,
    adp_q50             = adp_q50,
    adp_q75             = adp_q75,
    adp_q90             = adp_q90,
    adp_q95             = adp_q95,
    adp_q99             = adp_q99,
    prob_dp_gt10        = prob_dp_gt10,
    prob_dp_gt20        = prob_dp_gt20,
    prob_dp_gt30        = prob_dp_gt30,
    prob_dp_gt50        = prob_dp_gt50,
    mean_sd_re          = mean_sd_re,
    p90_sd_re           = p90_sd_re,
    mean_sd_dp          = mean_sd_dp,
    mean_prop_sign_rev  = mean_prop_sign_rev
  )
}


#' Summarise: average statistics across replications
#'
#' @param condition    One row of Design.
#' @param results      Matrix (replications x statistics).
#' @param fixed_objects Unused.
#' @return Named numeric vector (column means across replications).
Summarise <- function(condition, results, fixed_objects = NULL) {
  colMeans(results)
}


# ------------------------------------------------------------------------------
# Run simulation
# ------------------------------------------------------------------------------

smoke_test <- FALSE

out_file <- here("output", "res_prior_predictive_re.rds")

fixed_objects_re <- list(
  n_studies  = 1000,
  n_subjects = 50,
  dist_b0    = "logistic",
  sd_b0      = 0.75,
  dist_b1    = "normal",
  sd_b1      = 0.25
)

if (smoke_test) {
  res_test <- runSimulation(
    design        = Design[1, , drop = FALSE],
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
    replications  = 10,
    generate      = Generate,
    analyse       = Analyse,
    summarise     = Summarise,
    fixed_objects = fixed_objects_re,
    save_results  = TRUE,
    save_details  = list(
      safe                  = TRUE,
      out_rootdir           = here("output"),
      save_results_dirname  = "Simulation_PriorPredictive_RE",
      save_results_filename = "PriorPred_RE_Cond"
    ),
    parallel = TRUE,
    ncores   = parallel::detectCores() - 2,
    packages = character(0)
  )
  save(res, file = out_file)
  message("Simulation complete. Results saved to: ", out_file)

} else {
  message("Results file already exists. Load with: load('", out_file, "')")
}
