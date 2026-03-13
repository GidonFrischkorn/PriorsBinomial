# ==============================================================================
# Script 01: Prior Predictive Distribution Simulation
# ==============================================================================

library(SimDesign)
library(here)

# Source helper functions
source(here("R", "link_functions.R"))
source(here("R", "prior_sampling.R"))

# ------------------------------------------------------------------------------
# Design grid
# ------------------------------------------------------------------------------

Design <- createDesign(
  link    = c("logit", "probit"),
  dist_b0 = c("logistic", "normal", "cauchy"),
  sd_b0   = c(0.50, 0.75, 1.00, 1.50),
  dist_b1 = c("normal", "logistic", "cauchy"),
  sd_b1   = c(0.10, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.00, 1.50, 2.00)
)

# ------------------------------------------------------------------------------
# SimDesign functions
# ------------------------------------------------------------------------------

Generate <- function(condition, fixed_objects = NULL) {
  Attach(condition)
  n <- fixed_objects$n_draws

  list(
    b0 = sample_prior(n, dist = dist_b0, scale = sd_b0),
    b1 = sample_prior(n, dist = dist_b1, scale = sd_b1)
  )
}

Analyse <- function(condition, dat, fixed_objects = NULL) {
  Attach(condition)

  g_inv   <- apply_inverse_link
  p0      <- g_inv(dat$b0, link = link)
  p1      <- g_inv(dat$b0 + dat$b1, link = link)
  p2      <- g_inv(dat$b0 - dat$b1, link = link)
  delta_p <- p1 - p2

  c(
    # Intercept prior predictive (marginal p)
    p_q05  = as.numeric(quantile(p0, 0.05)),
    p_q25  = as.numeric(quantile(p0, 0.25)),
    p_q50  = as.numeric(quantile(p0, 0.50)),
    p_q75  = as.numeric(quantile(p0, 0.75)),
    p_q95  = as.numeric(quantile(p0, 0.95)),

    # Effect size on probability scale (absolute)
    adp_q50 = as.numeric(quantile(abs(delta_p), 0.50)),
    adp_q75 = as.numeric(quantile(abs(delta_p), 0.75)),
    adp_q90 = as.numeric(quantile(abs(delta_p), 0.90)),
    adp_q95 = as.numeric(quantile(abs(delta_p), 0.95)),
    adp_q99 = as.numeric(quantile(abs(delta_p), 0.99)),

    # Signed effect size quantiles (symmetry check)
    dp_q05  = as.numeric(quantile(delta_p, 0.05)),
    dp_q95  = as.numeric(quantile(delta_p, 0.95)),

    # Floor/ceiling diagnostic: P(p < 0.05 or p > 0.95)
    prob_floor_ceiling = mean(p0 < 0.05 | p0 > 0.95),

    # Effect size calibration
    prob_dp_gt10 = mean(abs(delta_p) > 0.10),
    prob_dp_gt20 = mean(abs(delta_p) > 0.20),
    prob_dp_gt30 = mean(abs(delta_p) > 0.30),
    prob_dp_gt50 = mean(abs(delta_p) > 0.50)
  )
}

Summarise <- function(condition, results, fixed_objects = NULL) {
  colMeans(results)
}


# ------------------------------------------------------------------------------
# Run simulation
# ------------------------------------------------------------------------------

out_file <- here("output", "res_prior_predictive.rds")

if (!file.exists(out_file)) {
  res <- runSimulation(
    design        = Design,
    replications  = 10,
    generate      = Generate,
    analyse       = Analyse,
    summarise     = Summarise,
    fixed_objects = list(n_draws = 10000),
    save_results  = TRUE,
    save_details  = list(
      safe                  = TRUE,
      out_rootdir           = here("output"),
      save_results_dirname  = "Simulation_PriorPredictive",
      save_results_filename = "PriorPred_Cond"
    ),
    parallel = TRUE,
    ncores   = parallel::detectCores() - 2,
    packages = character(0)   # helpers sourced above
  )
  save(res, file = out_file)
  message("Simulation complete. Results saved to: ", out_file)
} else {
  message("Results file already exists. Load with: load('", out_file, "')")
}
