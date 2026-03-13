# ==============================================================================
# Script 02b: Post-process RE Prior Predictive Simulation Results
# ==============================================================================

library(SimDesign)
library(dplyr)
library(tidyr)
library(here)

source(here("R", "link_functions.R"))
source(here("R", "prior_sampling.R"))

# ------------------------------------------------------------------------------
# Helper (must match 01b)
# ------------------------------------------------------------------------------

sample_sd_hyperprior <- function(n, hp_family, hp_mean) {
  switch(hp_family,
    exponential = rexp(n, rate = 1 / hp_mean),
    gamma       = rgamma(n, shape = 2, rate = 2 / hp_mean),
    half_normal = abs(rnorm(n, mean = 0, sd = hp_mean / sqrt(2 / pi))),
    half_t      = abs(rt(n, df = 3)) * (hp_mean / (2 * sqrt(3) / pi)),
    stop("Unknown hp_family: ", hp_family)
  )
}

# ------------------------------------------------------------------------------
# Load SimDesign results
# ------------------------------------------------------------------------------

load(here("output", "res_prior_predictive_re.rds"))

simdesign_meta <- c("REPLICATIONS", "SIM_TIME", "RAM_USED", "SEED", "COMPLETED")
summaries_re <- res |>
  select(-any_of(simdesign_meta))

saveRDS(summaries_re, here("output", "prior_predictive_re_summaries.rds"))
message("Saved: output/prior_predictive_re_summaries.rds  (",
        nrow(summaries_re), " conditions)")


# ------------------------------------------------------------------------------
# Re-simulate draws for visualisation
# ------------------------------------------------------------------------------
# For each of the 40 conditions, generate a large set of subject-level
# (p0, delta_p) draws for density/ridge plots in the report.

set.seed(42)

n_studies  <- 2000   # more studies for smoother marginal densities
n_subjects <- 50

draws_list <- purrr::pmap(
  list(
    hp_family    = summaries_re$hp_family,
    hp_mean      = summaries_re$hp_mean,
    re_structure = summaries_re$re_structure
  ),
  function(hp_family, hp_mean, re_structure) {

    # Fixed-effects priors (recommended values)
    b0 <- rlogis(n_studies, location = 0, scale = 0.75)
    b1 <- rlogis(n_studies, location = 0, scale = 0.25)

    # Study-level hyperprior draws
    sd_re_study    <- sample_sd_hyperprior(n_studies, hp_family, hp_mean)
    sd_slope_study <- if (re_structure == "intercept_slope") {
      sample_sd_hyperprior(n_studies, hp_family, hp_mean)
    } else {
      rep(0, n_studies)
    }

    # Expand to subject level
    b0_exp       <- rep(b0,             each = n_subjects)
    b1_exp       <- rep(b1,             each = n_subjects)
    sd_re_exp    <- rep(sd_re_study,    each = n_subjects)
    sd_slope_exp <- rep(sd_slope_study, each = n_subjects)

    n_total <- n_studies * n_subjects
    u <- rnorm(n_total, 0, sd_re_exp)
    v <- if (re_structure == "intercept_slope") {
      rnorm(n_total, 0, sd_slope_exp)
    } else {
      rep(0, n_total)
    }

    g_inv <- function(x) apply_inverse_link(x, link = "logit")

    data.frame(
      hp_family    = hp_family,
      hp_mean      = hp_mean,
      re_structure = re_structure,
      p_intercept  = g_inv(b0_exp + u),
      delta_p      = g_inv(b0_exp + u + b1_exp + v) -
                     g_inv(b0_exp + u - b1_exp - v),
      sd_re        = sd_re_exp
    )
  }
)

draws_re_long <- bind_rows(draws_list) |>
  mutate(
    abs_delta_p  = abs(delta_p),
    hp_family    = factor(hp_family,
                          levels = c("exponential", "gamma",
                                     "half_normal", "half_t"),
                          labels = c("Exponential", "Gamma(2, .)",
                                     "Half-Normal", "Half-t(3, 0, .)")),
    hp_mean      = factor(hp_mean),
    re_structure = factor(re_structure,
                          levels = c("intercept_only", "intercept_slope"),
                          labels = c("Intercept only", "Intercept + slope"))
  )

saveRDS(draws_re_long, here("output", "prior_predictive_re_draws.rds"))
message("Saved: output/prior_predictive_re_draws.rds  (",
        nrow(draws_re_long), " rows)")
