# ==============================================================================
# Script 02b: Post-process RE Prior Predictive Simulation Results
# ==============================================================================
#
# PURPOSE:
#   Load per-condition results from SimDesign (script 01b), combine into a
#   clean tidy data frame, and re-simulate subject-level draws for
#   visualisation in the Quarto report.
#
# INPUT:  output/res_prior_predictive_re.rds
# OUTPUT: output/prior_predictive_re_summaries.rds  (tidy summaries, 6 rows)
#         output/prior_predictive_re_draws.rds       (long-format subject draws)
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

sample_sd_hyperprior <- function(n, sd_prior_re) {
  switch(sd_prior_re,
    student_t   = abs(rt(n, df = 3)) * 2.5,
    gamma       = rgamma(n, shape = 2, rate = 4),
    exponential = rexp(n, rate = 4)
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
# For each of the 6 conditions, generate a large set of subject-level
# (p0, delta_p) draws for density/ridge plots in the report.

set.seed(42)

n_studies  <- 2000   # more studies for smoother marginal densities
n_subjects <- 50

draws_list <- purrr::pmap(
  list(
    sd_prior_re  = summaries_re$sd_prior_re,
    re_structure = summaries_re$re_structure
  ),
  function(sd_prior_re, re_structure) {

    # Fixed-effects priors (recommended values)
    b0 <- rlogis(n_studies, location = 0, scale = 0.75)
    b1 <- rnorm(n_studies,  mean = 0,     sd = 0.25)

    # Study-level hyperprior draws
    sd_re_study    <- sample_sd_hyperprior(n_studies, sd_prior_re)
    sd_slope_study <- if (re_structure == "intercept_slope") {
      sample_sd_hyperprior(n_studies, sd_prior_re)
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
      sd_prior_re  = sd_prior_re,
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
    sd_prior_re  = factor(sd_prior_re,
                           levels = c("student_t", "gamma", "exponential"),
                           labels = c("student_t(3, 0, 2.5)",
                                      "gamma(2, 4)",
                                      "exponential(4)")),
    re_structure = factor(re_structure,
                           levels = c("intercept_only", "intercept_slope"),
                           labels = c("Intercept only", "Intercept + slope"))
  )

saveRDS(draws_re_long, here("output", "prior_predictive_re_draws.rds"))
message("Saved: output/prior_predictive_re_draws.rds  (",
        nrow(draws_re_long), " rows)")
