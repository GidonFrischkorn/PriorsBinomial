# ==============================================================================
# Script 02: Post-process Prior Predictive Simulation Results
# ==============================================================================

library(SimDesign)
library(dplyr)
library(tidyr)
library(here)

source(here("R", "link_functions.R"))
source(here("R", "prior_sampling.R"))


# ------------------------------------------------------------------------------
# Load SimDesign results
# ------------------------------------------------------------------------------

load(here("output", "res_prior_predictive.rds"))

simdesign_meta <- c("REPLICATIONS", "SIM_TIME", "RAM_USED", "SEED", "COMPLETED")
summaries <- res |>
  select(-any_of(simdesign_meta))

saveRDS(summaries, here("output", "prior_predictive_summaries.rds"))
message("Saved: output/prior_predictive_summaries.rds  (", nrow(summaries), " conditions)")


# ------------------------------------------------------------------------------
# Generate long-format draws for density ridge plots
# (Sample a subset of conditions for visualization — not all 720)
# ------------------------------------------------------------------------------

fig1_conditions <- expand.grid(
  link    = c("logit", "probit"),
  dist_b0 = c("logistic", "normal", "cauchy"),
  sd_b0   = c(0.50, 0.75, 1.00, 1.50),
  dist_b1 = "logistic",
  sd_b1   = 0.25,
  stringsAsFactors = FALSE
)

fig2_conditions <- bind_rows(
  expand.grid(
    link    = "logit",
    dist_b0 = "logistic",
    sd_b0   = 0.75,
    dist_b1 = c("normal", "logistic", "cauchy"),
    sd_b1   = c(0.10, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.00, 1.50, 2.00),
    stringsAsFactors = FALSE
  ),
  expand.grid(
    link    = "probit",
    dist_b0 = "normal",
    sd_b0   = 0.75,
    dist_b1 = c("normal", "logistic", "cauchy"),
    sd_b1   = c(0.10, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.00, 1.50, 2.00),
    stringsAsFactors = FALSE
  )
)

viz_conditions <- bind_rows(fig1_conditions, fig2_conditions) |> distinct()

# Re-simulate draws for visualization (10,000 draws per visualization condition)
set.seed(42)
draws_list <- purrr::pmap(viz_conditions, function(link, dist_b0, sd_b0, dist_b1, sd_b1) {
  n <- 10000
  b0 <- sample_prior(n, dist = dist_b0, scale = sd_b0)
  b1 <- sample_prior(n, dist = dist_b1, scale = sd_b1)
  p0 <- apply_inverse_link(b0, link = link)
  p1 <- apply_inverse_link(b0 + b1, link = link)
  p2 <- apply_inverse_link(b0 - b1, link = link)
  data.frame(
    link        = link,
    dist_b0     = dist_b0,
    sd_b0       = sd_b0,
    dist_b1     = dist_b1,
    sd_b1       = sd_b1,
    p_intercept = p0,
    delta_p     = p1 - p2,
    abs_delta_p = abs(p1 - p2)
  )
})

draws_long <- bind_rows(draws_list)
saveRDS(draws_long, here("output", "prior_predictive_draws.rds"))
message("Saved: output/prior_predictive_draws.rds  (", nrow(draws_long), " rows)")
