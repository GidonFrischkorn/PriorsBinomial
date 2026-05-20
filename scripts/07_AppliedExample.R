# ==============================================================================
# Script 07: Applied Example — Oberauer (2019) Simple Span
# ==============================================================================
# Dataset: Oberauer, K. (2019). Working memory capacity limits memory for
#   bindings. Journal of Cognition, 2(1), 40. https://doi.org/10.5334/joc.86
#
# Data location: tutorial-m3-bmm/data/Oberauer_2019_SimpleSpan_agg.csv
#   (bmm tutorial project: https://osf.io/yb7wm/)
#
# Analysis: Binomial GLM with set size (ss_lin) as a continuous predictor of
# recall accuracy (corr / n_total), using the matched Logistic(0, 0.75)
# intercept prior and Logistic(0, 0.25) effect prior.
#
# Three steps of the workflow:
#   Step 1+2: Prior predictive check (matched vs. misfit Cauchy)
#   Step 3:   BF via Savage-Dickey density ratio + sensitivity analysis
# ==============================================================================

library(here)
library(dplyr)
library(ggplot2)
library(brms)
library(posterior)

source(here("R", "link_functions.R"))
source(here("R", "prior_sampling.R"))
source(here("R", "plotting.R"))

# ------------------------------------------------------------------------------
# Load and prepare data
# ------------------------------------------------------------------------------

data_path <- here("data", "Oberauer_2019_SimpleSpan_agg.csv")

d <- read.csv(data_path)
d_exp1 <- subset(d, exp == "closedset")
d_exp1$n_total  <- with(d_exp1, corr + other + npl)
# Standardise ss_lin to SD = 1 so Logistic(0, 0.25) applies directly.
# Raw ss_lin has SD ≈ 2.08; without scaling the prior would cover ~5×
# the intended effect range relative to the ±1 framework in the paper.
d_exp1$ss_lin_z <- d_exp1$ss_lin / sd(d_exp1$ss_lin)

message("Loaded Exp1 (closedset): ",
        length(unique(d_exp1$id)), " participants, ",
        nrow(d_exp1), " rows, set sizes: ",
        paste(sort(unique(d_exp1$setsize)), collapse = ", "))

# ------------------------------------------------------------------------------
# Step 1 + 2: Prior predictive check
# ------------------------------------------------------------------------------

set.seed(2024)
n_draws <- 20000

# Intercept: matched Logistic(0, 0.75)
b0_matched  <- rlogis(n_draws, location = 0, scale = 0.75)
p0_matched  <- ilogit(b0_matched)

# Effect: Logistic(0, 0.25) vs Cauchy(0, 2.5).
# Evaluate at ±1 SD of set size, matching the standardised predictor.
b1_matched <- rlogis(n_draws, location = 0, scale = 0.25)
b1_cauchy  <- rcauchy(n_draws, location = 0, scale = 2.5)

delta_p_matched <- ilogit(b0_matched + b1_matched) -
  ilogit(b0_matched - b1_matched)
delta_p_cauchy  <- ilogit(b0_matched + b1_cauchy) -
  ilogit(b0_matched - b1_cauchy)

pp_data <- bind_rows(
  data.frame(delta_p = delta_p_matched, prior = "Matched Logistic(0, 0.25)"),
  data.frame(delta_p = delta_p_cauchy,  prior = "Cauchy(0, 2.5)")
)

p_pp <- ggplot(pp_data, aes(x = delta_p, fill = prior, colour = prior)) +
  geom_density(alpha = 0.4, linewidth = 0.7) +
  geom_vline(xintercept = c(-0.20, -0.10, 0.10, 0.20),
             linetype = "dotted", colour = "grey40") +
  scale_x_continuous(expression(delta * italic(p)), limits = c(-1, 1)) +
  scale_fill_viridis_d("Effect prior", option = "D", end = 0.85) +
  scale_colour_viridis_d("Effect prior", option = "D", end = 0.85) +
  ggtitle("Prior predictive: matched vs. Cauchy effect prior") +
  prior_theme()

ggsave(here("figures", "applied_example_pp_check.pdf"),
       p_pp, width = 8, height = 5)
message("Saved: figures/applied_example_pp_check.pdf")

# ------------------------------------------------------------------------------
# Step 3: Fit H1 and H0, compute BF via Savage-Dickey density ratio
# ------------------------------------------------------------------------------

out_h1 <- here("output", "applied_example_fit_h1.rds")
out_h0 <- here("output", "applied_example_fit_h0.rds")

n_chains <- 4
n_iter   <- 4000

if (!file.exists(out_h1)) {
  fit_h1 <- brm(
    formula = corr | trials(n_total) ~ 1 + ss_lin_z + (1 | id),
    family  = binomial(link = "logit"),
    prior   = c(
      prior(logistic(0, 0.75), class = Intercept),
      prior(logistic(0, 0.25), class = b),
      prior(exponential(4),    class = sd)
    ),
    data         = d_exp1,
    chains       = n_chains,
    iter         = n_iter,
    sample_prior = TRUE,
    backend      = "cmdstanr",
    silent       = 2L,
    refresh      = 0
  )
  saveRDS(fit_h1, out_h1)
  message("Saved: ", out_h1)
} else {
  fit_h1 <- readRDS(out_h1)
  message("Loaded H1 fit from: ", out_h1)
}

if (!file.exists(out_h0)) {
  fit_h0 <- brm(
    formula = corr | trials(n_total) ~ 1 + (1 | id),
    family  = binomial(link = "logit"),
    prior   = c(
      prior(logistic(0, 0.75), class = Intercept),
      prior(exponential(4),    class = sd)
    ),
    data         = d_exp1,
    chains       = n_chains,
    iter         = n_iter,
    sample_prior = TRUE,
    backend      = "cmdstanr",
    silent       = 2L,
    refresh      = 0
  )
  saveRDS(fit_h0, out_h0)
  message("Saved: ", out_h0)
} else {
  fit_h0 <- readRDS(out_h0)
  message("Loaded H0 fit from: ", out_h0)
}

# ------------------------------------------------------------------------------
# Bridge sampling BF (robust when posterior is far from null)
# ------------------------------------------------------------------------------
# Savage-Dickey fails here: posterior sits ~30 SDs from 0, so KDE density at 0
# is numerically 0 regardless of chain length. Bridge sampling integrates the
# marginal likelihood directly and is unaffected by this.
# Requires save_pars(all = TRUE); separate fits without sample_prior.

library(bridgesampling)

out_h1_bs <- here("output", "applied_example_fit_h1_bs.rds")
out_h0_bs <- here("output", "applied_example_fit_h0_bs.rds")

if (!file.exists(out_h1_bs)) {
  fit_h1_bs <- brm(
    formula   = corr | trials(n_total) ~ 1 + ss_lin_z + (1 | id),
    family    = binomial(link = "logit"),
    prior     = c(
      prior(logistic(0, 0.75), class = Intercept),
      prior(logistic(0, 0.25), class = b),
      prior(exponential(4),    class = sd)
    ),
    data      = d_exp1,
    chains    = n_chains,
    iter      = n_iter,
    save_pars = save_pars(all = TRUE),
    backend   = "cmdstanr",
    silent    = 2L,
    refresh   = 0
  )
  saveRDS(fit_h1_bs, out_h1_bs)
  message("Saved: ", out_h1_bs)
} else {
  fit_h1_bs <- readRDS(out_h1_bs)
  message("Loaded H1-BS fit from: ", out_h1_bs)
}

if (!file.exists(out_h0_bs)) {
  fit_h0_bs <- brm(
    formula   = corr | trials(n_total) ~ 1 + (1 | id),
    family    = binomial(link = "logit"),
    prior     = c(
      prior(logistic(0, 0.75), class = Intercept),
      prior(exponential(4),    class = sd)
    ),
    data      = d_exp1,
    chains    = n_chains,
    iter      = n_iter,
    save_pars = save_pars(all = TRUE),
    backend   = "cmdstanr",
    silent    = 2L,
    refresh   = 0
  )
  saveRDS(fit_h0_bs, out_h0_bs)
  message("Saved: ", out_h0_bs)
} else {
  fit_h0_bs <- readRDS(out_h0_bs)
  message("Loaded H0-BS fit from: ", out_h0_bs)
}

lml_h1 <- bridge_sampler(fit_h1_bs, silent = TRUE)
lml_h0 <- bridge_sampler(fit_h0_bs, silent = TRUE)
bf_bs  <- bayes_factor(lml_h1, lml_h0)

message(sprintf("Bridge sampling BF10:          %.3e", bf_bs$bf))
message(sprintf("Bridge sampling log10(BF10):   %.2f", log10(bf_bs$bf)))

bs_result <- data.frame(
  method     = "bridge_sampling",
  BF10       = bf_bs$bf,
  log10_BF10 = log10(bf_bs$bf)
)
saveRDS(bs_result, here("output", "applied_example_bf_bridge.rds"))
message("Saved: output/applied_example_bf_bridge.rds")

# Savage-Dickey density ratio
sd_b1        <- 0.25
prior_at_0   <- dlogis(0, location = 0, scale = sd_b1)
post_samples <- as_draws_df(fit_h1)[["b_ss_lin_z"]]
x_range      <- range(post_samples)
x_pad        <- diff(x_range) * 0.1
kde          <- density(post_samples, bw = "nrd0",
                        from = min(x_range[1] - x_pad, 0),
                        to   = max(x_range[2] + x_pad, 0))
post_at_0    <- approxfun(kde$x, kde$y, rule = 2)(0)
bf10         <- prior_at_0 / post_at_0

message(sprintf("BF10 (sd_b1 = 0.25): %.2f", bf10))
message(sprintf("log10(BF10): %.2f", log10(bf10)))

# ------------------------------------------------------------------------------
# Savage-Dickey visualization: prior vs. posterior for b_ss_lin
# ------------------------------------------------------------------------------
# The posterior spike at -0.59 vs the prior at 0 makes the BF intuitive:
# the BF ratio is the prior-to-posterior density ratio at the null (x = 0).

x_seq      <- seq(-1.3, 0.9, length.out = 1000)
prior_line <- dlogis(x_seq, location = 0, scale = sd_b1)

kde_b1      <- density(post_samples, bw = "nrd0", n = 2048,
                       from = min(post_samples) - 0.1,
                       to   = max(post_samples) + 0.1)
post_line   <- approxfun(kde_b1$x, kde_b1$y, rule = 1)(x_seq)
post_line[is.na(post_line)] <- 0   # zero density outside the posterior range

sd_df <- data.frame(
  x            = rep(x_seq, 2),
  density      = c(prior_line, post_line),
  distribution = rep(c("Prior: Logistic(0, 0.25)", "Posterior"), each = 1000)
)

p_sd <- ggplot(sd_df, aes(x = x, y = density,
                          colour = distribution, linetype = distribution)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 0.05, y = max(prior_line) * 0.6, hjust = 0,
           label = "Null value", colour = "grey30", size = 4) +
  scale_x_continuous(expression(beta[1] ~ "(set size effect)")) +
  scale_y_continuous("Density") +
  scale_colour_viridis_d("", option = "D", end = 0.85) +
  scale_linetype_manual("", values = c("Prior: Logistic(0, 0.25)" = "solid",
                                       "Posterior" = "dashed")) +
  ggtitle("Savage-Dickey: prior vs. posterior for set size effect") +
  prior_theme()

ggsave(here("figures", "applied_example_savage_dickey.pdf"),
       p_sd, width = 8, height = 5)
message("Saved: figures/applied_example_savage_dickey.pdf")

# ------------------------------------------------------------------------------
# Sensitivity analysis over sd_b1
# ------------------------------------------------------------------------------

sensitivity_sds <- c(0.15, 0.25, 0.35)
bf_sensitivity  <- numeric(length(sensitivity_sds))

for (j in seq_along(sensitivity_sds)) {
  s          <- sensitivity_sds[j]
  p_at_0     <- dlogis(0, location = 0, scale = s)
  bf_sensitivity[j] <- p_at_0 / post_at_0
  message(sprintf("BF10 (sd_b1 = %.2f): %.3e  [log10: %.2f]",
                  s, bf_sensitivity[j], log10(bf_sensitivity[j])))
}

sensitivity_results <- data.frame(
  sd_b1 = sensitivity_sds,
  BF10  = bf_sensitivity,
  log10_BF10 = log10(bf_sensitivity)
)

saveRDS(sensitivity_results,
        here("output", "applied_example_sensitivity.rds"))
message("Saved: output/applied_example_sensitivity.rds")

# ------------------------------------------------------------------------------
# Posterior summary for b_ss_lin
# ------------------------------------------------------------------------------

post_summary <- summarise_draws(
  subset_draws(as_draws_df(fit_h1), variable = "b_ss_lin_z"),
  mean, median,
  ~quantile(.x, c(0.025, 0.975)),
  default_convergence_measures()
)

print(post_summary)

# Implied delta_p at each set size
# ss_lin values in data: -3.83, -1.83, 0.17, 2.17 (approx)
b0_post  <- mean(as_draws_df(fit_h1)[["b_Intercept"]])
b1_post  <- mean(as_draws_df(fit_h1)[["b_ss_lin_z"]])
ss_vals  <- sort(unique(d_exp1$ss_lin_z))

implied_p <- data.frame(
  setsize   = sort(unique(d_exp1$setsize)),
  ss_lin_z  = ss_vals,
  p_hat     = ilogit(b0_post + b1_post * ss_vals)
)

message("Implied recall probabilities at posterior means:")
print(implied_p)

# ------------------------------------------------------------------------------
# Fitted vs. observed: population-level predicted probability by set size
# ------------------------------------------------------------------------------

b0_draws <- as_draws_df(fit_h1)[["b_Intercept"]]
b1_draws <- as_draws_df(fit_h1)[["b_ss_lin_z"]]

# Predicted probability from fixed effects only (population-level)
pred_mat <- sapply(ss_vals, function(ss) ilogit(b0_draws + b1_draws * ss))

pred_df <- data.frame(
  setsize  = sort(unique(d_exp1$setsize)),
  ss_lin_z = ss_vals,
  p_mean  = colMeans(pred_mat),
  p_lo    = apply(pred_mat, 2, quantile, 0.025),
  p_hi    = apply(pred_mat, 2, quantile, 0.975)
)

d_exp1$p_corr <- d_exp1$corr / d_exp1$n_total

p_fit <- ggplot(pred_df, aes(x = setsize)) +
  geom_ribbon(aes(ymin = p_lo, ymax = p_hi), fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = p_mean), colour = "steelblue", linewidth = 0.9) +
  geom_jitter(data = d_exp1, aes(y = p_corr),
              width = 0.12, height = 0,
              alpha = 0.5, size = 1.8, colour = "grey30") +
  scale_x_continuous("Set size", breaks = sort(unique(d_exp1$setsize))) +
  scale_y_continuous("Recall accuracy (proportion correct)",
                     limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  ggtitle("Fitted vs. observed: binomial GLM with matched logistic priors") +
  prior_theme()

ggsave(here("figures", "applied_example_fitted_vs_obs.pdf"),
       p_fit, width = 7, height = 5)
message("Saved: figures/applied_example_fitted_vs_obs.pdf")

saveRDS(pred_df, here("output", "applied_example_pred_df.rds"))
message("Saved: output/applied_example_pred_df.rds")
