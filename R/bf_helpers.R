# bf_helpers.R
# Helper functions for Goal 2: BF calibration via brms.
# Each function fits a brms model for one simulated dataset and returns
# the Bayes factor BF10 (H1 vs H0).
#
# Required packages: brms, posterior

fit_and_get_bf <- function(dat,
                           link     = c("logit", "probit"),
                           dist_b0  = c("logistic", "normal"),
                           sd_b0    = 0.75,
                           dist_b1  = c("normal", "logistic", "cauchy"),
                           sd_b1    = 0.25,
                           n_chains = 4,
                           n_iter   = 2000,
                           n_cores  = 1,
                           silent   = TRUE) {
  link    <- match.arg(link)
  dist_b0 <- match.arg(dist_b0)
  dist_b1 <- match.arg(dist_b1)

  family <- switch(link,
                   logit  = brms::binomial(link = "logit"),
                   probit = brms::binomial(link = "probit")
  )

  prior_intercept <- make_brms_prior(dist_b0, sd_b0, class = "Intercept")
  prior_b1        <- make_brms_prior(dist_b1, sd_b1, class = "b", coef = "condition")

  # H1: includes effect of condition
  fit_h1 <- brms::brm(
    formula  = y | trials(n) ~ 1 + condition,
    data     = dat,
    family   = family,
    prior    = c(prior_intercept, prior_b1),
    chains   = n_chains,
    iter     = n_iter,
    cores    = n_cores,
    sample_prior = TRUE,
    save_pars = brms::save_pars(all = TRUE),  # required for bridge sampling
    silent   = if (silent) 2L else 0L,
    refresh  = 0
  )

  # H0: intercept only (condition effect constrained to 0)
  fit_h0 <- brms::brm(
    formula  = y | trials(n) ~ 1,
    data     = dat,
    family   = family,
    prior    = prior_intercept,
    chains   = n_chains,
    iter     = n_iter,
    cores    = n_cores,
    sample_prior = TRUE,
    save_pars = brms::save_pars(all = TRUE),
    silent   = if (silent) 2L else 0L,
    refresh  = 0
  )

  ml_h1 <- suppressMessages(bridgesampling::bridge_sampler(fit_h1, silent = silent))
  ml_h0 <- suppressMessages(bridgesampling::bridge_sampler(fit_h0, silent = silent))
  bf10 <- bridgesampling::bf(ml_h1, ml_h0)$bf

  as.numeric(bf10)
}

fit_and_get_bf_sd <- function(dat,
                              link     = c("logit", "probit"),
                              dist_b0  = c("logistic", "normal"),
                              sd_b0    = 0.75,
                              dist_b1  = c("normal", "logistic", "cauchy"),
                              sd_b1    = 0.25,
                              n_chains = 4,
                              n_iter   = 2000,
                              n_cores  = 1,
                              silent   = TRUE) {
  link    <- match.arg(link)
  dist_b0 <- match.arg(dist_b0)
  dist_b1 <- match.arg(dist_b1)

  family <- switch(link,
                   logit  = brms::binomial(link = "logit"),
                   probit = brms::binomial(link = "probit")
  )

  prior_intercept <- make_brms_prior(dist_b0, sd_b0, class = "Intercept")
  prior_b1        <- make_brms_prior(dist_b1, sd_b1,
                                     class = "b", coef = "condition")

  # Fit H1 only — no save_pars(all) or sample_prior needed for SD ratio
  fit_h1 <- brms::brm(
    formula = y | trials(n) ~ 1 + condition + (1 | gr(subject_id, by = "condition")),
    data    = dat,
    family  = family,
    prior   = c(prior_intercept, prior_b1),
    chains  = n_chains,
    iter    = n_iter,
    cores   = n_cores,
    silent  = if (silent) 2L else 0L,
    refresh = 0,
    backend = "cmdstanr"
  )

  prior_at_0 <- switch(dist_b1,
                       normal   = dnorm(0,   mean = 0,     sd    = sd_b1),
                       logistic = dlogis(0,  location = 0, scale = sd_b1),
                       cauchy   = dcauchy(0, location = 0, scale = sd_b1)
  )

  post_samples <- posterior::as_draws_df(fit_h1)[["b_condition"]]
  x_range  <- range(post_samples)
  x_pad    <- diff(x_range) * 0.1
  kde      <- density(post_samples, bw = "nrd0",
                      from = min(x_range[1] - x_pad, 0),
                      to   = max(x_range[2] + x_pad, 0))
  post_at_0 <- approxfun(kde$x, kde$y, rule = 2)(0)

  if (is.na(post_at_0) || is.nan(post_at_0)) return(NA_real_)
  if (post_at_0 <= 0) return(Inf)

  bf10 <- prior_at_0 / post_at_0
  if (is.na(bf10) || is.nan(bf10)) return(NA_real_)
  as.numeric(bf10)
}

make_brms_prior <- function(dist  = c("normal", "logistic", "cauchy"),
                            scale = 1,
                            class = "b",
                            coef  = NULL) {
  dist <- match.arg(dist)

  # brms prior string: normal(0, s), logistic(0, s), cauchy(0, s)
  prior_str <- switch(dist,
                      normal   = paste0("normal(0, ", scale, ")"),
                      logistic = paste0("logistic(0, ", scale, ")"),
                      cauchy   = paste0("cauchy(0, ", scale, ")")
  )

  if (!is.null(coef)) {
    brms::prior_string(prior_str, class = class, coef = coef)
  } else {
    brms::prior_string(prior_str, class = class)
  }
}

make_brms_sd_prior <- function(sd_prior_re = c("default", "gamma", "exponential"),
                               group       = "subject_id",
                               coef        = "Intercept") {
  sd_prior_re <- match.arg(sd_prior_re)

  prior_str <- switch(sd_prior_re,
                      default     = "student_t(3, 0, 2.5)",
                      gamma       = "gamma(2, 4)",
                      exponential = "exponential(4)"
  )

  brms::prior_string(prior_str, class = "sd", group = group, coef = coef)
}

fit_and_get_bf_sd_re <- function(dat,
                                 link         = c("logit", "probit"),
                                 dist_b0      = c("logistic", "normal"),
                                 sd_b0        = 0.75,
                                 dist_b1      = c("normal", "logistic", "cauchy"),
                                 sd_b1        = 0.25,
                                 sd_prior_re  = c("default", "gamma", "exponential"),
                                 random_slope = FALSE,
                                 n_chains     = 4,
                                 n_iter       = 2000,
                                 n_cores      = 1,
                                 adapt_delta  = 0.95,
                                 silent       = TRUE) {
  link        <- match.arg(link)
  dist_b0     <- match.arg(dist_b0)
  dist_b1     <- match.arg(dist_b1)
  sd_prior_re <- match.arg(sd_prior_re)

  family <- switch(link,
                   logit  = binomial(link = "logit"),
                   probit = binomial(link = "probit")
  )

  prior_intercept <- make_brms_prior(dist_b0, sd_b0, class = "Intercept")
  prior_b1        <- make_brms_prior(dist_b1, sd_b1,
                                     class = "b", coef = "condition")
  prior_sd_int    <- make_brms_sd_prior(sd_prior_re, group = "subject_id",
                                        coef = "Intercept")

  if (random_slope) {
    prior_sd_slope <- make_brms_sd_prior(sd_prior_re, group = "subject_id",
                                         coef = "condition")
    re_formula <- y | trials(n) ~ 1 + condition + (1 + condition | subject_id)
    all_priors <- c(prior_intercept, prior_b1, prior_sd_int, prior_sd_slope)
  } else {
    re_formula <- y | trials(n) ~ 1 + condition + (1 | subject_id)
    all_priors <- c(prior_intercept, prior_b1, prior_sd_int)
  }

  # Fit H1: fixed effects + random effects per subject
  fit_h1 <- brms::brm(
    formula = re_formula,
    data    = dat,
    family  = family,
    prior   = all_priors,
    chains  = n_chains,
    iter    = n_iter,
    cores   = n_cores,
    control = list(adapt_delta = adapt_delta),
    silent  = if (silent) 2L else 0L,
    refresh = 0,
    backend = "cmdstanr"
  )

  prior_at_0 <- switch(dist_b1,
                       normal   = dnorm(0,   mean = 0,     sd    = sd_b1),
                       logistic = dlogis(0,  location = 0, scale = sd_b1),
                       cauchy   = dcauchy(0, location = 0, scale = sd_b1)
  )

  post_samples <- posterior::as_draws_df(fit_h1)[["b_condition"]]
  x_range  <- range(post_samples)
  x_pad    <- diff(x_range) * 0.1
  kde      <- density(post_samples, bw = "nrd0",
                      from = min(x_range[1] - x_pad, 0),
                      to   = max(x_range[2] + x_pad, 0))
  post_at_0 <- approxfun(kde$x, kde$y, rule = 2)(0)

  if (is.na(post_at_0) || is.nan(post_at_0)) return(NA_real_)
  # post_at_0 == 0 means posterior has no mass at 0 — extreme evidence for H1
  if (post_at_0 <= 0) return(Inf)

  bf10 <- prior_at_0 / post_at_0
  if (is.na(bf10) || is.nan(bf10)) return(NA_real_)
  # Inf is a valid result (posterior far from 0); NA/NaN is not
  as.numeric(bf10)
}

fit_and_get_bf_bs_re <- function(dat,
                                 link         = c("logit", "probit"),
                                 dist_b0      = c("logistic", "normal"),
                                 sd_b0        = 0.75,
                                 dist_b1      = c("normal", "logistic", "cauchy"),
                                 sd_b1        = 0.25,
                                 sd_prior_re  = c("default", "gamma", "exponential"),
                                 n_chains     = 4,
                                 n_iter       = 2000,
                                 n_cores      = 1,
                                 adapt_delta  = 0.95,
                                 silent       = TRUE) {
  link        <- match.arg(link)
  dist_b0     <- match.arg(dist_b0)
  dist_b1     <- match.arg(dist_b1)
  sd_prior_re <- match.arg(sd_prior_re)

  family <- switch(link,
                   logit  = binomial(link = "logit"),
                   probit = binomial(link = "probit")
  )

  prior_intercept <- make_brms_prior(dist_b0, sd_b0, class = "Intercept")
  prior_b1        <- make_brms_prior(dist_b1, sd_b1,
                                     class = "b", coef = "condition")
  prior_sd_int    <- make_brms_sd_prior(sd_prior_re, group = "subject_id",
                                        coef = "Intercept")

  fit_h1 <- brms::brm(
    formula      = y | trials(n) ~ 1 + condition + (1 | subject_id),
    data         = dat,
    family       = family,
    prior        = c(prior_intercept, prior_b1, prior_sd_int),
    chains       = n_chains,
    iter         = n_iter,
    cores        = n_cores,
    control      = list(adapt_delta = adapt_delta),
    sample_prior = TRUE,
    save_pars    = brms::save_pars(all = TRUE),
    silent       = if (silent) 2L else 0L,
    refresh      = 0,
    backend      = "cmdstanr"
  )

  fit_h0 <- brms::brm(
    formula      = y | trials(n) ~ 1 + (1 | subject_id),
    data         = dat,
    family       = family,
    prior        = c(prior_intercept, prior_sd_int),
    chains       = n_chains,
    iter         = n_iter,
    cores        = n_cores,
    control      = list(adapt_delta = adapt_delta),
    sample_prior = TRUE,
    save_pars    = brms::save_pars(all = TRUE),
    silent       = if (silent) 2L else 0L,
    refresh      = 0,
    backend      = "cmdstanr"
  )

  # Bridge sampling for marginal likelihoods.
  # Wrap each call: high-dimensional posteriors (many random effects) can
  # produce NA log-densities that crash the bridge sampler's convergence loop.
  ml_h1 <- tryCatch(
    suppressMessages(bridgesampling::bridge_sampler(fit_h1, silent = silent)),
    error = function(e) NULL
  )
  ml_h0 <- tryCatch(
    suppressMessages(bridgesampling::bridge_sampler(fit_h0, silent = silent)),
    error = function(e) NULL
  )

  if (is.null(ml_h1) || is.null(ml_h0)) return(NA_real_)

  bf10 <- tryCatch(
    bridgesampling::bf(ml_h1, ml_h0)$bf,
    error = function(e) NA_real_
  )
  if (!is.finite(bf10)) return(NA_real_)
  as.numeric(bf10)
}

compute_bf_counts <- function(bf10, threshold = 3) {
  valid   <- !is.na(bf10)
  bf_val  <- bf10[valid]
  n_valid <- length(bf_val)

  n_h1           <- sum(bf_val > threshold)
  n_h0           <- sum(bf_val < 1 / threshold)
  n_inconclusive <- n_valid - n_h1 - n_h0

  c(n_valid = as.integer(n_valid),
    n_h1    = as.integer(n_h1),
    n_h0    = as.integer(n_h0),
    n_inconclusive = as.integer(n_inconclusive))
}

beta_posterior_summary <- function(k, n, alpha0 = 3, beta0 = 3) {
  k <- unname(k)
  n <- unname(n)
  post_alpha <- alpha0 + k
  post_beta  <- beta0 + (n - k)

  c(post_alpha  = post_alpha,
    post_beta   = post_beta,
    post_mean   = post_alpha / (post_alpha + post_beta),
    post_median = qbeta(0.5,   post_alpha, post_beta),
    ci_lower    = qbeta(0.025, post_alpha, post_beta),
    ci_upper    = qbeta(0.975, post_alpha, post_beta))
}

dirichlet_posterior_summary <- function(n_h1, n_h0, n_inconclusive,
                                        alpha0 = c(1, 3, 1)) {
  counts     <- c(unname(n_h1), unname(n_inconclusive), unname(n_h0))
  post_alpha <- alpha0 + counts
  alpha_sum  <- sum(post_alpha)

  # Marginal Beta summary for one category
  marginal_summary <- function(a_j, alpha_sum, prefix) {
    b_j <- alpha_sum - a_j
    setNames(
      c(a_j / alpha_sum,
        qbeta(0.5,   a_j, b_j),
        qbeta(0.025, a_j, b_j),
        qbeta(0.975, a_j, b_j)),
      paste0(prefix, c("_mean", "_median", "_ci_lower", "_ci_upper"))
    )
  }

  c(post_alpha_h1  = post_alpha[1],
    post_alpha_inc = post_alpha[2],
    post_alpha_h0  = post_alpha[3],
    marginal_summary(post_alpha[1], alpha_sum, "h1"),
    marginal_summary(post_alpha[2], alpha_sum, "inc"),
    marginal_summary(post_alpha[3], alpha_sum, "h0"))
}
