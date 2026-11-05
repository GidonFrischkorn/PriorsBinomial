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

# ==============================================================================
# Random-slopes BF calibration with one precompiled model per arm (Script 13)
# ==============================================================================
# The effect prior enters the Stan model as data: a user-defined
# effect_prior_lpdf() switches on an integer family code and takes the scale
# as a data value. One compiled model therefore serves every
# family x scale combination; refits go through update(recompile = FALSE).
# The second parameter is a SCALE, not an SD (Logistic SD = scale * pi / sqrt(3)).

#' Integer code of an effect-prior family
#'
#' @param dist_b1 Character, one of "normal", "logistic", "cauchy".
#' @return Integer 1 (normal), 2 (logistic) or 3 (cauchy), as used by
#'   `effect_prior_lpdf()` in the Stan functions block.
effect_prior_code <- function(dist_b1 = c("normal", "logistic", "cauchy")) {
  dist_b1 <- match.arg(dist_b1)
  switch(dist_b1, normal = 1L, logistic = 2L, cauchy = 3L)
}

#' Stanvars that pass the effect prior (family, scale) to Stan as data
#'
#' @param dist_b1 Character, effect-prior family ("normal", "logistic", "cauchy").
#' @param scale_b1 Positive numeric, scale of the effect prior (location 0).
#' @return A `brms` stanvars object: the `effect_prior_lpdf()` function plus the
#'   data variables `prior_fam_b1` and `prior_scale_b1`.
make_effect_prior_stanvars <- function(dist_b1 = c("normal", "logistic", "cauchy"),
                                       scale_b1 = 0.25) {
  stopifnot(is.numeric(scale_b1), length(scale_b1) == 1, scale_b1 > 0)
  fun_code <- "
  real effect_prior_lpdf(real x, int fam, real s) {
    if (fam == 1) return normal_lpdf(x | 0, s);
    if (fam == 2) return logistic_lpdf(x | 0, s);
    return cauchy_lpdf(x | 0, s);
  }
"
  brms::stanvar(scode = fun_code, block = "functions") +
    brms::stanvar(effect_prior_code(dist_b1), name = "prior_fam_b1",
                  scode = "int<lower=1, upper=3> prior_fam_b1;") +
    brms::stanvar(as.numeric(scale_b1), name = "prior_scale_b1",
                  scode = "real<lower=0> prior_scale_b1;")
}

#' Prior density of the effect at zero
#'
#' @param dist_b1 Character, effect-prior family.
#' @param scale_b1 Positive numeric, prior scale.
#' @return Numeric density of the zero-centred prior at 0.
prior_density_at_0 <- function(dist_b1 = c("normal", "logistic", "cauchy"),
                               scale_b1) {
  dist_b1 <- match.arg(dist_b1)
  switch(dist_b1,
         normal   = dnorm(0, mean = 0, sd = scale_b1),
         logistic = dlogis(0, location = 0, scale = scale_b1),
         cauchy   = dcauchy(0, location = 0, scale = scale_b1))
}

#' Posterior density of the effect at zero from MCMC draws
#'
#' @param draws Numeric vector of posterior draws of the effect.
#' @param method Character: "kde" (Gaussian KDE, bw = "nrd0", linear
#'   interpolation; identical to `fit_and_get_bf_sd_re()`), "logspline"
#'   (`logspline::logspline()` fit, evaluated with `dlogspline()`), or
#'   "normal" (normal approximation with the draws' mean and SD).
#' @return Numeric density at 0; `NA_real_` if the estimator fails.
posterior_density_at_0 <- function(draws, method = c("kde", "logspline", "normal")) {
  method <- match.arg(method)
  out <- tryCatch(
    switch(method,
           kde = {
             x_range <- range(draws)
             x_pad   <- diff(x_range) * 0.1
             kde     <- density(draws, bw = "nrd0",
                                from = min(x_range[1] - x_pad, 0),
                                to   = max(x_range[2] + x_pad, 0))
             approxfun(kde$x, kde$y, rule = 2)(0)
           },
           logspline = {
             ls_fit <- suppressWarnings(logspline::logspline(draws))
             logspline::dlogspline(0, ls_fit)
           },
           normal = dnorm(0, mean = mean(draws), sd = sd(draws))),
    error = function(e) NA_real_
  )
  as.numeric(out)
}

#' Build (compile) a random-slopes template model for one simulation arm
#'
#' Fits with `chains = 0`, so only compilation happens. The returned object is
#' passed to `fit_bf_sd_rs()` / `fit_bf_bs_rs()`, which refit it on new data
#' with `update(recompile = FALSE)`.
#'
#' @param dat Data frame with columns `y`, `n`, `subject_id` and the predictor.
#' @param link Character, "logit" or "probit".
#' @param dist_b0 Character, intercept-prior family ("logistic" or "normal").
#' @param scale_b0 Positive numeric, intercept-prior scale.
#' @param predictor Character, name of the predictor column ("condition" for
#'   the +/-1 within-subject factor, "x" for a continuous predictor).
#' @param hypothesis Character, "h1" (fixed effect of the predictor, prior set
#'   through stanvars) or "h0" (no fixed effect; random slope kept).
#' @param sd_prior_re Character, prior on both random-effect SDs, passed to
#'   `make_brms_sd_prior()`.
#' @param backend Character, "cmdstanr" (fast, Savage-Dickey) or "rstan"
#'   (needed for bridge sampling without per-fit recompilation).
#' @param random_slope Logical; `TRUE` (default, Script 13) adds a random
#'   slope of the predictor, `(1 + predictor | subject_id)`; `FALSE` fits
#'   random intercepts only, `(1 | subject_id)` (Script 16, between-subjects
#'   design with one row per subject).
#' @return A `brmsfit` without draws (compiled model only).
compile_bf_template_rs <- function(dat,
                                   link        = c("logit", "probit"),
                                   dist_b0     = c("logistic", "normal"),
                                   scale_b0    = 0.75,
                                   predictor   = "condition",
                                   hypothesis  = c("h1", "h0"),
                                   sd_prior_re = c("exponential", "gamma", "default"),
                                   backend     = c("cmdstanr", "rstan"),
                                   random_slope = TRUE) {
  link        <- match.arg(link)
  dist_b0     <- match.arg(dist_b0)
  hypothesis  <- match.arg(hypothesis)
  sd_prior_re <- match.arg(sd_prior_re)
  backend     <- match.arg(backend)

  re_part <- if (random_slope) paste0("(1 + ", predictor, " | subject_id)") else "(1 | subject_id)"
  rhs     <- if (hypothesis == "h1") paste("1 +", predictor, "+", re_part) else paste("1 +", re_part)
  form    <- brms::bf(stats::as.formula(paste("y | trials(n) ~", rhs)))

  priors <- c(make_brms_prior(dist_b0, scale_b0, class = "Intercept"),
              make_brms_sd_prior(sd_prior_re, group = "subject_id", coef = "Intercept"))
  if (random_slope) {
    priors <- c(priors, make_brms_sd_prior(sd_prior_re, group = "subject_id", coef = predictor))
  }
  stanvars <- NULL
  if (hypothesis == "h1") {
    priors <- c(priors,
                brms::prior_string("effect_prior(prior_fam_b1, prior_scale_b1)",
                                   class = "b", coef = predictor))
    stanvars <- make_effect_prior_stanvars("normal", 0.25)
  }

  suppressWarnings(suppressMessages(brms::brm(
    formula   = form,
    data      = dat,
    family    = binomial(link = link),
    prior     = priors,
    stanvars  = stanvars,
    chains    = 0,
    save_pars = brms::save_pars(all = TRUE),
    backend   = backend,
    silent    = 2L,
    refresh   = 0
  )))
}

#' Refit a precompiled template on new data
#'
#' @param template A `brmsfit` template from `compile_bf_template_rs()`.
#' @param dat Data frame for this replication.
#' @param stanvars Optional stanvars (effect prior) for an H1 template.
#' @param seed Integer seed passed to Stan.
#' @param n_chains,n_warmup,n_sampling Integers, MCMC settings (per chain).
#' @param adapt_delta Numeric, NUTS target acceptance rate.
#' @return A `brmsfit` with draws.
refit_template_rs <- function(template, dat, stanvars = NULL, seed,
                              n_chains = 4, n_warmup = 1000, n_sampling = 1000,
                              adapt_delta = 0.95) {
  args <- list(template, newdata = dat, recompile = FALSE, chains = n_chains,
               cores = 1, warmup = n_warmup, iter = n_warmup + n_sampling,
               control = list(adapt_delta = adapt_delta), seed = seed,
               silent = 2L, refresh = 0)
  if (!is.null(stanvars)) args$stanvars <- stanvars
  suppressWarnings(suppressMessages(do.call(stats::update, args)))
}

#' Convergence diagnostics of a random-slopes fit
#'
#' @param fit A `brmsfit` with draws.
#' @param predictor Character, predictor name.
#' @return A list: `b_draws` (numeric draws of the effect), and named numeric
#'   `diag` with `n_div`, `rhat_max`, `ess_bulk_min` and `ess_tail_min` over
#'   the effect and the random-effect SDs present in the fit (both for a
#'   random-slopes model, the intercept SD only for random intercepts).
fit_diagnostics_rs <- function(fit, predictor = "condition") {
  pars  <- c(paste0("b_", predictor), "sd_subject_id__Intercept",
             paste0("sd_subject_id__", predictor))
  pars  <- intersect(pars, brms::variables(fit))
  dr    <- posterior::subset_draws(posterior::as_draws_array(fit), variable = pars)
  mats  <- lapply(pars, function(p) posterior::extract_variable_matrix(dr, p))
  np    <- brms::nuts_params(fit)
  list(b_draws = as.numeric(mats[[1]]),
       diag = c(n_div        = sum(np$Value[np$Parameter == "divergent__"]),
                rhat_max     = max(vapply(mats, posterior::rhat, numeric(1))),
                ess_bulk_min = min(vapply(mats, posterior::ess_bulk, numeric(1))),
                ess_tail_min = min(vapply(mats, posterior::ess_tail, numeric(1)))))
}

#' Savage-Dickey BF10 for one effect prior, refitting a precompiled template
#'
#' @param template A `brmsfit` from `compile_bf_template_rs(hypothesis = "h1")`.
#' @param dat Data frame for this replication.
#' @param dist_b1 Character, effect-prior family.
#' @param scale_b1 Positive numeric, effect-prior scale.
#' @param predictor Character, predictor name used in the template.
#' @param bf_methods Character vector of posterior-density estimators passed to
#'   `posterior_density_at_0()`.
#' @param seed Integer seed for Stan.
#' @param n_chains,n_warmup,n_sampling Integers, MCMC settings (per chain).
#' @param adapt_delta Numeric, NUTS target acceptance rate.
#' @return Named numeric vector: `bf10_<method>` for each method, posterior
#'   mean and SD of the effect, `n_div`, `rhat_max`, `ess_bulk_min`,
#'   `ess_tail_min` (effect and both random-effect SDs) and `wall_s`.
fit_bf_sd_rs <- function(template, dat,
                         dist_b1    = c("normal", "logistic", "cauchy"),
                         scale_b1   = 0.25,
                         predictor  = "condition",
                         bf_methods = c("kde", "logspline", "normal"),
                         seed       = sample.int(.Machine$integer.max, 1),
                         n_chains   = 4,
                         n_warmup   = 1000,
                         n_sampling = 1000,
                         adapt_delta = 0.95) {
  dist_b1 <- match.arg(dist_b1)
  t_start <- Sys.time()
  fit <- refit_template_rs(template, dat,
                           make_effect_prior_stanvars(dist_b1, scale_b1),
                           seed = seed, n_chains = n_chains, n_warmup = n_warmup,
                           n_sampling = n_sampling, adapt_delta = adapt_delta)
  wall_s <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  dg <- fit_diagnostics_rs(fit, predictor)

  prior_0 <- prior_density_at_0(dist_b1, scale_b1)
  bfs <- vapply(bf_methods, function(m) {
    post_0 <- posterior_density_at_0(dg$b_draws, m)
    if (is.na(post_0)) return(NA_real_)
    if (post_0 <= 0) return(Inf)
    prior_0 / post_0
  }, numeric(1))
  names(bfs) <- paste0("bf10_", bf_methods)

  c(bfs, post_mean = mean(dg$b_draws), post_sd = sd(dg$b_draws), dg$diag,
    wall_s = wall_s)
}

#' Bridge-sampled marginal likelihood of H0 (no fixed effect, random slope kept)
#'
#' @param template_h0 rstan `brmsfit` template from
#'   `compile_bf_template_rs(hypothesis = "h0", backend = "rstan")`.
#' @param dat Data frame for this replication.
#' @param seed Integer seed, used for Stan and (offset by 1) for the bridge
#'   sampler's random draws.
#' @param ... MCMC settings passed to `refit_template_rs()`.
#' @return A list with `ml_h0` (`bridge` object) and `wall_h0` (seconds).
#'   Errors (fit or bridge sampling) are not caught here.
fit_h0_bridge_rs <- function(template_h0, dat, seed, ...) {
  t_start <- Sys.time()
  fit_h0 <- refit_template_rs(template_h0, dat, seed = seed, ...)
  ml_h0  <- withr::with_seed(seed + 1L, suppressWarnings(suppressMessages(
    bridgesampling::bridge_sampler(fit_h0, silent = TRUE))))
  list(ml_h0 = ml_h0,
       wall_h0 = as.numeric(difftime(Sys.time(), t_start, units = "secs")))
}

#' Bridge-sampling BF10 for one effect prior with random slopes
#'
#' Refits the precompiled rstan H1 template (no recompilation for
#' `bridge_sampler()`) and divides by a supplied H0 marginal likelihood from
#' `fit_h0_bridge_rs()`. A logspline Savage-Dickey BF10 is computed from the
#' same H1 draws as a cross-check at no extra sampling cost. Errors in the H1
#' fit or bridge sampling are not caught here; `fit_prior_grid_rs()` records
#' them.
#'
#' @param template_h1 rstan `brmsfit` H1 template.
#' @param dat Data frame for this replication.
#' @param dist_b1 Character, effect-prior family.
#' @param scale_b1 Positive numeric, effect-prior scale.
#' @param ml_h0 `bridge` object for H0 on `dat`.
#' @param predictor Character, predictor name used in the template.
#' @param seed Integer seed, used for Stan and (offset by 1) for the bridge
#'   sampler's random draws.
#' @param n_chains,n_warmup,n_sampling Integers, MCMC settings (per chain).
#' @param adapt_delta Numeric, NUTS target acceptance rate.
#' @return Named numeric vector: `bf10_bridge`, `bf10_logspline`,
#'   `pct_error_h1` (100 x bridge-sampling coefficient of variation of the H1
#'   marginal likelihood), posterior mean/SD of the effect, `n_div`,
#'   `rhat_max`, `ess_bulk_min`, `ess_tail_min`, and `wall_s` (H1 refit +
#'   bridge sampling, seconds).
fit_bf_bs_rs <- function(template_h1, dat,
                         dist_b1    = c("normal", "logistic", "cauchy"),
                         scale_b1   = 0.25,
                         ml_h0,
                         predictor  = "condition",
                         seed       = sample.int(.Machine$integer.max - 1L, 1),
                         n_chains   = 4,
                         n_warmup   = 1000,
                         n_sampling = 5000,
                         adapt_delta = 0.95) {
  dist_b1 <- match.arg(dist_b1)
  t_start <- Sys.time()
  fit_h1 <- refit_template_rs(template_h1, dat,
                              make_effect_prior_stanvars(dist_b1, scale_b1),
                              seed = seed, n_chains = n_chains,
                              n_warmup = n_warmup, n_sampling = n_sampling,
                              adapt_delta = adapt_delta)
  ml_h1 <- withr::with_seed(seed + 1L, suppressWarnings(suppressMessages(
    bridgesampling::bridge_sampler(fit_h1, silent = TRUE))))
  bf10 <- exp(ml_h1$logml - ml_h0$logml)
  if (is.na(bf10)) stop("bridge sampling returned NA logml")
  pct_err <- tryCatch(100 * bridgesampling::error_measures(ml_h1)$cv,
                      error = function(e) NA_real_)
  dg <- fit_diagnostics_rs(fit_h1, predictor)
  post_0 <- posterior_density_at_0(dg$b_draws, "logspline")
  bf10_ls <- if (is.na(post_0)) NA_real_ else if (post_0 <= 0) Inf else
    prior_density_at_0(dist_b1, scale_b1) / post_0

  c(bf10_bridge = bf10, bf10_logspline = bf10_ls, pct_error_h1 = pct_err,
    post_mean = mean(dg$b_draws), post_sd = sd(dg$b_draws), dg$diag,
    wall_s = as.numeric(difftime(Sys.time(), t_start, units = "secs")))
}

#' Output names of the per-prior fit functions (used to fill failed fits)
#'
#' @param engine Character, "sd" or "bridge".
#' @param bf_methods Character vector of density estimators (engine "sd").
#' @param ... Ignored (absorbs other fit arguments).
#' @return Character vector of names.
fit_bf_rs_names <- function(engine = c("sd", "bridge"),
                            bf_methods = c("kde", "logspline", "normal"), ...) {
  engine <- match.arg(engine)
  bf <- if (engine == "sd") paste0("bf10_", bf_methods) else
    c("bf10_bridge", "bf10_logspline", "pct_error_h1")
  c(bf, "post_mean", "post_sd", "n_div", "rhat_max", "ess_bulk_min",
    "ess_tail_min", "wall_s")
}

#' Fit the full effect-prior grid to one dataset (paired design)
#'
#' Failures are recorded, not swallowed: each failed prior yields `NA` values
#' (never a capped BF), the number of failures is returned as `n_fail`, and the
#' failure messages are attached as attribute `fail_msgs` and raised as one
#' warning (SimDesign logs it with its message). If H0 bridge sampling fails,
#' or more than half of the priors fail, the replication is stopped with an
#' error so SimDesign logs it and redraws the dataset.
#'
#' @param template A cmdstanr `brmsfit` H1 template (engine "sd").
#' @param dat Data frame for this replication.
#' @param prior_grid Data frame with columns `dist_b1` and `scale_b1`.
#' @param predictor Character, predictor name used in the templates.
#' @param engine Character, "sd" (Savage-Dickey via `fit_bf_sd_rs()`) or
#'   "bridge" (`fit_bf_bs_rs()`, with H0 fitted once per dataset).
#' @param template_h1_bs,template_h0_bs rstan templates for engine "bridge".
#' @param seeds Integer vector of length `nrow(prior_grid) + 1`: one seed per
#'   prior, the last one for H0 (drawn in Analyse with `sample.int()`).
#' @param ... Further arguments passed to the per-prior fit function.
#' @return Named numeric vector: per-prior outputs suffixed with
#'   `_<dist_b1>_<scale_b1>`, plus `n_fail` and (engine "bridge") `wall_s_h0`.
#'   Attribute `fail_msgs`: character vector of failure messages.
fit_prior_grid_rs <- function(template, dat, prior_grid, predictor = "condition",
                              engine = c("sd", "bridge"),
                              template_h1_bs = NULL, template_h0_bs = NULL,
                              seeds = sample.int(.Machine$integer.max - 1L,
                                                 nrow(prior_grid) + 1),
                              ...) {
  engine <- match.arg(engine)
  n_pri  <- nrow(prior_grid)
  stopifnot(length(seeds) == n_pri + 1)
  dots <- list(...)

  ml_h0 <- NULL
  wall_h0 <- NA_real_
  if (engine == "bridge") {
    mcmc <- dots[intersect(names(dots), c("n_chains", "n_warmup", "n_sampling", "adapt_delta"))]
    h0 <- tryCatch(do.call(fit_h0_bridge_rs,
                           c(list(template_h0_bs, dat, seed = seeds[n_pri + 1]), mcmc)),
                   error = function(e) e)
    if (inherits(h0, "error")) stop("H0 fit/bridge sampling failed: ", conditionMessage(h0))
    ml_h0   <- h0$ml_h0
    wall_h0 <- h0$wall_h0
  }

  out_names <- fit_bf_rs_names(engine, ...)
  fail_msgs <- character(0)
  out <- vector("list", n_pri)
  for (i in seq_len(n_pri)) {
    dist  <- as.character(prior_grid$dist_b1[i])
    scale <- prior_grid$scale_b1[i]
    res <- tryCatch(
      if (engine == "sd") {
        fit_bf_sd_rs(template, dat, dist_b1 = dist, scale_b1 = scale,
                     predictor = predictor, seed = seeds[i], ...)
      } else {
        fit_bf_bs_rs(template_h1_bs, dat, dist_b1 = dist, scale_b1 = scale,
                     ml_h0 = ml_h0, predictor = predictor, seed = seeds[i], ...)
      },
      error = function(e) e)
    if (inherits(res, "error")) {
      fail_msgs <- c(fail_msgs, paste0(dist, "(", scale, "): ", conditionMessage(res)))
      res <- setNames(rep(NA_real_, length(out_names)), out_names)
    }
    out[[i]] <- setNames(res, paste0(names(res), "_", dist, "_", format(scale, nsmall = 2)))
  }

  n_fail <- length(fail_msgs)
  if (n_fail > n_pri / 2) {
    stop(sprintf("%d of %d prior fits failed; first: %s", n_fail, n_pri, fail_msgs[1]))
  }
  if (n_fail > 0) {
    warning(sprintf("%d of %d prior fits failed; first: %s", n_fail, n_pri, fail_msgs[1]),
            call. = FALSE)
  }
  res <- c(unlist(out), n_fail = n_fail)
  if (engine == "bridge") res <- c(res, wall_s_h0 = wall_h0)
  attr(res, "fail_msgs") <- fail_msgs
  res
}

# ------------------------------------------------------------------------------
# Design, prior grids and run configuration for Script 13 / 14 (shared so the
# merge step checks against the same Design)
# ------------------------------------------------------------------------------

#' Data-generating Design for the random-slopes calibration (Script 13)
#'
#' Approved grid (2026-10-05): logit_factor fully crossed (true_b1 x
#' n_subjects x n_trials {20, 50} x slope SD {0.10, 0.30, 0.50}, 54 rows);
#' probit_factor and logit_continuous with n_trials = 20 and slope SD 0.30
#' (9 rows each). Row order is fixed; `design_row` is the stable row ID used
#' for seeds, chunking and merging.
#'
#' @param full_crossing Logical; `TRUE` returns all arms crossed with all
#'   factor levels (used to pick test cells), with `design_row = NA`.
#' @return A SimDesign design (tibble) with columns `arm`, `true_b1`,
#'   `n_subjects`, `n_trials`, `true_sd_slope`, `design_row`.
make_design_rs <- function(full_crossing = FALSE) {
  d <- SimDesign::createDesign(
    arm           = c("logit_factor", "probit_factor", "logit_continuous"),
    true_b1       = c(0.00, 0.10, 0.20),
    n_subjects    = c(20, 30, 50),
    n_trials      = c(20, 50),
    true_sd_slope = c(0.10, 0.30, 0.50)
  )
  if (full_crossing) {
    d$design_row <- NA_integer_
    return(d)
  }
  d <- d[d$arm == "logit_factor" | (d$n_trials == 20 & d$true_sd_slope == 0.30), ]
  d$design_row <- seq_len(nrow(d))
  d
}

#' Effect-prior grids per arm (Script 13)
#'
#' @return Named list of data frames (`dist_b1`, `scale_b1`): logit_factor
#'   3 families x 7 scales; probit_factor Normal x 7 scales; logit_continuous
#'   3 families x {0.15, 0.25, 0.30, 0.50}.
make_prior_grids_rs <- function() {
  families    <- c("normal", "logistic", "cauchy")
  scales_full <- c(0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.707)
  list(
    logit_factor     = expand.grid(dist_b1 = families, scale_b1 = scales_full,
                                   stringsAsFactors = FALSE),
    probit_factor    = expand.grid(dist_b1 = "normal", scale_b1 = scales_full,
                                   stringsAsFactors = FALSE),
    logit_continuous = expand.grid(dist_b1 = families,
                                   scale_b1 = c(0.15, 0.25, 0.30, 0.50),
                                   stringsAsFactors = FALSE)
  )
}

#' Parse a row specification such as "1-12", "3,5,7-9" or "all"
#'
#' @param spec Character scalar; empty or "all" selects all rows.
#' @param n_rows Integer, number of Design rows.
#' @return Sorted unique integer vector of row IDs within 1..n_rows.
parse_rows_spec <- function(spec, n_rows) {
  if (is.null(spec) || !nzchar(spec) || identical(spec, "all")) return(seq_len(n_rows))
  parts <- strsplit(gsub("\\s", "", spec), ",")[[1]]
  rows <- unlist(lapply(parts, function(p) {
    if (grepl("^[0-9]+-[0-9]+$", p)) {
      ab <- as.integer(strsplit(p, "-")[[1]])
      seq(ab[1], ab[2])
    } else if (grepl("^[0-9]+$", p)) {
      as.integer(p)
    } else {
      stop("Cannot parse row specification part: '", p, "'")
    }
  }))
  rows <- sort(unique(rows))
  if (any(rows < 1 | rows > n_rows)) {
    stop("Row specification '", spec, "' outside 1..", n_rows)
  }
  rows
}

#' Rows of chunk k out of n_chunks (contiguous, near-equal sizes)
#'
#' @param k Integer chunk index (1-based).
#' @param n_chunks Integer number of chunks.
#' @param n_rows Integer number of Design rows.
#' @return Integer vector of row IDs.
chunk_rows <- function(k, n_chunks, n_rows) {
  stopifnot(k >= 1, k <= n_chunks, n_chunks <= n_rows)
  split(seq_len(n_rows), cut(seq_len(n_rows), n_chunks, labels = FALSE))[[k]]
}

#' Number of worker cores from the environment
#'
#' @return Integer: `SLURM_CPUS_PER_TASK` if set, else `PB_NCORES`, else
#'   `parallel::detectCores() - 2` (at least 1).
resolve_ncores <- function() {
  for (v in c("SLURM_CPUS_PER_TASK", "PB_NCORES")) {
    x <- suppressWarnings(as.integer(Sys.getenv(v, "")))
    if (!is.na(x) && x >= 1) return(x)
  }
  max(1L, parallel::detectCores() - 2L)
}
