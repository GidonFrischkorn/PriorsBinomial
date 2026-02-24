# bf_helpers.R
# Helper functions for Goal 2: BF calibration via brms.
# Each function fits a brms model for one simulated dataset and returns
# the Bayes factor BF10 (H1 vs H0).
#
# Required packages: brms, posterior


#' Fit H1 and H0 brms models and return BF10 via bridge sampling
#'
#' H1: k ~ Binomial(n, p),  logit(p) = b0 + b1 * condition
#' H0: k ~ Binomial(n, p),  logit(p) = b0
#'
#' Uses bridge sampling via bridgesampling::bridge_sampler() to estimate
#' the marginal likelihood for each model, then returns BF10 = ml_H1 / ml_H0.
#'
#' @param dat         data.frame. Columns: y (successes), n (trials),
#'                    condition (sum-to-zero contrast, values +1 and -1).
#' @param link        Character. "logit" or "probit".
#' @param dist_b0     Character. Intercept prior family: "logistic" or "normal".
#' @param sd_b0       Numeric. Intercept prior scale. Default 0.75.
#' @param dist_b1     Character. Effect prior family: "normal", "logistic", or "cauchy".
#' @param sd_b1       Numeric. Effect prior scale.
#' @param n_chains    Integer. MCMC chains. Default 4.
#' @param n_iter      Integer. MCMC iterations per chain. Default 2000.
#' @param n_cores     Integer. Parallel chains. Default 1 (parallelism handled
#'                    at the simulation level by SimDesign).
#' @param silent      Logical. Suppress brms messages. Default TRUE.
#'
#' @return Numeric scalar. BF10 (Bayes factor for H1 over H0).
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


  # Bridge sampling for marginal likelihoods
  ml_h1 <- bridgesampling::bridge_sampler(fit_h1, silent = silent)
  ml_h0 <- bridgesampling::bridge_sampler(fit_h0, silent = silent)

  bf10 <- bridgesampling::bf(ml_h1, ml_h0)$bf

  as.numeric(bf10)
}


#' Fit H1 brms model and return BF10 via the Savage-Dickey density ratio
#'
#' H1: k ~ Binomial(n, p),  g(p) = b0 + b1 * condition
#' H0: k ~ Binomial(n, p),  g(p) = b0  [implicit — b1 = 0]
#'
#' The Savage-Dickey density ratio applies because H0 is nested in H1 and
#' the priors on b0 are identical in both models (independence of priors).
#' BF10 = pi(b1 = 0) / p(b1 = 0 | data):
#'   BF10 > 1  =>  evidence for H1 (b1 != 0)
#'   BF10 < 1  =>  evidence for H0 (b1 = 0)
#'
#' Prior density is computed analytically (dnorm/dlogis/dcauchy at 0).
#' Posterior density is estimated via KDE (stats::density, Silverman's
#' bandwidth) interpolated with stats::approxfun().
#'
#' @param dat         data.frame. Columns: y (successes), n (trials),
#'                    condition (sum-to-zero contrast, values +1 and -1).
#' @param link        Character. "logit" or "probit".
#' @param dist_b0     Character. Intercept prior: "logistic" or "normal".
#' @param sd_b0       Numeric. Intercept prior scale. Default 0.75.
#' @param dist_b1     Character. Effect prior: "normal", "logistic",
#'                    or "cauchy".
#' @param sd_b1       Numeric. Effect prior scale.
#' @param n_chains    Integer. MCMC chains. Default 4.
#' @param n_iter      Integer. MCMC iterations per chain. Default 2000.
#' @param n_cores     Integer. Parallel chains. Default 1.
#' @param silent      Logical. Suppress brms/Stan messages. Default TRUE.
#'
#' @return Numeric scalar. BF10. Returns NA_real_ on error or when KDE
#'   evaluation yields a non-finite density.
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

  # Analytical prior density at b1 = 0
  prior_at_0 <- switch(dist_b1,
                       normal   = dnorm(0,   mean = 0,     sd    = sd_b1),
                       logistic = dlogis(0,  location = 0, scale = sd_b1),
                       cauchy   = dcauchy(0, location = 0, scale = sd_b1)
  )

  # KDE posterior density at b1 = 0
  # Extend KDE range to always cover 0, guarding against NA from approxfun()
  # when all posterior mass is far from zero (large true-effect conditions).
  post_samples <- posterior::as_draws_df(fit_h1)[["b_condition"]]
  x_range  <- range(post_samples)
  x_pad    <- diff(x_range) * 0.1
  kde      <- density(post_samples, bw = "nrd0",
                      from = min(x_range[1] - x_pad, 0),
                      to   = max(x_range[2] + x_pad, 0))
  post_at_0 <- approxfun(kde$x, kde$y, rule = 2)(0)

  if (!is.finite(post_at_0) || post_at_0 <= 0) return(NA_real_)

  bf10 <- prior_at_0 / post_at_0
  if (!is.finite(bf10)) return(NA_real_)
  as.numeric(bf10)
}


#' Build a brms prior specification for a given distribution family
#'
#' @param dist   Character. Prior family: "normal", "logistic", or "cauchy".
#' @param scale  Numeric. Scale parameter (sd for normal, scale for others).
#' @param class  Character. brms prior class (e.g., "Intercept", "b").
#' @param coef   Character or NULL. Coefficient name for class = "b".
#'
#' @return A brms prior object.
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


#' Build a brms prior specification for the random-effect SD
#'
#' Constructs the brms prior for class = "sd" for a named grouping factor.
#' Three options ranging from very wide (default brms) to strongly constrained
#' priors that concentrate mass below 0.5 on the logit/probit scale.
#'
#' Scale interpretation (logit scale around p = 0.5):
#'   sd_re = 0.25 => ~6% between-subject SD in probability
#'   sd_re = 0.5  => ~12% between-subject SD in probability
#'   student_t(3, 0, 2.5) places 95% mass up to ~7 logit units (implausibly wide)
#'
#' @param sd_prior_re Character. One of:
#'   "default"     — student_t(3, 0, 2.5), the brms default; very wide.
#'   "gamma"       — gamma(2, 4); mean = 0.5, P(sd < 0.5) ~= 0.71.
#'   "exponential" — exponential(4); mean = 0.25, P(sd < 0.5) ~= 0.86.
#' @param group Character. Grouping factor name in the brms formula.
#'   Default "subject_id".
#' @param coef Character. Coefficient name for the SD prior (e.g.,
#'   "Intercept" for random intercepts, "condition" for random slopes).
#'   Default "Intercept".
#'
#' @return A brms prior object for class = "sd", group = group, coef = coef.
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


#' Fit within-subjects H1 brms model and return BF10 via Savage-Dickey ratio
#'
#' Extends fit_and_get_bf_sd() to a within-subjects design with subject-level
#' random intercepts, optionally also including random slopes for condition.
#' Each subject appears in both conditions.
#'
#' The Savage-Dickey density ratio remains valid: H0 (b_condition = 0) is
#' nested in H1, and all other priors — including RE structure — are identical
#' in both models. Only one brms fit (H1) is needed.
#'
#' H1 (intercept only):  y | trials(n) ~ 1 + condition + (1 | subject_id)
#' H1 (intercept+slope): y | trials(n) ~ 1 + condition + (1 + condition | subject_id)
#' H0: implicit — b_condition = 0 (Savage-Dickey density ratio)
#'
#' @param dat          data.frame. Columns: y (successes), n (trials),
#'                     condition (sum-to-zero contrast, values +1 and -1),
#'                     subject_id (integer; each value appears exactly twice,
#'                     once per condition).
#' @param link         Character. "logit" or "probit".
#' @param dist_b0      Character. Intercept prior: "logistic" or "normal".
#' @param sd_b0        Numeric. Intercept prior scale. Default 0.75.
#' @param dist_b1      Character. Effect prior: "normal", "logistic",
#'                     or "cauchy".
#' @param sd_b1        Numeric. Effect prior scale.
#' @param sd_prior_re  Character. SD prior label applied to both the random
#'                     intercept SD and (when random_slope = TRUE) the random
#'                     slope SD: "default", "gamma", or "exponential".
#'                     See make_brms_sd_prior() for details.
#' @param random_slope Logical. If TRUE, also include a by-subject random slope
#'                     for condition: (1 + condition | subject_id). An LKJ(1)
#'                     (uniform) prior is used for the intercept-slope
#'                     correlation. Default FALSE.
#' @param n_chains     Integer. MCMC chains. Default 4.
#' @param n_iter       Integer. MCMC iterations per chain. Default 2500
#'                     (slightly higher than fixed-effects version for RE
#'                     model stability).
#' @param n_cores      Integer. Parallel chains. Default 1.
#' @param silent       Logical. Suppress brms/Stan messages. Default TRUE.
#'
#' @return Numeric scalar. BF10. Returns NA_real_ on error or when KDE
#'   evaluation yields a non-finite density.
fit_and_get_bf_sd_re <- function(dat,
                                 link         = c("logit", "probit"),
                                 dist_b0      = c("logistic", "normal"),
                                 sd_b0        = 0.75,
                                 dist_b1      = c("normal", "logistic", "cauchy"),
                                 sd_b1        = 0.25,
                                 sd_prior_re  = c("default", "gamma", "exponential"),
                                 random_slope = FALSE,
                                 n_chains     = 4,
                                 n_iter       = 2500,
                                 n_cores      = 1,
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
    silent  = if (silent) 2L else 0L,
    refresh = 0,
    backend = "cmdstanr"
  )

  # Analytical prior density at b_condition = 0
  prior_at_0 <- switch(dist_b1,
                       normal   = dnorm(0,   mean = 0,     sd    = sd_b1),
                       logistic = dlogis(0,  location = 0, scale = sd_b1),
                       cauchy   = dcauchy(0, location = 0, scale = sd_b1)
  )

  # KDE posterior density at b_condition = 0
  post_samples <- posterior::as_draws_df(fit_h1)[["b_condition"]]
  x_range  <- range(post_samples)
  x_pad    <- diff(x_range) * 0.1
  kde      <- density(post_samples, bw = "nrd0",
                      from = min(x_range[1] - x_pad, 0),
                      to   = max(x_range[2] + x_pad, 0))
  post_at_0 <- approxfun(kde$x, kde$y, rule = 2)(0)

  if (!is.finite(post_at_0) || post_at_0 <= 0) return(NA_real_)

  bf10 <- prior_at_0 / post_at_0
  if (!is.finite(bf10)) return(NA_real_)
  as.numeric(bf10)
}
