# link_functions.R
# Forward and inverse link functions for binomial GLMs, plus exact analytical
# density formulas for the "matched prior" principle.
#
# The matched prior principle (probability integral transform):
#   If X ~ D with CDF F, then F(X) ~ Uniform(0, 1).
#   - Logit link: g^{-1} = ilogit = logistic CDF
#     => Logistic(0, 1) intercept prior => Uniform(0, 1) on p [exact]
#   - Probit link: g^{-1} = pnorm = normal CDF
#     => Normal(0, 1) intercept prior => Uniform(0, 1) on p [exact]


# ==============================================================================
# Link functions
# ==============================================================================

logit <- function(p) log(p / (1 - p))
ilogit <- function(x) 1 / (1 + exp(-x))

probit <- function(p) qnorm(p)
iprobit <- function(x) pnorm(x)

apply_link <- function(p, link = c("logit", "probit")) {
  link <- match.arg(link)
  switch(link,
    logit  = logit(p),
    probit = probit(p)
  )
}

apply_inverse_link <- function(x, link = c("logit", "probit")) {
  link <- match.arg(link)
  switch(link,
    logit  = ilogit(x),
    probit = iprobit(x)
  )
}


# ==============================================================================
# Analytical density formulas (matched prior cases)
# ==============================================================================

d_logistic_on_p <- function(p, scale = 1) {
  dlogis(logit(p), location = 0, scale = scale) / (p * (1 - p))
}

d_normal_on_p <- function(p, sigma = 1) {
  z <- qnorm(p)
  dnorm(z, mean = 0, sd = sigma) / dnorm(z, mean = 0, sd = 1)
}


# ==============================================================================
# Validation helpers
# ==============================================================================

validate_matched_prior <- function(n_draws = 100000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Logit + Logistic(0, 1) => Uniform(0, 1)
  b0_logit <- rlogis(n_draws, location = 0, scale = 1)
  p_logit  <- ilogit(b0_logit)
  ks_logit <- ks.test(p_logit, "punif")

  # Probit + Normal(0, 1) => Uniform(0, 1)
  b0_probit <- rnorm(n_draws, mean = 0, sd = 1)
  p_probit  <- iprobit(b0_probit)
  ks_probit <- ks.test(p_probit, "punif")

  list(
    logit_ks  = ks_logit,
    probit_ks = ks_probit
  )
}


# ==============================================================================
# Orthonormal contrast helpers
# ==============================================================================

orthonormal_contrasts <- function(J) {
  stopifnot(J >= 2L)
  # Helmert contrasts, then normalise each column to unit length
  H <- contr.helmert(J)
  Q <- apply(H, 2, function(col) col / sqrt(sum(col^2)))

  # Predictor SD: sqrt(mean(q^2)) = sqrt(1/J) since sum(q^2) = 1 and mean(q) = 0
  pred_sd <- sqrt(1 / J)

  list(
    Q              = Q,
    predictor_sd   = pred_sd,
    scaling_factor = sqrt(J)
  )
}


# ==============================================================================
# Conditional transforms: prior on delta_p at a fixed baseline (Reviewer 3)
# ==============================================================================
# Model: g(p) = eta0 + b1 * x, x in {-1, +1}, with eta0 = b0 + u held fixed.
# delta_p(b1) = g^{-1}(eta0 + b1) - g^{-1}(eta0 - b1) is odd and strictly
# increasing in b1, with range (-1, 1). Densities are obtained by change of
# variables (Jacobian), CDFs and quantiles by monotone mapping. No integration
# over b0 or u. The effect prior is family(location = 0, scale = scale); for
# the Logistic family `scale` is the logistic scale parameter, not an SD.

#' Density of the inverse link (derivative of g^{-1})
#'
#' @param x Numeric vector on the linear-predictor scale.
#' @param link Character, "logit" or "probit".
#' @return Numeric vector: d g^{-1}(x) / dx.
d_inverse_link <- function(x, link = c("logit", "probit")) {
  link <- match.arg(link)
  switch(link,
    logit  = dlogis(x),
    probit = dnorm(x)
  )
}

#' Density, CDF and quantile function of a zero-centred effect prior
#'
#' @param x Numeric vector (density/CDF: value of b1; quantile: probability).
#' @param family Character, "normal", "logistic" or "cauchy".
#' @param scale Positive scale parameter (sd for normal; scale otherwise).
#' @param type Character, "d" (density), "p" (CDF) or "q" (quantile).
#' @return Numeric vector of the requested function values.
effect_prior_fun <- function(x,
                             family = c("normal", "logistic", "cauchy"),
                             scale  = 1,
                             type   = c("d", "p", "q")) {
  family <- match.arg(family)
  type   <- match.arg(type)
  fun <- switch(family,
    normal   = switch(type, d = dnorm,  p = pnorm,  q = qnorm),
    logistic = switch(type, d = dlogis, p = plogis, q = qlogis),
    cauchy   = switch(type, d = dcauchy, p = pcauchy, q = qcauchy)
  )
  fun(x, 0, scale)
}

#' delta_p as a function of b1 at a fixed baseline linear predictor
#'
#' @param b1 Numeric vector of effect values (linear-predictor scale, +-1 coding).
#' @param eta0 Numeric scalar, baseline linear predictor b0 + u.
#' @param link Character, "logit" or "probit".
#' @return Numeric vector g^{-1}(eta0 + b1) - g^{-1}(eta0 - b1).
delta_p_from_b1 <- function(b1, eta0, link = c("logit", "probit")) {
  link <- match.arg(link)
  apply_inverse_link(eta0 + b1, link) - apply_inverse_link(eta0 - b1, link)
}

#' Inverse of delta_p_from_b1: the b1 that produces a given delta_p
#'
#' @param d Numeric vector of delta_p values in (-1, 1).
#' @param eta0 Numeric scalar, baseline linear predictor b0 + u.
#' @param link Character, "logit" or "probit".
#' @param tol Numeric tolerance passed to uniroot().
#' @return Numeric vector of b1 values (NaN outside (-1, 1); +-Inf at +-1).
b1_from_delta_p <- function(d, eta0, link = c("logit", "probit"),
                            tol = 1e-12) {
  link <- match.arg(link)
  vapply(d, function(di) {
    if (is.na(di) || abs(di) > 1) return(NaN)
    if (di == 0) return(0)
    if (abs(di) == 1) return(sign(di) * Inf)
    # odd function: solve for |d| and restore the sign
    a <- abs(di)
    f <- function(b) delta_p_from_b1(b, eta0, link) - a
    root <- uniroot(f, lower = 0, upper = 1, extendInt = "upX",
                    tol = tol, maxiter = 1000)$root
    sign(di) * root
  }, numeric(1))
}

#' Density of delta_p at a fixed baseline (change of variables)
#'
#' @param d Numeric vector of delta_p values in (-1, 1).
#' @param eta0 Numeric scalar, baseline linear predictor b0 + u
#'   (equivalently g(p0) for a baseline probability p0).
#' @param link Character, "logit" or "probit".
#' @param family Character, effect prior family: "normal", "logistic", "cauchy".
#' @param scale Positive scale parameter of the effect prior.
#' @return Numeric vector of densities pi_b1(b1(d)) / |d delta_p / d b1|;
#'   0 for |d| >= 1 (outside the support), NA for NA input.
#' @details The density can diverge at d -> +-1: the Jacobian falls off with
#'   the tails of the inverse-link density while a heavier-tailed prior falls
#'   off more slowly. Checked numerically at eta0 = 0, scales 0.25 / 0.707 /
#'   1.5 (logs/t3b_ddeltap_check.log): divergent for probit + Cauchy and
#'   probit + Logistic at all three scales, probit + Normal at scale 1.5,
#'   logit + Cauchy at all three scales, logit + Logistic at scale 1.5;
#'   vanishing in the other cases. The singularity is integrable (the density
#'   is a change of variables of a proper prior; total mass = 1 via
#'   `p_delta_p()`), but much of the mass can sit extremely close to +-1
#'   (probit + Cauchy(0, 0.707): 6.3% within 1e-12 of +-1), so numerical
#'   integration of this density underestimates it. Use `p_delta_p()` for
#'   probability masses.
d_delta_p <- function(d, eta0, link = c("logit", "probit"),
                      family = c("normal", "logistic", "cauchy"),
                      scale = 1) {
  link   <- match.arg(link)
  family <- match.arg(family)
  b1 <- b1_from_delta_p(d, eta0, link)
  jac <- d_inverse_link(eta0 + b1, link) + d_inverse_link(eta0 - b1, link)
  num <- effect_prior_fun(b1, family, scale, "d")
  out <- num / jac
  # 0/0 (prior density and Jacobian both underflow): the prior term is 0
  out[is.nan(out) & !is.na(num) & num == 0] <- 0
  # delta_p has support (-1, 1): density 0 at and beyond +-1; NA input stays NA
  out[!is.na(d) & (abs(d) >= 1 | !is.finite(b1))] <- 0
  out
}

#' CDF of delta_p at a fixed baseline
#'
#' @param d Numeric vector of delta_p values.
#' @param eta0 Numeric scalar, baseline linear predictor b0 + u.
#' @param link Character, "logit" or "probit".
#' @param family Character, effect prior family.
#' @param scale Positive scale parameter of the effect prior.
#' @return Numeric vector P(delta_p <= d) = F_b1(b1(d)).
p_delta_p <- function(d, eta0, link = c("logit", "probit"),
                      family = c("normal", "logistic", "cauchy"),
                      scale = 1) {
  link   <- match.arg(link)
  family <- match.arg(family)
  dd <- pmin(pmax(d, -1), 1)
  effect_prior_fun(b1_from_delta_p(dd, eta0, link), family, scale, "p")
}

#' Quantile function of delta_p at a fixed baseline
#'
#' @param prob Numeric vector of probabilities in [0, 1].
#' @param eta0 Numeric scalar, baseline linear predictor b0 + u.
#' @param link Character, "logit" or "probit".
#' @param family Character, effect prior family.
#' @param scale Positive scale parameter of the effect prior.
#' @return Numeric vector delta_p(Q_b1(prob)) (exact, monotone mapping).
q_delta_p <- function(prob, eta0, link = c("logit", "probit"),
                      family = c("normal", "logistic", "cauchy"),
                      scale = 1) {
  link   <- match.arg(link)
  family <- match.arg(family)
  delta_p_from_b1(effect_prior_fun(prob, family, scale, "q"), eta0, link)
}

#' CDF of |delta_p| at a fixed baseline
#'
#' @param a Numeric vector of non-negative |delta_p| values.
#' @param eta0 Numeric scalar, baseline linear predictor b0 + u.
#' @param link Character, "logit" or "probit".
#' @param family Character, effect prior family.
#' @param scale Positive scale parameter of the effect prior.
#' @return Numeric vector P(|delta_p| <= a) = 2 F_b1(b1(a)) - 1.
p_abs_delta_p <- function(a, eta0, link = c("logit", "probit"),
                          family = c("normal", "logistic", "cauchy"),
                          scale = 1) {
  link   <- match.arg(link)
  family <- match.arg(family)
  aa <- pmin(pmax(a, 0), 1)
  2 * effect_prior_fun(b1_from_delta_p(aa, eta0, link), family, scale, "p") - 1
}

#' Quantile function of |delta_p| at a fixed baseline
#'
#' @param prob Numeric vector of probabilities in [0, 1].
#' @param eta0 Numeric scalar, baseline linear predictor b0 + u.
#' @param link Character, "logit" or "probit".
#' @param family Character, effect prior family.
#' @param scale Positive scale parameter of the effect prior.
#' @return Numeric vector delta_p(Q_b1((1 + prob) / 2)); exact because
#'   |delta_p| = delta_p(|b1|) and the prior is symmetric about 0.
q_abs_delta_p <- function(prob, eta0, link = c("logit", "probit"),
                          family = c("normal", "logistic", "cauchy"),
                          scale = 1) {
  link   <- match.arg(link)
  family <- match.arg(family)
  delta_p_from_b1(effect_prior_fun((1 + prob) / 2, family, scale, "q"),
                  eta0, link)
}

#' Scale giving the same IQR of b1 as a Logistic(0, s) prior
#'
#' IQR of Logistic(0, s) = 2 s log(3); Normal(0, sigma) = 2 sigma qnorm(0.75);
#' Cauchy(0, gamma) = 2 gamma.
#'
#' @param s_logistic Numeric vector of logistic scale parameters.
#' @param family Character, target family: "normal", "logistic", "cauchy".
#' @return Numeric vector of target-family scales with matched IQR.
matched_iqr_scale <- function(s_logistic,
                              family = c("normal", "logistic", "cauchy")) {
  family <- match.arg(family)
  switch(family,
    normal   = s_logistic * log(3) / qnorm(0.75),
    logistic = s_logistic,
    cauchy   = s_logistic * log(3)
  )
}

#' Exact floor/ceiling mass of an intercept prior on the probability scale
#'
#' @param scale Numeric vector of positive scale parameters of the intercept prior.
#' @param family Character, "normal", "logistic" or "cauchy".
#' @param link Character, "logit" or "probit".
#' @param lower Numeric, lower probability bound (e.g. 0.05).
#' @param upper Numeric, upper probability bound (e.g. 0.95).
#' @return Numeric vector P(p < lower) + P(p > upper) = F(g(lower)) + 1 - F(g(upper)).
floor_ceiling_mass <- function(scale, family = c("normal", "logistic", "cauchy"),
                               link = c("logit", "probit"),
                               lower = 0.05, upper = 0.95) {
  family <- match.arg(family)
  link   <- match.arg(link)
  vapply(scale, function(s) {
    effect_prior_fun(apply_link(lower, link), family, s, "p") +
      1 - effect_prior_fun(apply_link(upper, link), family, s, "p")
  }, numeric(1))
}
