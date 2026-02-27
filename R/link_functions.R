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

#' Logit transform (log odds)
#'
#' @param p Numeric vector of probabilities in (0, 1).
#' @return Numeric vector on the real line.
logit <- function(p) log(p / (1 - p))

#' Inverse logit (logistic function)
#'
#' @param x Numeric vector on the real line.
#' @return Numeric vector of probabilities in (0, 1).
ilogit <- function(x) 1 / (1 + exp(-x))

#' Probit transform (inverse normal CDF)
#'
#' @param p Numeric vector of probabilities in (0, 1).
#' @return Numeric vector on the real line.
probit <- function(p) qnorm(p)

#' Inverse probit (normal CDF)
#'
#' @param x Numeric vector on the real line.
#' @return Numeric vector of probabilities in (0, 1).
iprobit <- function(x) pnorm(x)

#' Apply link function by name (probability → linear predictor)
#'
#' @param p    Numeric vector of probabilities in (0, 1).
#' @param link Character. One of "logit" or "probit".
#' @return Numeric vector on the real line.
apply_link <- function(p, link = c("logit", "probit")) {
  link <- match.arg(link)
  switch(link,
    logit  = logit(p),
    probit = probit(p)
  )
}

#' Apply inverse link function by name
#'
#' @param x    Numeric vector on the real line.
#' @param link Character. One of "logit" or "probit".
#' @return Numeric vector of probabilities in (0, 1).
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

#' Exact density of p = ilogit(b0) when b0 ~ Logistic(0, scale)
#'
#' Derived by change of variables: f_p(p) = f_b(logit(p)) / (p * (1 - p))
#' where f_b = dlogis(x, 0, scale).
#'
#' Special case: at scale = 1, this equals dunif(p, 0, 1) = 1 for all p in (0,1).
#' This is the exact analytical expression of the matched prior principle for
#' the logit link.
#'
#' @param p     Numeric vector of probabilities in (0, 1).
#' @param scale Positive numeric. Scale parameter of the Logistic prior on b0.
#'              Default 1 gives Uniform(0, 1) on the probability scale.
#' @return Numeric vector of density values.
d_logistic_on_p <- function(p, scale = 1) {
  dlogis(logit(p), location = 0, scale = scale) / (p * (1 - p))
}

#' Exact density of p = iprobit(b0) when b0 ~ Normal(0, sigma)
#'
#' Derived by change of variables: f_p(p) = f_b(qnorm(p)) / dnorm(qnorm(p))
#' where f_b = dnorm(x, 0, sigma).
#'
#' Simplifies to: f_p(p) = (1 / sigma) * phi(Phi^{-1}(p) / sigma) / phi(Phi^{-1}(p))
#' where phi = standard normal density, Phi^{-1} = qnorm.
#'
#' Special case: at sigma = 1, this equals dunif(p, 0, 1) = 1 for all p in (0,1).
#' This is the exact analytical expression of the matched prior principle for
#' the probit link.
#'
#' @param p     Numeric vector of probabilities in (0, 1).
#' @param sigma Positive numeric. Standard deviation of the Normal prior on b0.
#'              Default 1 gives Uniform(0, 1) on the probability scale.
#' @return Numeric vector of density values.
d_normal_on_p <- function(p, sigma = 1) {
  z <- qnorm(p)
  dnorm(z, mean = 0, sd = sigma) / dnorm(z, mean = 0, sd = 1)
}


# ==============================================================================
# Validation helpers
# ==============================================================================

#' Validate matched prior analytical densities via KS test
#'
#' Checks that simulated p = g^{-1}(b0) values, when b0 is drawn from the
#' matched prior with scale/sigma = 1, are consistent with Uniform(0, 1).
#'
#' @param n_draws Integer. Number of Monte Carlo draws. Default 100000.
#' @param seed    Integer or NULL. Random seed.
#' @return List with KS test results for logit and probit cases.
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

#' Construct orthonormal contrasts for a J-level factor
#'
#' Returns a J x (J-1) orthonormal contrast matrix Q and the predictor SD
#' of each column (which equals 1/sqrt(J) for balanced designs).
#'
#' The construction follows Rouder et al. (2012, p. 363): Helmert contrasts
#' are normalised to unit length. The specific basis does not matter for BF
#' computation when a spherically symmetric prior is used.
#'
#' @param J Integer >= 2. Number of factor levels.
#' @return List with components:
#'   \item{Q}{J x (J-1) orthonormal contrast matrix}
#'   \item{predictor_sd}{Numeric scalar, SD of each contrast column = 1/sqrt(J)}
#'   \item{scaling_factor}{Numeric scalar, sd_b1 multiplier = sqrt(J)}
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
