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
