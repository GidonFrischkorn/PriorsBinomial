# prior_sampling.R
# Functions for drawing samples from prior distributions on the logit/probit
# scale. Supports normal, logistic, and Cauchy distributions.


#' Draw samples from a prior distribution
#'
#' @param n        Integer. Number of draws.
#' @param dist     Character. Prior family: "normal", "logistic", or "cauchy".
#' @param location Numeric. Location parameter (mean for normal, center for
#'                 logistic/cauchy). Default 0.
#' @param scale    Numeric. Scale parameter (SD for normal, scale for
#'                 logistic/cauchy). Default 1.
#'
#' @return Numeric vector of length n.
sample_prior <- function(n,
                          dist     = c("normal", "logistic", "cauchy"),
                          location = 0,
                          scale    = 1) {
  dist <- match.arg(dist)
  switch(dist,
    normal   = rnorm(n,    mean     = location, sd    = scale),
    logistic = rlogis(n,   location = location, scale = scale),
    cauchy   = rcauchy(n,  location = location, scale = scale)
  )
}
