# prior_sampling.R
# Functions for drawing samples from prior distributions on the logit/probit
# scale. Supports normal, logistic, and Cauchy distributions.

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
