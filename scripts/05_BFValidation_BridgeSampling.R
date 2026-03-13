# ==============================================================================
# Script 05: Validate Savage-Dickey BFs Against Bridge Sampling
# ==============================================================================

library(SimDesign)
library(brms)
library(bridgesampling)
library(here)

source(here("R", "link_functions.R"))
source(here("R", "bf_helpers.R"))

# ------------------------------------------------------------------------------
# Design grid (targeted subset)
# ------------------------------------------------------------------------------

Design <- createDesign(
  dist_b1    = c("normal", "logistic", "cauchy"),
  sd_b1      = c(0.25, 0.50),
  true_b1    = c(0.00, 0.10, 0.50)
)

# ------------------------------------------------------------------------------
# SimDesign functions
# ------------------------------------------------------------------------------

#' Generate: identical to Script 03 but with fixed n_subjects and n_trials
Generate <- function(condition, fixed_objects = NULL) {
  Attach(condition)

  n_subjects <- fixed_objects$n_subjects
  n_trials   <- fixed_objects$n_trials
  true_sd_re <- fixed_objects$true_sd_re

  true_p0 <- runif(1, fixed_objects$b0_range[1], fixed_objects$b0_range[2])
  true_b0 <- apply_link(true_p0, "logit")

  g_inv <- function(x) apply_inverse_link(x, link = "logit")

  n_total <- 2L * n_subjects
  cond    <- rep(c(1L, -1L), each = n_subjects)
  u       <- rnorm(n_total, mean = 0, sd = true_sd_re)

  p_subj <- g_inv(true_b0 + true_b1 * cond + u)

  data.frame(
    subject_id = seq_len(n_total),
    y          = rbinom(n_total, n_trials, p_subj),
    n          = n_trials,
    condition  = cond,
    true_p0    = true_p0
  )
}

#' Analyse: compute BF10 via both Savage-Dickey and bridge sampling
#'
#' Returns both BF estimates for the same dataset, enabling a paired comparison.
Analyse <- function(condition, dat, fixed_objects = NULL) {
  Attach(condition)
  sd_b0       <- fixed_objects$sd_b0
  sd_prior_re <- fixed_objects$sd_prior_re
  dist_b0     <- fixed_objects$dist_b0
  true_p0     <- dat$true_p0[1]

  n_cores <- fixed_objects$n_cores

  # Savage-Dickey density ratio (same as Script 03)
  sd_err <- NULL
  bf10_sd <- tryCatch(
    fit_and_get_bf_sd_re(
      dat         = dat,
      link        = "logit",
      dist_b0     = dist_b0,
      sd_b0       = sd_b0,
      dist_b1     = dist_b1,
      sd_b1       = sd_b1,
      sd_prior_re = sd_prior_re,
      n_cores     = n_cores
    ),
    error = function(e) {
      sd_err <<- conditionMessage(e)
      NA_real_
    }
  )
  gc()

  # Bridge sampling (fits H1 + H0 separately)
  bs_err <- NULL
  bf10_bs <- tryCatch(
    fit_and_get_bf_bs_re(
      dat         = dat,
      link        = "logit",
      dist_b0     = dist_b0,
      sd_b0       = sd_b0,
      dist_b1     = dist_b1,
      sd_b1       = sd_b1,
      sd_prior_re = sd_prior_re,
      n_cores     = n_cores
    ),
    error = function(e) {
      bs_err <<- conditionMessage(e)
      NA_real_
    }
  )
  gc()

  # If both methods failed, re-throw so SimDesign logs the actual message
  if (is.na(bf10_sd) && is.na(bf10_bs)) {
    stop("SD: ", sd_err %||% "returned NA", " | BS: ", bs_err %||% "returned NA")
  }

  c(BF10_SD = bf10_sd, BF10_BS = bf10_bs, true_p0 = true_p0)
}

#' Summarise: compare SD and BS estimates
Summarise <- function(condition, results, fixed_objects = NULL) {
  bf_sd <- results[, "BF10_SD"]
  bf_bs <- results[, "BF10_BS"]

  # Paired comparison on log10 scale (more stable)
  both_valid <- !is.na(bf_sd) & !is.na(bf_bs) & bf_sd > 0 & bf_bs > 0
  log_sd     <- log10(bf_sd[both_valid])
  log_bs     <- log10(bf_bs[both_valid])

  # Correlation and agreement
  r_log10 <- if (length(log_sd) > 2) cor(log_sd, log_bs) else NA_real_
  mae_log10 <- mean(abs(log_sd - log_bs))
  bias_log10 <- mean(log_sd - log_bs)

  # Detection-rate agreement at BF > 3 threshold
  detect_sd <- bf_sd[both_valid] > 3
  detect_bs <- bf_bs[both_valid] > 3
  agree_h1  <- mean(detect_sd == detect_bs)

  # Detection-rate agreement at BF < 1/3 threshold (H0 support)
  null_sd   <- bf_sd[both_valid] < 1/3
  null_bs   <- bf_bs[both_valid] < 1/3
  agree_h0  <- mean(null_sd == null_bs)

  c(
    n_both_valid  = sum(both_valid),
    n_sd_only     = sum(!is.na(bf_sd) & is.na(bf_bs)),
    n_bs_only     = sum(is.na(bf_sd) & !is.na(bf_bs)),
    r_log10       = r_log10,
    mae_log10     = mae_log10,
    bias_log10    = bias_log10,
    P_BF10_gt3_SD = mean(bf_sd > 3, na.rm = TRUE),
    P_BF10_gt3_BS = mean(bf_bs > 3, na.rm = TRUE),
    P_BF01_gt3_SD = mean(bf_sd < 1/3, na.rm = TRUE),
    P_BF01_gt3_BS = mean(bf_bs < 1/3, na.rm = TRUE),
    agree_h1      = agree_h1,
    agree_h0      = agree_h0,
    mean_true_p0  = mean(results[, "true_p0"], na.rm = TRUE)
  )
}


# ------------------------------------------------------------------------------
# Run simulation
# ------------------------------------------------------------------------------

smoke_test <- FALSE

out_file <- here("output", "res_bf_validation.rds")

if (smoke_test) {
  test_design <- Design[1, , drop = FALSE]
  res_test <- runSimulation(
    design        = test_design,
    replications  = 2,
    generate      = Generate,
    analyse       = Analyse,
    summarise     = Summarise,
    fixed_objects = list(
      b0_range    = c(0.4, 0.9),
      sd_b0       = 0.75,
      true_sd_re  = 0.25,
      sd_prior_re = "exponential",
      dist_b0     = "logistic",
      n_subjects  = 60,
      n_trials    = 50,
      n_cores     = 4
    ),
    parallel = FALSE
  )
  message("Smoke test complete.")
  print(res_test)

} else if (!file.exists(out_file)) {

  # --- Parallel cluster setup ---
  ncores_total <- parallel::detectCores() - 2
  ncores_brms  <- 2L                                   # cores per brms fit
  ncores_sim   <- max(1L, floor(ncores_total / ncores_brms))

  # Resolve paths in main process (avoids here() issues on workers)
  src_link <- here("R", "link_functions.R")
  src_bf   <- here("R", "bf_helpers.R")

  cl <- parallel::makeCluster(ncores_sim)
  parallel::clusterExport(cl, c("src_link", "src_bf"), envir = environment())
  parallel::clusterEvalQ(cl, {
    source(src_link)
    source(src_bf)
  })

  res <- runSimulation(
    design        = Design,
    replications  = 100,
    generate      = Generate,
    analyse       = Analyse,
    summarise     = Summarise,
    fixed_objects = list(
      b0_range    = c(0.4, 0.9),
      sd_b0       = 0.75,
      true_sd_re  = 0.25,
      sd_prior_re = "exponential",
      dist_b0     = "logistic",
      n_subjects  = 60,
      n_trials    = 50,
      n_cores     = ncores_brms
    ),
    save_results  = TRUE,
    save_details  = list(
      safe                  = TRUE,
      out_rootdir           = here("output"),
      save_results_dirname  = "Simulation_BFValidation",
      save_results_filename = "BFValid_Cond"
    ),
    cl       = cl,
    packages = c("brms", "bridgesampling", "posterior", "cmdstanr")
  )
  parallel::stopCluster(cl)
  save(res, file = out_file)
  message("Simulation complete. Results saved to: ", out_file)

} else {
  message("Results file already exists. Load with: load('", out_file, "')")
}
