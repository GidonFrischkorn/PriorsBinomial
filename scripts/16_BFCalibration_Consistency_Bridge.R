# ==============================================================================
# Script 16: Consistency simulation re-run with bridge sampling (PsychMethods R1)
# ==============================================================================
# Purpose: Information consistency (BF grows with the effect at fixed N) and
# Lindley consistency (BF grows with N at fixed effect) for the Logistic(0, 0.25)
# effect prior, as in Script 03c, but with the bridge-sampling engine of
# R/bf_helpers.R (one compiled rstan model per hypothesis, refits via
# update(recompile = FALSE); H0 fitted once per dataset).
#
# Design (identical to Script 03c): between-subjects, +/-1 condition coding,
# one aggregated row per subject, random intercepts (true SD 0.25,
# Exponential(4) prior), Logistic(0, 0.75) intercept prior, logit link,
# baseline probability ~ Uniform(0.4, 0.9); true_b1 in {0, .05, .1, .2, .3, .5},
# n_subjects in {20, 30, 50, 100, 200} per condition, n_trials = 20;
# 50 replications per cell. The second prior parameter is a SCALE
# (Logistic SD = scale * pi / sqrt(3)); the design column is `scale_b1`.
#
# A failed fit (H0 or H1, sampling or bridge) is recorded as NA and counted in
# n_fail (allow_na = TRUE); it is never mapped to a number.
#
# Run modes (first command-line argument; default "full"):
#   smoke   one cell (n_subjects = 20, true_b1 = 0.20), 1 replication, 2 cores
#           (PB_NCORES=1 runs it serially, e.g. where workers cannot open a
#           local socket);
#           writes output/res_bf_calibration_consistency_bridge_smoke.rds and
#           output/bf_calibration_consistency_bridge_raw_smoke.rds
#   full    all 30 cells x PB_REPS (default 50) replications; guarded (skipped
#           when the result file exists) and resumable through SimDesign's temp
#           file in the output directory
# Environment: PB_SEED (default 20261104; per-row seeds via genSeeds() on the
#   Design), PB_NCORES (workers; SLURM_CPUS_PER_TASK wins; default
#   detectCores() - 2), PB_REPS, PB_N_SAMPLING (post-warmup draws per chain,
#   default 5000 as in Script 13), CMDSTAN (unused here: rstan backend).
#
# Output: output/Simulation_BFCalibration_Consistency_Bridge/
#           res_bf_calibration_consistency_bridge.rds   SimDesign object
#           SIMDESIGN-TEMPFILE_consistency_bridge        resume file (deleted on completion)
#         output/bf_calibration_consistency_bridge_raw.rds   one row per replication:
#           dist_b0, dist_b1, scale_b1, true_b1, n_subjects, n_trials, design_row, rep,
#           bf10_bridge, bf10_logspline (cross-check from the same H1 draws),
#           log10_bf10 (bridge), pct_error_h1, post_mean, post_sd, n_div, rhat_max,
#           ess_bulk_min, ess_tail_min, wall_s (H1 fit + bridge), wall_s_h0,
#           n_fail (0/1), rep_wall_s, true_p0, data_sum_y
#         output/bf_calibration_consistency_bridge_summaries.rds  one row per cell
#           (Summarise() output: median_log10_bf10, P_bf10_gt3, P_bf10_gt10,
#           P_bf01_gt3, P_bf01_gt10, n_failed, n_valid, diagnostics)
# ==============================================================================

Sys.setenv(OMP_NUM_THREADS = "1")

library(SimDesign)
library(brms)
library(here)

source(here("R", "link_functions.R"))
source(here("R", "bf_helpers.R"))

args     <- commandArgs(trailingOnly = TRUE)
run_mode <- if (length(args) > 0) args[1] else "full"
stopifnot(run_mode %in% c("smoke", "full"))

# ------------------------------------------------------------------------------
# Design
# ------------------------------------------------------------------------------

Design <- createDesign(
  dist_b0    = "logistic",
  dist_b1    = "logistic",
  scale_b1   = 0.25,
  true_b1    = c(0.00, 0.05, 0.10, 0.20, 0.30, 0.50),
  n_subjects = c(20, 30, 50, 100, 200),
  n_trials   = 20
)
Design$design_row <- seq_len(nrow(Design))

b0_range    <- c(0.4, 0.9)
true_sd_re  <- 0.25
scale_b0    <- 0.75
sd_prior_re <- "exponential"

# ------------------------------------------------------------------------------
# Data generation (as in Script 03c)
# ------------------------------------------------------------------------------

#' Simulate one between-subjects dataset with random intercepts
#'
#' @param n_subjects Subjects per condition.
#' @param n_trials Binomial trials per subject.
#' @param true_b1 True effect on the logit scale (+/-1 coding).
#' @param b0_range Range of the baseline probability (uniform).
#' @param true_sd_re SD of the random intercepts.
#' @return Data frame with `subject_id`, `y`, `n`, `condition`, `true_p0`.
simulate_consistency_data <- function(n_subjects, n_trials, true_b1, b0_range, true_sd_re) {
  true_p0 <- runif(1, b0_range[1], b0_range[2])
  true_b0 <- apply_link(true_p0, "logit")
  n_total <- 2L * n_subjects
  cond    <- rep(c(1L, -1L), each = n_subjects)
  u       <- rnorm(n_total, mean = 0, sd = true_sd_re)
  p_subj  <- apply_inverse_link(true_b0 + true_b1 * cond + u, link = "logit")
  data.frame(subject_id = seq_len(n_total),
             y          = rbinom(n_total, n_trials, p_subj),
             n          = n_trials,
             condition  = cond,
             true_p0    = true_p0)
}

# ------------------------------------------------------------------------------
# SimDesign functions
# ------------------------------------------------------------------------------

Generate <- function(condition, fixed_objects = NULL) {
  Attach(condition)
  fixed_objects$simulate_consistency_data(n_subjects, n_trials, true_b1,
                                          fixed_objects$b0_range,
                                          fixed_objects$true_sd_re)
}

Analyse <- function(condition, dat, fixed_objects = NULL) {
  Attach(condition)
  if (!exists("fit_bf_bs_rs", mode = "function")) {
    for (f in fixed_objects$helper_files) sys.source(f, envir = globalenv())
  }
  t_start <- Sys.time()
  mcmc  <- fixed_objects$mcmc
  seeds <- sample.int(.Machine$integer.max - 1L, 2)   # H1, H0
  out_names <- c(fit_bf_rs_names("bridge"), "wall_s_h0")
  res <- tryCatch({
    h0 <- do.call(fit_h0_bridge_rs,
                  c(list(fixed_objects$templates$h0, dat, seed = seeds[2]), mcmc))
    r  <- do.call(fit_bf_bs_rs,
                  c(list(fixed_objects$templates$h1, dat, dist_b1 = as.character(dist_b1),
                         scale_b1 = scale_b1, ml_h0 = h0$ml_h0, predictor = "condition",
                         seed = seeds[1]), mcmc))
    c(r, wall_s_h0 = h0$wall_h0)
  }, error = function(e) {
    warning("fit failed: ", conditionMessage(e), call. = FALSE)
    setNames(rep(NA_real_, length(out_names)), out_names)
  })
  c(res,
    n_fail     = as.integer(is.na(res[["bf10_bridge"]])),
    rep_wall_s = as.numeric(difftime(Sys.time(), t_start, units = "secs")),
    true_p0    = dat$true_p0[1],
    data_sum_y = sum(dat$y))
}

Summarise <- function(condition, results, fixed_objects = NULL) {
  results <- as.data.frame(results)
  bf10 <- results$bf10_bridge
  safe_max  <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  safe_min  <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
  safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  c(median_log10_bf10 = stats::median(log10(bf10)[is.finite(log10(bf10))]),
    P_bf10_gt3   = safe_mean(bf10 > 3),
    P_bf10_gt10  = safe_mean(bf10 > 10),
    P_bf01_gt3   = safe_mean(1 / bf10 > 3),
    P_bf01_gt10  = safe_mean(1 / bf10 > 10),
    n_failed     = sum(is.na(bf10)),
    n_valid      = sum(!is.na(bf10)),
    n_inf        = sum(is.infinite(bf10)),
    max_rhat     = safe_max(results$rhat_max),
    min_ess_bulk = safe_min(results$ess_bulk_min),
    min_ess_tail = safe_min(results$ess_tail_min),
    share_any_div    = safe_mean(results$n_div > 0),
    max_pct_error_h1 = safe_max(results$pct_error_h1),
    mean_wall_s      = safe_mean(results$wall_s),
    mean_wall_s_h0   = safe_mean(results$wall_s_h0),
    mean_rep_wall_s  = mean(results$rep_wall_s),
    mean_true_p0     = mean(results$true_p0))
}

# ------------------------------------------------------------------------------
# Templates: one compiled rstan model per hypothesis, random intercepts only
# ------------------------------------------------------------------------------

t0 <- Sys.time()
template_dat <- withr::with_seed(1, simulate_consistency_data(20, 20, 0.1, b0_range, true_sd_re))
templates <- lapply(c(h1 = "h1", h0 = "h0"), function(h) {
  compile_bf_template_rs(template_dat, link = "logit", dist_b0 = "logistic",
                         scale_b0 = scale_b0, predictor = "condition", hypothesis = h,
                         sd_prior_re = sd_prior_re, backend = "rstan",
                         random_slope = FALSE)
})
message("rstan templates compiled in ",
        round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s")
# The templates must not contain a random slope (between-subjects design)
stopifnot(!any(grepl("sd_subject_id__condition", brms::variables(templates$h1))))

fixed_cfg <- list(
  templates    = templates,
  simulate_consistency_data = simulate_consistency_data,
  helper_files = c(here("R", "link_functions.R"), here("R", "bf_helpers.R")),
  b0_range     = b0_range,
  true_sd_re   = true_sd_re,
  mcmc         = list(n_chains = 4, n_warmup = 1000,
                      n_sampling = as.integer(Sys.getenv("PB_N_SAMPLING", "5000")),
                      adapt_delta = 0.95)
)
sim_packages <- c("brms", "posterior", "logspline", "bridgesampling", "withr")
sim_control  <- list(allow_na = TRUE, allow_nan = TRUE)
global_seed  <- as.integer(Sys.getenv("PB_SEED", "20261104"))
seeds_all    <- genSeeds(Design, iseed = global_seed)

#' Flatten a SimDesign result to one row per replication
#'
#' @param res SimDesign result with stored replications.
#' @return Data frame with the design columns, `rep` and all Analyse outputs
#'   plus `log10_bf10`.
flatten_consistency_bridge <- function(res) {
  raw <- as.data.frame(SimResults(res))
  raw <- raw[order(raw$design_row), ]
  raw$rep <- stats::ave(seq_len(nrow(raw)), raw$design_row, FUN = seq_along)
  raw$log10_bf10 <- log10(raw$bf10_bridge)
  design_cols <- c("dist_b0", "dist_b1", "scale_b1", "true_b1", "n_subjects",
                   "n_trials", "design_row", "rep")
  raw[, c(design_cols, setdiff(names(raw), design_cols))]
}

# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------

if (run_mode == "smoke") {
  sel <- which(Design$n_subjects == 20 & Design$true_b1 == 0.20)
  nc_smoke <- suppressWarnings(as.integer(Sys.getenv("PB_NCORES", "2")))
  if (is.na(nc_smoke) || nc_smoke < 1) nc_smoke <- 2L
  message("Smoke: design row ", sel, ", 1 replication, ", nc_smoke, " core(s)")
  res_s <- runSimulation(
    design = Design[sel, ], replications = 1, generate = Generate,
    analyse = Analyse, summarise = Summarise, fixed_objects = fixed_cfg,
    parallel = nc_smoke > 1, ncores = nc_smoke, packages = sim_packages,
    control = sim_control, seed = seeds_all[sel], save = FALSE
  )
  saveRDS(res_s, here("output", "res_bf_calibration_consistency_bridge_smoke.rds"))
  raw_s <- flatten_consistency_bridge(res_s)
  saveRDS(raw_s, here("output", "bf_calibration_consistency_bridge_raw_smoke.rds"))
  cat(sprintf(paste0("SMOKE_CONSISTENCY_BRIDGE ns=%d b1=%.2f | log10 BF10 bridge=%.3f ",
                     "logspline=%.3f | pct_err=%.2f | n_fail=%d div=%d rhat_max=%.4f ",
                     "ess_bulk_min=%.0f ess_tail_min=%.0f | wall H1=%.1fs H0=%.1fs rep=%.1fs\n"),
              raw_s$n_subjects, raw_s$true_b1, raw_s$log10_bf10, log10(raw_s$bf10_logspline),
              raw_s$pct_error_h1, raw_s$n_fail, raw_s$n_div, raw_s$rhat_max,
              raw_s$ess_bulk_min, raw_s$ess_tail_min, raw_s$wall_s, raw_s$wall_s_h0,
              raw_s$rep_wall_s))
  print(names(raw_s))
  print(SimExtract(res_s, what = "errors"))
  print(SimExtract(res_s, what = "warnings"))
  message("SMOKE CONSISTENCY BRIDGE COMPLETE")

} else {
  out_dir  <- here("output", "Simulation_BFCalibration_Consistency_Bridge")
  res_file <- file.path(out_dir, "res_bf_calibration_consistency_bridge.rds")
  n_reps   <- as.integer(Sys.getenv("PB_REPS", "50"))
  ncores   <- resolve_ncores()
  if (file.exists(res_file)) {
    message("Result exists, simulation skipped: ", res_file)
    res <- readRDS(res_file)
  } else {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    message(sprintf("Consistency (bridge): %d cells x %d reps, %d cores, seed %d; brms fits = %d",
                    nrow(Design), n_reps, ncores, global_seed, 2L * nrow(Design) * n_reps))
    res <- runSimulation(
      design        = Design,
      replications  = n_reps,
      generate      = Generate,
      analyse       = Analyse,
      summarise     = Summarise,
      fixed_objects = fixed_cfg,
      seed          = seeds_all,
      save          = TRUE,
      resume        = TRUE,
      store_results = TRUE,
      save_results  = FALSE,
      save_details  = list(safe = TRUE, out_rootdir = out_dir,
                           tmpfilename = "SIMDESIGN-TEMPFILE_consistency_bridge"),
      parallel = ncores > 1,
      ncores   = ncores,
      packages = sim_packages,
      control  = sim_control,
      verbose  = TRUE,
      progress = FALSE
    )
    saveRDS(res, res_file)
    message("Simulation complete: ", res_file)
  }
  raw <- flatten_consistency_bridge(res)
  saveRDS(raw, here("output", "bf_calibration_consistency_bridge_raw.rds"))
  saveRDS(as.data.frame(res), here("output", "bf_calibration_consistency_bridge_summaries.rds"))
  message(sprintf("Saved raw (%d replications, %d failed) and summaries (%d cells) to output/",
                  nrow(raw), sum(raw$n_fail), nrow(res)))
  message("CONSISTENCY BRIDGE COMPLETE")
}
