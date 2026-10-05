# ==============================================================================
# Script 13: BF Calibration with Random Slopes (revision, PsychMethods R1)
# ==============================================================================
# Purpose: Re-calibrate BF10 for the effect prior when the analysis model
# matches the random-effects structure of a within-subject design:
# y | trials(n) ~ 1 + condition + (1 + condition | subject_id).
#
# Paired design: in every replication ALL effect priors (family x scale) of an
# arm are fitted to the SAME dataset (as in 03d). The effect prior is passed to
# Stan as data (stanvars), so one compiled model per arm serves all priors;
# refits use update(recompile = FALSE).
#
# Arms (Design factor `arm`):
#   logit_factor     logit link, Logistic(0, 0.75) intercept, +/-1 factor
#   probit_factor    probit link, Normal(0, 0.75) intercept, +/-1 factor
#   logit_continuous logit link, Logistic(0, 0.75) intercept, trial-level
#                    x ~ N(0, 1), Bernoulli data (not aggregated)
# Exponential(4) on both random-effect SDs in all arms.
#
# Full grid (approved 2026-10-05: grid C + logit slope SD 0.50, 200 reps),
# defined in make_design_rs() / make_prior_grids_rs() in R/bf_helpers.R:
#   logit_factor     21 priors (3 families x 7 scales); true_b1 x n_subjects x
#                    n_trials {20, 50} x slope SD {0.10, 0.30, 0.50} = 54 rows
#   probit_factor    Normal x 7 scales; n_trials = 20, slope SD 0.30 = 9 rows
#   logit_continuous 3 families x {0.15, 0.25, 0.30, 0.50}; n_trials = 20,
#                    slope SD 0.30 = 9 rows
# BF engine: bridge sampling on 4 x 5,000 post-warmup draws (H0 fitted once per
# dataset), logspline Savage-Dickey from the same H1 draws as a cross-check.
#
# The second prior parameter is a SCALE (Logistic SD = scale * pi / sqrt(3)).
#
# Run modes (first command-line argument; default "full"):
#   smoke         Savage-Dickey engine, 2 reps, 4 cells, 4 cores
#   smoke_bridge  bridge engine, 1 rep of one cell per arm, 3 priors per arm
#   chunkcheck    rows 1-2 as one chunk vs. two one-row chunks, same seed:
#                 data and BFs must be identical
#   errcheck      injected prior failures: one -> warning, most -> error
#   benchmark     BF-method benchmark (Savage-Dickey variants vs bridge)
#   full          full simulation or one chunk of it (guarded; resumable)
#
# Environment variables for `full` (all optional):
#   PB_ROWS        Design rows to run, e.g. "1-12" or "3,5,60-72" (or 2nd arg)
#   PB_CHUNK, PB_N_CHUNKS  run contiguous chunk PB_CHUNK of PB_N_CHUNKS
#                  (SLURM_ARRAY_TASK_ID is used for PB_CHUNK if set)
#   PB_REPS        replications per row (default 200)
#   PB_SEED        global seed (default 20261005); per-row seeds come from
#                  genSeeds() on the FULL Design, so results do not depend on
#                  chunking
#   PB_NCORES      workers (SLURM_CPUS_PER_TASK wins; default detectCores()-2)
#   PB_OUT_ROOT    output root (default output/Simulation_BFCalibration_RS)
#   CMDSTAN        CmdStan path (default cmdstanr::cmdstan_path())
# Each chunk writes to output/Simulation_BFCalibration_RS/rows_<spec>/;
# scripts/14_Merge_BFCalibration_RS.R combines the chunks.
# ==============================================================================

# One OpenMP thread per Stan chain (see Script 05).
Sys.setenv(OMP_NUM_THREADS = "1")

library(SimDesign)
library(brms)
library(here)

source(here("R", "link_functions.R"))
source(here("R", "bf_helpers.R"))

args     <- commandArgs(trailingOnly = TRUE)
run_mode <- if (length(args) > 0) args[1] else "full"
stopifnot(run_mode %in% c("smoke", "smoke_bridge", "chunkcheck", "errcheck",
                          "benchmark", "full"))

if (nzchar(Sys.getenv("CMDSTAN"))) cmdstanr::set_cmdstan_path(Sys.getenv("CMDSTAN"))
message("CmdStan ", cmdstanr::cmdstan_version(), " at ", cmdstanr::cmdstan_path())

# ------------------------------------------------------------------------------
# Design and effect-prior grids (shared with Script 14)
# ------------------------------------------------------------------------------

Design      <- make_design_rs()
Design_all  <- make_design_rs(full_crossing = TRUE)  # for test cells only
prior_grids <- make_prior_grids_rs()

arm_cfg <- list(
  logit_factor     = list(link = "logit",  dist_b0 = "logistic", predictor = "condition"),
  probit_factor    = list(link = "probit", dist_b0 = "normal",   predictor = "condition"),
  logit_continuous = list(link = "logit",  dist_b0 = "logistic", predictor = "x")
)

# ------------------------------------------------------------------------------
# Data generation
# ------------------------------------------------------------------------------

#' Simulate one within-subject dataset with random intercepts and slopes
#'
#' @param n_subjects,n_trials,true_b1,true_sd_slope Design values.
#' @param link,predictor Arm settings ("condition" or "x").
#' @param b0_range Range of the baseline probability (uniform).
#' @param true_sd_re SD of the random intercepts.
#' @return Data frame with `subject_id`, `y`, `n`, the predictor, `true_p0`.
simulate_rs_data <- function(n_subjects, n_trials, true_b1, true_sd_slope,
                             link, predictor, b0_range, true_sd_re) {
  true_p0 <- runif(1, b0_range[1], b0_range[2])
  true_b0 <- apply_link(true_p0, link)
  u_int   <- rnorm(n_subjects, mean = 0, sd = true_sd_re)
  u_slope <- rnorm(n_subjects, mean = 0, sd = true_sd_slope)

  if (predictor == "condition") {
    # One aggregated row per subject x condition (as in 03d)
    subj <- rep(seq_len(n_subjects), times = 2)
    pred <- rep(c(1L, -1L), each = n_subjects)
    n    <- n_trials
  } else {
    # Trial-level standardized predictor, Bernoulli outcome
    subj <- rep(seq_len(n_subjects), each = 2L * n_trials)
    pred <- rnorm(length(subj))
    n    <- 1L
  }
  eta <- true_b0 + u_int[subj] + (true_b1 + u_slope[subj]) * pred
  p   <- apply_inverse_link(eta, link = link)

  out <- data.frame(subject_id = subj,
                    y          = rbinom(length(subj), n, p),
                    n          = n,
                    true_p0    = true_p0)
  out[[predictor]] <- pred
  out
}

# ------------------------------------------------------------------------------
# SimDesign functions
# ------------------------------------------------------------------------------

Generate <- function(condition, fixed_objects = NULL) {
  Attach(condition)
  # Workers only receive globals when runSimulation() is called at top level;
  # source the helpers on the worker if they are missing (e.g. run_chunk())
  if (!exists("fit_prior_grid_rs", mode = "function")) {
    for (f in fixed_objects$helper_files) sys.source(f, envir = globalenv())
  }
  cfg <- fixed_objects$arm_cfg[[arm]]
  # Generator passed via fixed_objects so workers have it in any calling context
  fixed_objects$simulate_rs_data(n_subjects    = n_subjects,
                                 n_trials      = n_trials,
                                 true_b1       = true_b1,
                                 true_sd_slope = true_sd_slope,
                                 link          = cfg$link,
                                 predictor     = cfg$predictor,
                                 b0_range      = fixed_objects$b0_range,
                                 true_sd_re    = fixed_objects$true_sd_re)
}

Analyse <- function(condition, dat, fixed_objects = NULL) {
  Attach(condition)
  if (!exists("fit_prior_grid_rs", mode = "function")) {
    for (f in fixed_objects$helper_files) sys.source(f, envir = globalenv())
  }
  cfg     <- fixed_objects$arm_cfg[[arm]]
  grid    <- fixed_objects$prior_grids[[arm]]
  t_start <- Sys.time()
  # One Stan/bridge seed per prior plus one for H0, drawn from the
  # replication's SimDesign RNG stream (reproducible under runSimulation(seed))
  seeds <- sample.int(.Machine$integer.max - 1L, nrow(grid) + 1)
  args <- list(
    template    = fixed_objects$templates[[arm]],
    dat         = dat,
    prior_grid  = grid,
    predictor   = cfg$predictor,
    engine      = fixed_objects$bf_engine,
    seeds       = seeds,
    n_chains    = fixed_objects$n_chains,
    n_warmup    = fixed_objects$n_warmup,
    n_sampling  = fixed_objects$n_sampling,
    adapt_delta = fixed_objects$adapt_delta
  )
  if (fixed_objects$bf_engine == "sd") {
    args$bf_methods <- fixed_objects$bf_methods
  } else {
    args$template_h1_bs <- fixed_objects$templates_bs[[arm]]$h1
    args$template_h0_bs <- fixed_objects$templates_bs[[arm]]$h0
  }
  res <- do.call(fit_prior_grid_rs, args)  # warns / stops on failures
  c(res,
    rep_wall_s = as.numeric(difftime(Sys.time(), t_start, units = "secs")),
    true_p0    = dat$true_p0[1],
    data_sum_y = sum(dat$y))          # fingerprint of the simulated data
}

Summarise <- function(condition, results, fixed_objects = NULL) {
  results <- as.data.frame(results)
  safe_max  <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  safe_min  <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
  safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  pick <- function(prefix) grep(paste0("^", prefix, "_(normal|logistic|cauchy)_"),
                                names(results), value = TRUE)
  stat <- function(cols, fun, label) {
    if (length(cols) == 0) return(numeric(0))  # e.g. no pct_error_h1 for engine "sd"
    out <- vapply(cols, function(col) fun(results[[col]]), numeric(1))
    names(out) <- paste0(label, "_", cols)
    out
  }
  bf_cols <- grep("^bf10_", names(results), value = TRUE)
  summ_bf <- unlist(lapply(bf_cols, function(col) {
    bf10 <- results[[col]]
    setNames(c(median(log10(bf10), na.rm = TRUE),
               safe_mean(bf10 > 3),
               safe_mean(bf10 < 1 / 3),
               sum(is.na(bf10)),
               sum(is.infinite(bf10))),
             paste0(c("median_log10_", "P_gt3_", "P_lt1_3_", "n_na_", "n_inf_"), col))
  }))
  summ_diag <- c(
    stat(pick("rhat_max"),     safe_max,                         "max"),
    stat(pick("rhat_max"),     function(x) safe_mean(x > 1.01),  "share_gt1.01"),
    stat(pick("ess_bulk_min"), safe_min,                         "min"),
    stat(pick("ess_tail_min"), safe_min,                         "min"),
    stat(pick("n_div"),        function(x) safe_mean(x > 0),     "share_any"),
    stat(pick("pct_error_h1"), safe_max,                         "max"),
    stat(pick("pct_error_h1"), function(x) safe_mean(x > 5),     "share_gt5"),
    stat(pick("wall_s"),       safe_mean,                        "mean")
  )
  c(summ_bf, summ_diag,
    total_n_fail    = sum(results$n_fail),
    share_reps_fail = mean(results$n_fail > 0),
    mean_wall_s_h0  = if (is.null(results$wall_s_h0)) NA_real_ else safe_mean(results$wall_s_h0),
    mean_rep_wall_s = mean(results$rep_wall_s),
    max_rep_wall_s  = max(results$rep_wall_s),
    mean_true_p0    = mean(results$true_p0))
}

# ------------------------------------------------------------------------------
# Templates: one compiled model per arm (compiled once, reused by workers)
# ------------------------------------------------------------------------------

b0_range   <- c(0.4, 0.9)
true_sd_re <- 0.25

build_templates <- function(arms, backend = "cmdstanr", hypothesis = "h1") {
  setNames(lapply(arms, function(a) {
    cfg <- arm_cfg[[a]]
    # Template data only fix the data structure; values are irrelevant
    dat <- withr::with_seed(1, simulate_rs_data(20, 20, 0.1, 0.1, cfg$link,
                                                cfg$predictor, b0_range, true_sd_re))
    compile_bf_template_rs(dat, link = cfg$link, dist_b0 = cfg$dist_b0,
                           scale_b0 = 0.75, predictor = cfg$predictor,
                           hypothesis = hypothesis, sd_prior_re = "exponential",
                           backend = backend)
  }), arms)
}

# Compiled in every session: the executables live in the session tempdir, so
# a cached template .rds goes stale once the compiling R session ends.
bf_engine <- if (run_mode %in% c("smoke", "errcheck")) "sd" else "bridge"

t0 <- Sys.time()
templates <- build_templates(names(arm_cfg))
message("cmdstanr templates compiled in ",
        round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s")

templates_bs <- NULL
if (bf_engine == "bridge") {
  t0 <- Sys.time()
  templates_bs <- setNames(lapply(names(arm_cfg), function(a) {
    list(h1 = build_templates(a, backend = "rstan", hypothesis = "h1")[[1]],
         h0 = build_templates(a, backend = "rstan", hypothesis = "h0")[[1]])
  }), names(arm_cfg))
  message("rstan templates compiled in ",
          round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s")
}

fixed_cfg <- list(
  arm_cfg          = arm_cfg,
  prior_grids      = prior_grids,
  templates        = templates,
  templates_bs     = templates_bs,
  simulate_rs_data = simulate_rs_data,
  helper_files     = c(here("R", "link_functions.R"), here("R", "bf_helpers.R")),
  bf_engine        = bf_engine,
  b0_range         = b0_range,
  true_sd_re       = true_sd_re,
  bf_methods       = c("kde", "logspline", "normal"),
  n_chains         = 4,
  n_warmup         = 1000,
  n_sampling       = if (bf_engine == "bridge") 5000 else 1000,
  adapt_delta      = 0.95
)

sim_packages <- c("brms", "posterior", "logspline", "bridgesampling", "withr")
# NA must be a valid Analyse result: a failed prior fit is recorded as NA (and
# counted in n_fail). SimDesign's default treats any NA as an error and redraws
# the whole dataset, which would silently discard and replace hard datasets.
sim_control <- list(allow_na = TRUE, allow_nan = TRUE)

global_seed <- as.integer(Sys.getenv("PB_SEED", "20261005"))
# Per-row L'Ecuyer seeds from the FULL Design: row k always gets the same seed,
# whichever chunk it runs in
seeds_all <- genSeeds(Design, iseed = global_seed)

#' Run a set of Design rows as one chunk with its own output directory
#'
#' @param rows Integer vector of `design_row` IDs.
#' @param reps Integer replications per row.
#' @param root Directory under which `rows_<label>/` is created.
#' @param fixed SimDesign fixed_objects.
#' @param ncores Integer number of workers.
#' @return The SimDesign result (also saved as `res_rows_<label>.rds`), or
#'   `NULL` if that file already exists.
run_chunk <- function(rows, reps, root, fixed, ncores) {
  label     <- if (length(rows) == 1) sprintf("%03d", rows) else
    if (all(diff(rows) == 1)) sprintf("%03d-%03d", min(rows), max(rows)) else
      paste(sprintf("%03d", rows), collapse = "_")
  chunk_dir <- file.path(root, paste0("rows_", label))
  res_file  <- file.path(chunk_dir, paste0("res_rows_", label, ".rds"))
  # Guard. An interrupted chunk resumes: with save = TRUE / resume = TRUE
  # SimDesign keeps finished conditions in its temp file (tmpfilename below,
  # in chunk_dir), so re-running the same command skips them.
  if (file.exists(res_file)) {
    message("Chunk already complete: ", res_file)
    return(invisible(NULL))
  }
  dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)
  d <- Design[match(rows, Design$design_row), ]
  message(sprintf("Chunk %s: %d rows, %d reps, %d cores, expected brms fits (H1 + H0) = %d",
                  label, nrow(d), reps, ncores,
                  sum(vapply(fixed$prior_grids, nrow, integer(1))[d$arm] + 1L) * reps))
  res <- runSimulation(
    design        = d,
    replications  = reps,
    generate      = Generate,
    analyse       = Analyse,
    summarise     = Summarise,
    fixed_objects = fixed,
    seed          = seeds_all[rows],
    save          = TRUE,
    resume        = TRUE,
    # Raw replications are stored inside the result object (SimResults()).
    # save_results = TRUE is not used: SimResults() then re-reads the files
    # through setwd(<wd>/<out_rootdir>) and fails for an absolute out_rootdir.
    store_results = TRUE,
    save_results  = FALSE,
    save_details  = list(
      safe        = TRUE,
      out_rootdir = chunk_dir,
      tmpfilename = paste0("SIMDESIGN-TEMPFILE_rows_", label)
    ),
    parallel = ncores > 1,
    ncores   = ncores,
    packages = sim_packages,
    control  = sim_control,
    verbose  = TRUE,       # one progress line per Design row in the log
    progress = FALSE
  )
  saveRDS(res, res_file)
  message("Chunk complete: ", res_file)
  invisible(res)
}

#' One log line per replication of a small test run
#'
#' @param res A SimDesign result object.
#' @param tag Character prefix for the log line.
#' @return Invisibly, the raw results.
print_rep_lines <- function(res, tag) {
  raw <- SimResults(res)
  col <- function(r, pre) unlist(r[grep(paste0("^", pre, "_(normal|logistic|cauchy)"), names(r))])
  for (i in seq_len(nrow(raw))) {
    # Arms have different prior grids, so SimResults() pads other arms'
    # columns with NA; keep this row's own columns only
    r <- raw[i, ]
    r <- r[, vapply(r, function(v) !all(is.na(v)), logical(1))]
    bf_main <- if (any(grepl("^bf10_bridge_", names(r)))) "bf10_bridge" else "bf10_kde"
    cat(sprintf(paste0("%s arm=%s ns=%d nt=%d b1=%.2f slopeSD=%.2f | rep_wall=%.1fs ",
                       "fit_wall mean=%.1fs | n_fail=%d div=%d rhat_max=%.4f ",
                       "ess_bulk_min=%.0f ess_tail_min=%.0f pct_err_max=%s | ",
                       "log10 %s: %s\n"),
                tag, r$arm, r$n_subjects, r$n_trials, r$true_b1, r$true_sd_slope,
                r$rep_wall_s, mean(col(r, "wall_s")), as.integer(r$n_fail),
                as.integer(sum(col(r, "n_div"))), max(col(r, "rhat_max")),
                min(col(r, "ess_bulk_min")), min(col(r, "ess_tail_min")),
                if (bf_main == "bf10_bridge") sprintf("%.2f", max(col(r, "pct_error_h1"))) else "NA",
                bf_main, paste(round(log10(col(r, bf_main)), 3), collapse = " ")))
  }
  invisible(raw)
}

#' Three priors spread over an arm's grid (first, middle, last row)
#'
#' @param g Prior-grid data frame.
#' @return Subset of `g` with up to three rows.
spread_priors <- function(g) g[unique(round(seq(1, nrow(g), length.out = 3))), ]

# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------

if (run_mode == "smoke") {
  # Savage-Dickey engine; smallest and largest logit cells, one probit cell,
  # the largest continuous cell (full crossing, so not all are in Design)
  sel <- with(Design_all,
    (arm == "logit_factor" & n_subjects == 20 & n_trials == 20 & true_b1 == 0 & true_sd_slope == 0.1) |
    (arm == "logit_factor" & n_subjects == 50 & n_trials == 50 & true_b1 == 0.2 & true_sd_slope == 0.3) |
    (arm == "probit_factor" & n_subjects == 30 & n_trials == 50 & true_b1 == 0.1 & true_sd_slope == 0.3) |
    (arm == "logit_continuous" & n_subjects == 50 & n_trials == 50 & true_b1 == 0.1 & true_sd_slope == 0.3))
  res_test <- runSimulation(
    design = Design_all[sel, ], replications = 2, generate = Generate,
    analyse = Analyse, summarise = Summarise, fixed_objects = fixed_cfg,
    parallel = TRUE, ncores = 4, packages = sim_packages, control = sim_control,
    seed = genSeeds(Design_all[sel, ], iseed = 1301), save = FALSE
  )
  saveRDS(res_test, here("output", "res_bf_calibration_rs_smoke.rds"))
  print_rep_lines(res_test, "SMOKE")
  print(SimExtract(res_test, what = "errors"))
  print(SimExtract(res_test, what = "warnings"))
  message("SMOKE TEST COMPLETE")

} else if (run_mode == "smoke_bridge") {
  # Bridge engine through SimDesign: 1 replication of the largest Design row
  # per arm (n_subjects = 50, n_trials = 50 logit / 20 probit and
  # continuous), 3 priors per arm
  fixed_sb <- fixed_cfg
  fixed_sb$prior_grids <- lapply(prior_grids, spread_priors)
  sel <- with(Design, n_subjects == 50 & true_b1 == 0.1 & true_sd_slope == 0.3 &
                (arm != "logit_factor" | n_trials == 50))
  res_sb <- runSimulation(
    design = Design[sel, ], replications = 1, generate = Generate,
    analyse = Analyse, summarise = Summarise, fixed_objects = fixed_sb,
    parallel = TRUE, ncores = 3, packages = sim_packages, control = sim_control,
    seed = seeds_all[Design$design_row[sel]], save = FALSE
  )
  saveRDS(res_sb, here("output", "res_bf_calibration_rs_smoke_bridge.rds"))
  raw <- print_rep_lines(res_sb, "SMOKE_BRIDGE")
  cat("H0 wall (s):", round(raw$wall_s_h0, 1), "\n")
  cat("Summarise diagnostic columns (first 12 of",
      length(grep("^(max|min|share_any|share_gt)", names(res_sb))), "):\n")
  print(head(grep("^(max|min|share_any|share_gt)", names(res_sb), value = TRUE), 12))
  print(as.data.frame(res_sb[, c("arm", "total_n_fail", "share_reps_fail",
                                 grep("^max_rhat_max_|^min_ess_tail_min_", names(res_sb), value = TRUE)[1:2])]))
  print(SimExtract(res_sb, what = "errors"))
  print(SimExtract(res_sb, what = "warnings"))
  cat("\n---- sessionInfo ----\n")
  print(sessionInfo())
  cat("CmdStan version:", cmdstanr::cmdstan_version(), "\n")
  message("SMOKE BRIDGE COMPLETE")

} else if (run_mode == "chunkcheck") {
  # Rows 1-2 run (A) as one chunk and (B) as two one-row chunks, same global
  # seed, 2 reps, 3 priors per arm, bridge engine. Data fingerprints and BFs
  # must be identical. Merge with Script 14 afterwards.
  fixed_c <- fixed_cfg
  fixed_c$prior_grids <- lapply(prior_grids, spread_priors)
  root_a <- here("output", "chunkcheck_RS", "A")
  root_b <- here("output", "chunkcheck_RS", "B")
  nc <- 2
  run_chunk(1:2, 2, root_a, fixed_c, nc)
  run_chunk(1, 2, root_b, fixed_c, nc)
  run_chunk(2, 2, root_b, fixed_c, nc)
  get_raw <- function(f) { r <- SimResults(readRDS(f)); r[order(r$design_row), ] }
  ra <- get_raw(file.path(root_a, "rows_001-002", "res_rows_001-002.rds"))
  rb <- rbind(get_raw(file.path(root_b, "rows_001", "res_rows_001.rds")),
              get_raw(file.path(root_b, "rows_002", "res_rows_002.rds")))
  cols <- c("design_row", "data_sum_y", "true_p0", grep("^bf10_", names(ra), value = TRUE))
  print(as.data.frame(ra[, cols]), digits = 10)
  print(as.data.frame(rb[, cols]), digits = 10)
  same_data <- identical(ra$data_sum_y, rb$data_sum_y) && identical(ra$true_p0, rb$true_p0)
  bf_a <- as.matrix(ra[, grep("^bf10_", names(ra))])
  bf_b <- as.matrix(rb[, grep("^bf10_", names(rb))])
  cat(sprintf("CHUNKCHECK identical data: %s | identical BFs: %s | max |d log10 BF| = %.3g | n BF values = %d\n",
              same_data, identical(bf_a, bf_b), max(abs(log10(bf_a) - log10(bf_b))), length(bf_a)))
  message("CHUNKCHECK COMPLETE")

} else if (run_mode == "errcheck") {
  # Injected failures (negative scale fails inside make_effect_prior_stanvars):
  # run A, 1 of 3 priors fails -> SimDesign warning with message, n_fail = 1;
  # run B, 2 of 3 fail -> replication stopped -> SimDesign error with message
  sel  <- with(Design, arm == "logit_factor" & n_subjects == 20 & n_trials == 20 &
                 true_b1 == 0 & true_sd_slope == 0.3)
  bad  <- data.frame(dist_b1 = c("normal", "logistic", "cauchy"),
                     scale_b1 = c(0.25, -1, 0.707))
  bad2 <- data.frame(dist_b1 = c("normal", "logistic", "cauchy"),
                     scale_b1 = c(0.25, -1, -1))
  for (g in list(A = bad, B = bad2)) {
    fixed_e <- fixed_cfg
    fixed_e$prior_grids$logit_factor <- g
    res_e <- tryCatch(runSimulation(
      design = Design[sel, ], replications = 1, generate = Generate,
      analyse = Analyse, summarise = Summarise, fixed_objects = fixed_e,
      seed = 99L, max_errors = 3, save = FALSE, packages = sim_packages,
      control = sim_control
    ), error = function(e) e)
    if (inherits(res_e, "error")) {
      cat("ERRCHECK runSimulation error:", conditionMessage(res_e), "\n")
    } else {
      cat("ERRCHECK n_fail in results:", SimResults(res_e)$n_fail,
          "| BF NA count:", sum(is.na(unlist(SimResults(res_e)[grep("^bf10_", names(SimResults(res_e)))]))), "\n")
      print(SimExtract(res_e, what = "errors"))
      print(SimExtract(res_e, what = "warnings"))
    }
  }
  message("ERRCHECK COMPLETE")

} else if (run_mode == "benchmark") {
  # BF-method benchmark: Savage-Dickey (KDE / logspline / normal approximation)
  # at 4,000 and 20,000 post-warmup draws versus bridge sampling (20,000 draws,
  # rstan templates). logit_factor arm, n_subjects = 30, n_trials = 50,
  # true slope SD = 0.30.
  bench_file <- here("output", "bf_rs_benchmark.rds")
  tmpl_bs    <- templates_bs$logit_factor
  cfg        <- arm_cfg$logit_factor
  bench_priors <- data.frame(dist_b1  = c("normal", "logistic", "cauchy"),
                             scale_b1 = c(0.15, 0.25, 0.707))
  bench_b1 <- c(0, 0.10, 0.20, 0.35)
  rows <- list()
  for (b1 in bench_b1) for (k in 1:2) {
    ds_seed <- 2026 + 10 * k + round(100 * b1)
    dat <- withr::with_seed(ds_seed,
                            simulate_rs_data(30, 50, b1, 0.30, cfg$link, cfg$predictor,
                                             b0_range, true_sd_re))
    h0 <- fit_h0_bridge_rs(tmpl_bs$h0, dat, seed = ds_seed, n_sampling = 5000)
    for (j in seq_len(nrow(bench_priors))) {
      d <- bench_priors$dist_b1[j]; s <- bench_priors$scale_b1[j]
      fit_seed <- ds_seed * 10L + j
      sd4  <- fit_bf_sd_rs(templates$logit_factor, dat, d, s, seed = fit_seed, n_sampling = 1000)
      sd20 <- fit_bf_sd_rs(templates$logit_factor, dat, d, s, seed = fit_seed, n_sampling = 5000)
      bs   <- fit_bf_bs_rs(tmpl_bs$h1, dat, d, s, ml_h0 = h0$ml_h0,
                           seed = fit_seed, n_sampling = 5000)
      row <- data.frame(true_b1 = b1, rep = k, dist_b1 = d, scale_b1 = s,
                        lbf_bs = log10(bs[["bf10_bridge"]]),
                        bs_pct_err = bs[["pct_error_h1"]],
                        bs_wall_h1 = bs[["wall_s"]], bs_wall_h0 = h0$wall_h0,
                        lbf_kde_4k = log10(sd4[["bf10_kde"]]),
                        lbf_ls_4k = log10(sd4[["bf10_logspline"]]),
                        lbf_norm_4k = log10(sd4[["bf10_normal"]]),
                        lbf_kde_20k = log10(sd20[["bf10_kde"]]),
                        lbf_ls_20k = log10(sd20[["bf10_logspline"]]),
                        lbf_norm_20k = log10(sd20[["bf10_normal"]]),
                        lbf_ls_bsdraws = log10(bs[["bf10_logspline"]]),
                        wall_4k = sd4[["wall_s"]], wall_20k = sd20[["wall_s"]],
                        div_4k = sd4[["n_div"]], div_20k = sd20[["n_div"]],
                        div_bs = bs[["n_div"]])
      print(row, digits = 4)
      rows[[length(rows) + 1]] <- row
    }
  }
  bench <- do.call(rbind, rows)
  saveRDS(bench, bench_file)
  ok  <- is.finite(bench$lbf_bs)
  big <- abs(bench$lbf_bs) >= 1
  for (m in c("kde_4k", "ls_4k", "norm_4k", "kde_20k", "ls_20k", "norm_20k", "ls_bsdraws")) {
    dd <- bench[[paste0("lbf_", m)]] - bench$lbf_bs
    cat(sprintf(paste0("BENCH %-10s vs bridge: all mean|d|=%.3f max|d|=%.3f bias=%+.3f | ",
                       "|log10BF|<1 (n=%d) mean|d|=%.3f max=%.3f | >=1 (n=%d) mean|d|=%.3f max=%.3f\n"),
                m, mean(abs(dd[ok])), max(abs(dd[ok])), mean(dd[ok]),
                sum(ok & !big), mean(abs(dd[ok & !big])), max(abs(dd[ok & !big])),
                sum(ok & big), mean(abs(dd[ok & big])), max(abs(dd[ok & big]))))
  }
  cat(sprintf("BENCH wall: SD 4k mean=%.2fs, SD 20k mean=%.2fs, bridge H1 mean=%.2fs, bridge H0 mean=%.2fs\n",
              mean(bench$wall_4k), mean(bench$wall_20k), mean(bench$bs_wall_h1),
              mean(unique(bench$bs_wall_h0))))
  cat(sprintf("BENCH bridge pct error H1 (100 x CV): median=%.3f%% max=%.3f%%; bridge NA=%d; divergences SD4k/SD20k/bridge = %d/%d/%d\n",
              median(bench$bs_pct_err, na.rm = TRUE), max(bench$bs_pct_err, na.rm = TRUE),
              sum(!ok), sum(bench$div_4k), sum(bench$div_20k), sum(bench$div_bs)))
  message("BENCHMARK COMPLETE")

} else if (run_mode == "full") {
  # Rows: PB_ROWS (or 2nd argument) > PB_CHUNK / SLURM_ARRAY_TASK_ID with
  # PB_N_CHUNKS > all rows
  rows_spec <- if (length(args) > 1) args[2] else Sys.getenv("PB_ROWS", "")
  chunk_id  <- Sys.getenv("PB_CHUNK", Sys.getenv("SLURM_ARRAY_TASK_ID", ""))
  n_chunks  <- as.integer(Sys.getenv("PB_N_CHUNKS", "1"))
  rows <- if (nzchar(rows_spec)) {
    parse_rows_spec(rows_spec, nrow(Design))
  } else if (nzchar(chunk_id)) {
    chunk_rows(as.integer(chunk_id), n_chunks, nrow(Design))
  } else {
    seq_len(nrow(Design))
  }
  n_reps <- as.integer(Sys.getenv("PB_REPS", "200"))
  ncores <- resolve_ncores()
  message("Design rows in full grid: ", nrow(Design), " (",
          paste(names(table(Design$arm)), table(Design$arm), sep = " = ", collapse = ", "),
          "); this run: ", length(rows), " rows; replications: ", n_reps,
          "; global seed: ", global_seed)
  out_root <- Sys.getenv("PB_OUT_ROOT", file.path("output", "Simulation_BFCalibration_RS"))
  if (!grepl("^/", out_root)) out_root <- here(out_root)
  run_chunk(rows, n_reps, out_root, fixed_cfg, ncores)
}
