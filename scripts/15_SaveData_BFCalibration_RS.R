# ==============================================================================
# Script 15: Post-process the random-slopes BF calibration (after Script 14)
# ==============================================================================
# Input : output/Simulation_BFCalibration_RS/res_bf_calibration_rs_raw.rds
#         (Script 14; one row per replication, one column per
#         <stat>_<family>_<scale>). If it does not exist the script falls back,
#         with a message, to output/res_bf_calibration_rs_smoke_bridge.rds
#         (1 replication of one cell per arm, 3 priors per arm) and writes
#         every output with the suffix "_smoke" so test output is never taken
#         for the real run.
#
# Conventions: the second prior parameter is a SCALE (Logistic SD = scale *
# pi / sqrt(3)). BF10 is the bridge-sampling estimate; the logspline
# Savage-Dickey value from the same draws is a cross-check only. A failed fit
# is NA and counted in n_fail; it is never mapped to a number. Inf (extreme
# evidence) is counted separately (n_inf) and counts as BF10 > threshold.
# Arms are never pooled: scales of the logit, probit and continuous arms are
# not comparable.
#
# Outputs (output/<name>[_smoke].rds; all tibbles):
#
# (a) bf_calibration_rs_rates.rds — one row per arm x cell x prior
#   arm, true_b1, n_subjects, n_trials, true_sd_slope  design cell
#   family, scale, prior_label                          effect prior
#   n_reps          replications in the cell
#   n_fail          failed fits of this prior (bf10 NA)
#   n_valid         n_reps - n_fail (denominator of all rates)
#   n_inf           bf10 = Inf among valid fits; n_zero: bf10 = 0
#   k_bf10_gt3, P_bf10_gt3, P_bf10_gt3_lo/_hi   count, rate, Wilson 95% CI
#   k_bf01_gt3, P_bf01_gt3, ...                  BF01 = 1/BF10 > 3
#   k_bf10_gt10, P_bf10_gt10, ...; k_bf01_gt10, P_bf01_gt10, ...
#   n_finite        finite log10 BF10 values used for the quantiles
#   median_log10_bf10, q25_log10_bf10, q75_log10_bf10, iqr_log10_bf10
#   median_log10_bf10_logspline      cross-check estimator, same draws
#   max_abs_dlog10_bridge_logspline  largest |log10 bridge - log10 logspline|
#   max_rhat, min_ess_bulk, min_ess_tail, share_any_div, max_pct_error_h1,
#   share_pct_error_h1_gt5           MCMC / bridge diagnostics over valid fits
#
# (b) bf_calibration_rs_paired.rds — one row per arm x cell x prior pair
#   arm, true_b1, n_subjects, n_trials, true_sd_slope
#   comparison      "vs_cauchy_0.707" (prior A vs the reference Cauchy(0, 0.707)),
#                   "logistic_vs_cauchy_same_scale", "logistic_vs_normal_same_scale",
#                   "logistic_vs_normal_equal_iqr", "logistic_vs_cauchy_equal_iqr"
#   family_a, scale_a, label_a; family_b, scale_b, label_b   A minus B
#   scale_b_implied  for equal-IQR pairs: the exact scale with the same IQR as
#                    Logistic(0, scale_a); scale_b is the nearest grid scale
#   iqr_rel_mismatch |scale_b - scale_b_implied| / scale_b_implied (0 otherwise)
#   n_reps, n_pairs  replications; replications where both fits succeeded
#   P_bf10_gt3_a, P_bf10_gt3_b, diff_P_bf10_gt3, diff_P_bf10_gt3_lo/_hi
#   P_bf01_gt3_a, P_bf01_gt3_b, diff_P_bf01_gt3, diff_P_bf01_gt3_lo/_hi
#   disc_bf10_a_only, disc_bf10_b_only, mcnemar_p_bf10 (exact binomial),
#   disc_bf01_a_only, disc_bf01_b_only, mcnemar_p_bf01
#   n_finite_pairs   pairs with finite log10 BF10 under both priors
#   dlog10_median, dlog10_median_lo/_hi, dlog10_q25, dlog10_q75
#   Bootstrap CIs: replication ids resampled within the cell, n_boot = 1000,
#   percentile 2.5/97.5; seed = 20261104 + design_row (withr::with_seed).
#
# (c) bf_calibration_rs_roc_auc.rds — one row per arm x cell (true_b1 != 0) x prior
#   arm, true_b1, n_subjects, n_trials, true_sd_slope, family, scale, prior_label
#   n_signal, n_noise      valid (non-NA) BF10 values in the effect / null cell
#   n_fail_signal, n_fail_noise
#   auc, auc_lo, auc_hi    Wilcoxon AUC of log10 BF10 (effect vs null) with
#                          500-resample bootstrap CI (independent resampling of
#                          the two cells; seed = 20261104 + design_row).
#                          AUC is identical for H1 and H0 detection.
#     bf_calibration_rs_roc_curves.rds — long ROC object
#   arm, true_b1, n_subjects, n_trials, true_sd_slope, family, scale, prior_label
#   direction   "H1" (score BF10, signal = effect cell, noise = null cell) or
#               "H0" (score BF01, signal = null cell, noise = effect cell)
#   kind        "curve" (fpr on a 101-point grid, tpr interpolated from the
#               empirical ROC; tpr_lo/tpr_hi 500-resample bootstrap band) or
#               "operating_point" (threshold_bf = 3 or 10: fpr = P(score > t | noise),
#               tpr = P(score > t | signal); tpr_lo/tpr_hi NA)
#   threshold_bf, fpr, tpr, tpr_lo, tpr_hi
#
# (d) bf_calibration_rs_sd_decomp.rds — one row per arm x cell
#   Savage-Dickey identity log BF10 = log p(0) - log p(0 | y) (natural log):
#   the prior ordinate is analytic, the posterior ordinate is implied by the
#   bridge BF, log_post_0 = log_prior_0 - log bf10_bridge. No |log BF| cutoff
#   is applied: the bridge benchmark (output/bf_rs_benchmark.rds) showed a
#   bridge CV of 1.2-2.4% over log10 BF10 in [-1.3, 4.3], so unlike the KDE
#   estimates of Script 11 the bridge values are not unreliable in the tail.
#   Fits with pct_error_h1 > 5 are flagged (counted), not excluded.
#   n_reps, n_reps_all_finite   replications with finite BF10 under ALL priors
#   n_priors
#   iqr_log_post_0_median/_q25/_q75   per replication, IQR of log_post_0 across
#                                     the priors of the arm; distribution over reps
#   range_log_post_0_median           per-replication range across priors
#   iqr_log_bf_median                 per-replication IQR of log BF10 across priors
#   spearman_log_bf_min/_median/_max  Spearman correlation of log BF10 across
#                                     replications, over all prior pairs (NA if < 3 reps)
#   n_flag_pct_error_gt5, n_flag_logspline_dlog10_gt0.5   flagged fits (not excluded)
#     bf_calibration_rs_sd_decomp_reps.rds — one row per replication x prior
#   arm, true_b1, n_subjects, n_trials, true_sd_slope, design_row, rep,
#   family, scale, prior_label, log_prior_0, log_bf10, log10_bf10, log_post_0,
#   pct_error_h1, flag_pct_error_gt5, flag_logspline_dlog10_gt0.5
#
# Every object carries attr "input_source" ("merged" or "smoke_bridge").
# Test hook: PB_RS_RAW=<file> reads any raw-format file instead (a SimDesign
# object or the data frame of Script 14) and tags outputs with
# "_" + PB_RS_TAG (default "custom").
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(here)
})

source(here("R", "bf_helpers.R"))
source(here("R", "rs_postprocess.R"))

check_equal_iqr_formulas()
check_roc_helpers()
message("Formula and helper self-checks passed")

n_boot_paired <- 1000
n_boot_roc    <- 500
boot_seed     <- 20261104
iqr_rel_tol   <- 0.03   # judgment call: equal-IQR partner must be a grid scale
                        # within 3% of the implied scale (see message below)
cell_vars <- c("arm", "true_b1", "n_subjects", "n_trials", "true_sd_slope")

# ---- Input -------------------------------------------------------------------

raw_file <- here("output", "Simulation_BFCalibration_RS", "res_bf_calibration_rs_raw.rds")
if (nzchar(Sys.getenv("PB_RS_RAW"))) {
  # Test hook: any raw-format file (e.g. a synthetic one); outputs are tagged
  raw <- readRDS(Sys.getenv("PB_RS_RAW"))
  if (inherits(raw, "SimDesign")) raw <- as.data.frame(SimDesign::SimResults(raw))
  input_source <- paste0("custom:", Sys.getenv("PB_RS_RAW"))
  out_tag <- paste0("_", Sys.getenv("PB_RS_TAG", "custom"))
  message("Input overridden by PB_RS_RAW; outputs get the suffix '", out_tag, "'")
} else if (file.exists(raw_file)) {
  raw <- readRDS(raw_file)
  input_source <- "merged"
  out_tag <- ""
} else {
  smoke_file <- here("output", "res_bf_calibration_rs_smoke_bridge.rds")
  message("Merged results not found (", raw_file, "); falling back to the ",
          "smoke-test object ", smoke_file, ". Outputs get the suffix '_smoke'.")
  raw <- as.data.frame(SimDesign::SimResults(readRDS(smoke_file)))
  input_source <- "smoke_bridge"
  out_tag <- "_smoke"
}
out_path <- function(name) here("output", paste0(name, out_tag, ".rds"))
save_out <- function(x, name) {
  attr(x, "input_source") <- input_source
  saveRDS(x, out_path(name))
  message("Saved ", out_path(name), " (", nrow(x), " rows)")
}

long <- rs_raw_to_long(raw)
message("Long format: ", nrow(long), " replication x prior rows; ",
        n_distinct(long$design_row), " design rows; failed fits: ",
        sum(is.na(long$bf10_bridge)))

# ---- (a) Rates with Wilson CIs ------------------------------------------------

rate_cols <- function(k, n, name) {
  ci <- t(mapply(wilson_ci, k, n))
  tibble::tibble(!!paste0("k_", name) := k,
                 !!paste0("P_", name) := ifelse(n > 0, k / n, NA_real_),
                 !!paste0("P_", name, "_lo") := ci[, 1],
                 !!paste0("P_", name, "_hi") := ci[, 2])
}
safe_q <- function(x, p) if (length(x)) unname(stats::quantile(x, p)) else NA_real_
safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

rates <- long |>
  summarise(
    n_reps  = n(),
    n_fail  = sum(is.na(bf10_bridge)),
    n_valid = sum(!is.na(bf10_bridge)),
    n_inf   = sum(is.infinite(bf10_bridge)),
    n_zero  = sum(bf10_bridge == 0, na.rm = TRUE),
    k_bf10_gt3  = sum(bf10_bridge > 3,       na.rm = TRUE),
    k_bf01_gt3  = sum(1 / bf10_bridge > 3,   na.rm = TRUE),
    k_bf10_gt10 = sum(bf10_bridge > 10,      na.rm = TRUE),
    k_bf01_gt10 = sum(1 / bf10_bridge > 10,  na.rm = TRUE),
    n_finite          = sum(is.finite(log10_bf10)),
    median_log10_bf10 = safe_q(log10_bf10[is.finite(log10_bf10)], 0.5),
    q25_log10_bf10    = safe_q(log10_bf10[is.finite(log10_bf10)], 0.25),
    q75_log10_bf10    = safe_q(log10_bf10[is.finite(log10_bf10)], 0.75),
    median_log10_bf10_logspline =
      safe_q(log10(bf10_logspline)[is.finite(log10(bf10_logspline))], 0.5),
    max_abs_dlog10_bridge_logspline = safe_max(abs(log10_bf10 - log10(bf10_logspline))),
    max_rhat      = safe_max(rhat_max),
    min_ess_bulk  = safe_min(ess_bulk_min),
    min_ess_tail  = safe_min(ess_tail_min),
    share_any_div = safe_mean(n_div > 0),
    max_pct_error_h1       = safe_max(pct_error_h1),
    share_pct_error_h1_gt5 = safe_mean(pct_error_h1 > 5),
    .by = all_of(c(cell_vars, "family", "scale", "prior_label"))
  ) |>
  mutate(iqr_log10_bf10 = q75_log10_bf10 - q25_log10_bf10)

rates <- bind_cols(
  select(rates, -starts_with("k_")),
  rate_cols(rates$k_bf10_gt3,  rates$n_valid, "bf10_gt3"),
  rate_cols(rates$k_bf01_gt3,  rates$n_valid, "bf01_gt3"),
  rate_cols(rates$k_bf10_gt10, rates$n_valid, "bf10_gt10"),
  rate_cols(rates$k_bf01_gt10, rates$n_valid, "bf01_gt10")
) |>
  relocate(n_finite, starts_with("median_"), starts_with("q25_"), starts_with("q75_"),
           iqr_log10_bf10, .after = P_bf01_gt10_hi) |>
  arrange(arm, true_sd_slope, n_trials, n_subjects, true_b1, family, scale)
save_out(rates, "bf_calibration_rs_rates")

# ---- (b) Paired differences (same dataset) ------------------------------------

prior_grids <- make_prior_grids_rs()

build_pairs <- function(grid) {
  g <- tibble::tibble(family = as.character(grid$dist_b1), scale = grid$scale_b1)
  has <- function(f, s) any(g$family == f & abs(g$scale - s) < 1e-9)
  pairs <- list()
  add <- function(comparison, fa, sa, fb, sb, implied = NA_real_) {
    pairs[[length(pairs) + 1]] <<- tibble::tibble(
      comparison = comparison, family_a = fa, scale_a = sa, family_b = fb, scale_b = sb,
      scale_b_implied = implied,
      iqr_rel_mismatch = if (is.na(implied)) 0 else abs(sb - implied) / implied)
  }
  if (has("cauchy", 0.707)) {
    for (i in seq_len(nrow(g))) {
      if (g$family[i] == "cauchy" && abs(g$scale[i] - 0.707) < 1e-9) next
      add("vs_cauchy_0.707", g$family[i], g$scale[i], "cauchy", 0.707)
    }
  }
  for (s in g$scale[g$family == "logistic"]) {
    if (has("cauchy", s)) add("logistic_vs_cauchy_same_scale", "logistic", s, "cauchy", s)
    if (has("normal", s)) add("logistic_vs_normal_same_scale", "logistic", s, "normal", s)
    for (fb in c("normal", "cauchy")) {
      implied <- equal_iqr_scale(s, fb)
      cand <- g$scale[g$family == fb]
      if (!length(cand)) next
      near <- cand[which.min(abs(cand - implied))]
      if (abs(near - implied) / implied <= iqr_rel_tol) {
        add(paste0("logistic_vs_", fb, "_equal_iqr"), "logistic", s, fb, near, implied)
      }
    }
  }
  if (!length(pairs)) return(NULL)
  list_rbind(pairs) |>
    mutate(label_a = prior_label_rs(family_a, scale_a),
           label_b = prior_label_rs(family_b, scale_b))
}
pairs_by_arm <- imap(prior_grids, function(g, a) {
  p <- build_pairs(g)
  if (!is.null(p)) p$arm <- a
  p
})
pairs_all <- list_rbind(pairs_by_arm)

message("Equal-IQR pairs present in the grids (tolerance ", iqr_rel_tol * 100, "%):")
eq <- filter(pairs_all, grepl("equal_iqr", comparison))
if (nrow(eq)) {
  for (i in seq_len(nrow(eq))) message(sprintf(
    "  %-16s %s vs %s (implied scale %.4f, mismatch %.1f%%)", eq$arm[i], eq$label_a[i],
    eq$label_b[i], eq$scale_b_implied[i], 100 * eq$iqr_rel_mismatch[i]))
} else message("  none")
for (a in names(prior_grids)) {
  ls <- prior_grids[[a]]$scale_b1[prior_grids[[a]]$dist_b1 == "logistic"]
  for (fb in c("normal", "cauchy")) {
    missing <- ls[!ls %in% eq$scale_a[eq$arm == a & eq$family_b == fb]]
    if (length(missing)) message(sprintf(
      "  %-16s no %s grid scale within tolerance for Logistic scale(s) %s (implied: %s)",
      a, fb, paste(missing, collapse = ", "),
      paste(round(equal_iqr_scale(missing, fb), 3), collapse = ", ")))
  }
}

boot_q <- function(x) unname(stats::quantile(x, c(0.025, 0.975), na.rm = TRUE))
mcnemar_p <- function(b, c) if (b + c == 0) NA_real_ else stats::binom.test(b, b + c, 0.5)$p.value

paired_cell <- function(df, pairs, design_row) {
  wide <- df |>
    select(rep, prior_label, bf10_bridge) |>
    pivot_wider(names_from = prior_label, values_from = bf10_bridge) |>
    arrange(rep)
  n_reps <- nrow(wide)
  idx <- withr::with_seed(as.integer(boot_seed + design_row),
                          matrix(sample.int(n_reps, n_boot_paired * n_reps, replace = TRUE),
                                 nrow = n_boot_paired))
  boot_stat <- function(d, fun) {
    # d has NA where a pair is unusable; a resample with no usable pair gives NA
    m <- matrix(d[idx], nrow = n_boot_paired)
    boot_q(apply(m, 1, fun))
  }
  pmap(pairs, function(comparison, family_a, scale_a, family_b, scale_b,
                       scale_b_implied, iqr_rel_mismatch, label_a, label_b, arm) {
    if (!all(c(label_a, label_b) %in% names(wide))) return(NULL)
    a <- wide[[label_a]]; b <- wide[[label_b]]
    ok <- !is.na(a) & !is.na(b)
    h1a <- a > 3; h1b <- b > 3; h0a <- 1 / a > 3; h0b <- 1 / b > 3
    d_h1 <- ifelse(ok, as.numeric(h1a) - as.numeric(h1b), NA_real_)
    d_h0 <- ifelse(ok, as.numeric(h0a) - as.numeric(h0b), NA_real_)
    fin  <- ok & is.finite(log10(a)) & is.finite(log10(b))
    d_lg <- ifelse(fin, log10(a) - log10(b), NA_real_)
    ci_h1 <- boot_stat(d_h1, function(x) mean(x, na.rm = TRUE))
    ci_h0 <- boot_stat(d_h0, function(x) mean(x, na.rm = TRUE))
    ci_lg <- boot_stat(d_lg, function(x) stats::median(x, na.rm = TRUE))
    tibble::tibble(
      comparison = comparison, family_a = family_a, scale_a = scale_a, label_a = label_a,
      family_b = family_b, scale_b = scale_b, label_b = label_b,
      scale_b_implied = scale_b_implied, iqr_rel_mismatch = iqr_rel_mismatch,
      n_reps = n_reps, n_pairs = sum(ok),
      P_bf10_gt3_a = safe_mean(h1a[ok]), P_bf10_gt3_b = safe_mean(h1b[ok]),
      diff_P_bf10_gt3 = safe_mean(d_h1), diff_P_bf10_gt3_lo = ci_h1[1], diff_P_bf10_gt3_hi = ci_h1[2],
      P_bf01_gt3_a = safe_mean(h0a[ok]), P_bf01_gt3_b = safe_mean(h0b[ok]),
      diff_P_bf01_gt3 = safe_mean(d_h0), diff_P_bf01_gt3_lo = ci_h0[1], diff_P_bf01_gt3_hi = ci_h0[2],
      disc_bf10_a_only = sum(ok & h1a & !h1b), disc_bf10_b_only = sum(ok & !h1a & h1b),
      mcnemar_p_bf10 = mcnemar_p(sum(ok & h1a & !h1b), sum(ok & !h1a & h1b)),
      disc_bf01_a_only = sum(ok & h0a & !h0b), disc_bf01_b_only = sum(ok & !h0a & h0b),
      mcnemar_p_bf01 = mcnemar_p(sum(ok & h0a & !h0b), sum(ok & !h0a & h0b)),
      n_finite_pairs = sum(fin),
      dlog10_median = safe_q(d_lg[fin], 0.5), dlog10_median_lo = ci_lg[1], dlog10_median_hi = ci_lg[2],
      dlog10_q25 = safe_q(d_lg[fin], 0.25), dlog10_q75 = safe_q(d_lg[fin], 0.75))
  }) |> list_rbind()
}

paired <- long |>
  nest(data = -all_of(c(cell_vars, "design_row"))) |>
  mutate(res = pmap(list(data, arm, design_row), function(d, a, dr) {
    p <- pairs_by_arm[[a]]
    if (is.null(p)) NULL else paired_cell(d, p, dr)
  })) |>
  select(-data, -design_row) |>
  unnest(res) |>
  arrange(arm, true_sd_slope, n_trials, n_subjects, true_b1, comparison, family_a, scale_a)
save_out(paired, "bf_calibration_rs_paired")

# ---- (c) ROC / AUC: effect cell vs matching null cell -------------------------

fpr_grid <- seq(0, 1, length.out = 101)
null_vars <- c("arm", "n_subjects", "n_trials", "true_sd_slope", "family", "scale")

effect_cells <- long |> filter(true_b1 != 0) |>
  nest(sig = -all_of(c(cell_vars, "design_row", "family", "scale", "prior_label")))
null_cells <- long |> filter(true_b1 == 0) |>
  nest(noi = -all_of(c(null_vars, "prior_label")))
roc_input <- inner_join(effect_cells, select(null_cells, -prior_label),
                        by = join_by(arm, n_subjects, n_trials, true_sd_slope, family, scale))
message("ROC/AUC: ", nrow(roc_input), " effect-cell x prior combinations with a ",
        "matching null cell (", nrow(effect_cells) - nrow(roc_input), " without)")

roc_one <- function(sig, noi, design_row) {
  s_all <- sig$log10_bf10; n_all <- noi$log10_bf10
  s <- s_all[!is.na(s_all)]; n <- n_all[!is.na(n_all)]
  seed <- as.integer(boot_seed + design_row)
  boot_pairs <- withr::with_seed(seed, replicate(n_boot_roc, list(
    s = s[sample.int(length(s), replace = TRUE)],
    n = n[sample.int(length(n), replace = TRUE)]), simplify = FALSE))
  auc <- auc_wilcoxon(s, n)
  auc_b <- map_dbl(boot_pairs, function(b) auc_wilcoxon(b$s, b$n))
  auc_row <- tibble::tibble(
    n_signal = length(s), n_noise = length(n),
    n_fail_signal = sum(is.na(s_all)), n_fail_noise = sum(is.na(n_all)),
    auc = auc, auc_lo = boot_q(auc_b)[1], auc_hi = boot_q(auc_b)[2])
  curve_dir <- function(direction, signal, noise, boots) {
    if (!length(signal) || !length(noise)) return(NULL)
    tpr <- roc_at_grid(noise, signal, fpr_grid)
    bm  <- vapply(boots, function(b) roc_at_grid(b$noise, b$signal, fpr_grid),
                  numeric(length(fpr_grid)))
    curve <- tibble::tibble(direction = direction, kind = "curve", threshold_bf = NA_real_,
                            fpr = fpr_grid, tpr = tpr,
                            tpr_lo = apply(bm, 1, stats::quantile, 0.025, na.rm = TRUE),
                            tpr_hi = apply(bm, 1, stats::quantile, 0.975, na.rm = TRUE))
    ops <- map(c(3, 10), function(t) tibble::tibble(
      direction = direction, kind = "operating_point", threshold_bf = t,
      fpr = mean(noise > log10(t)), tpr = mean(signal > log10(t)),
      tpr_lo = NA_real_, tpr_hi = NA_real_)) |> list_rbind()
    bind_rows(curve, ops)
  }
  curves <- bind_rows(
    curve_dir("H1", s, n, map(boot_pairs, \(b) list(signal = b$s, noise = b$n))),
    curve_dir("H0", -n, -s, map(boot_pairs, \(b) list(signal = -b$n, noise = -b$s))))
  list(auc = auc_row, curves = curves)
}

# Zero-row templates fix the column contract even when no cell has a null partner
key_template <- tibble::tibble(arm = character(), true_b1 = numeric(), n_subjects = numeric(),
                               n_trials = numeric(), true_sd_slope = numeric(),
                               family = character(), scale = numeric(), prior_label = character())
auc_template <- bind_cols(key_template, tibble::tibble(
  n_signal = integer(), n_noise = integer(), n_fail_signal = integer(), n_fail_noise = integer(),
  auc = numeric(), auc_lo = numeric(), auc_hi = numeric()))
curves_template <- bind_cols(key_template, tibble::tibble(
  direction = character(), kind = character(), threshold_bf = numeric(), fpr = numeric(),
  tpr = numeric(), tpr_lo = numeric(), tpr_hi = numeric()))

roc_res <- roc_input |>
  mutate(res = pmap(list(sig, noi, design_row), roc_one)) |>
  select(-sig, -noi, -design_row)
unnest_part <- function(part, template) {
  if (nrow(roc_res) == 0) return(template)
  roc_res |> mutate(part = map(res, part)) |> select(-res) |> unnest(part)
}
roc_auc <- bind_rows(auc_template, unnest_part("auc", auc_template)) |>
  arrange(arm, true_sd_slope, n_trials, n_subjects, true_b1, family, scale)
roc_curves <- bind_rows(curves_template, unnest_part("curves", curves_template)) |>
  arrange(arm, true_sd_slope, n_trials, n_subjects, true_b1, family, scale, direction, kind, fpr)
stopifnot(identical(names(roc_auc), names(auc_template)),
          identical(names(roc_curves), names(curves_template)))
save_out(roc_auc, "bf_calibration_rs_roc_auc")
save_out(roc_curves, "bf_calibration_rs_roc_curves")

# ---- (d) Savage-Dickey decomposition ------------------------------------------

decomp_reps <- long |>
  mutate(
    log_prior_0 = log(map2_dbl(family, scale, prior_density_at_0)),
    log_bf10    = log(bf10_bridge),
    log_post_0  = log_prior_0 - log_bf10,
    flag_pct_error_gt5 = !is.na(pct_error_h1) & pct_error_h1 > 5,
    flag_logspline_dlog10_gt0.5 =
      is.finite(log10_bf10) & is.finite(log10(bf10_logspline)) &
      abs(log10_bf10 - log10(bf10_logspline)) > 0.5
  ) |>
  select(all_of(c(cell_vars, "design_row", "rep", "family", "scale", "prior_label")),
         log_prior_0, log_bf10, log10_bf10, log_post_0, pct_error_h1,
         flag_pct_error_gt5, flag_logspline_dlog10_gt0.5)
stopifnot(isTRUE(all.equal(decomp_reps$log_bf10,
                           decomp_reps$log_prior_0 - decomp_reps$log_post_0)))

per_rep <- decomp_reps |>
  summarise(n_priors = n(),
            all_finite = all(is.finite(log_bf10)),
            iqr_log_post_0   = if (all(is.finite(log_bf10))) IQR(log_post_0) else NA_real_,
            range_log_post_0 = if (all(is.finite(log_bf10))) diff(range(log_post_0)) else NA_real_,
            iqr_log_bf       = if (all(is.finite(log_bf10))) IQR(log_bf10) else NA_real_,
            .by = all_of(c(cell_vars, "design_row", "rep")))

spearman_cell <- function(df) {
  w <- df |> select(rep, prior_label, log_bf10) |>
    pivot_wider(names_from = prior_label, values_from = log_bf10) |>
    select(-rep) |> as.matrix()
  w[!is.finite(w)] <- NA
  if (nrow(w) < 3 || ncol(w) < 2) return(c(NA_real_, NA_real_, NA_real_))
  r <- suppressWarnings(stats::cor(w, method = "spearman", use = "pairwise.complete.obs"))
  v <- r[upper.tri(r)]
  if (all(is.na(v))) return(c(NA_real_, NA_real_, NA_real_))
  c(min(v, na.rm = TRUE), stats::median(v, na.rm = TRUE), max(v, na.rm = TRUE))
}

decomp <- per_rep |>
  summarise(n_reps = n(), n_reps_all_finite = sum(all_finite), n_priors = max(n_priors),
            iqr_log_post_0_median = safe_q(iqr_log_post_0[all_finite], 0.5),
            iqr_log_post_0_q25    = safe_q(iqr_log_post_0[all_finite], 0.25),
            iqr_log_post_0_q75    = safe_q(iqr_log_post_0[all_finite], 0.75),
            range_log_post_0_median = safe_q(range_log_post_0[all_finite], 0.5),
            iqr_log_bf_median       = safe_q(iqr_log_bf[all_finite], 0.5),
            .by = all_of(c(cell_vars, "design_row"))) |>
  left_join(
    decomp_reps |> nest(d = -all_of(c(cell_vars, "design_row"))) |>
      mutate(sp = map(d, spearman_cell),
             spearman_log_bf_min = map_dbl(sp, 1), spearman_log_bf_median = map_dbl(sp, 2),
             spearman_log_bf_max = map_dbl(sp, 3)) |> select(-d, -sp),
    by = join_by(arm, true_b1, n_subjects, n_trials, true_sd_slope, design_row)) |>
  left_join(
    decomp_reps |> summarise(n_flag_pct_error_gt5 = sum(flag_pct_error_gt5),
                             n_flag_logspline_dlog10_gt0.5 = sum(flag_logspline_dlog10_gt0.5),
                             .by = all_of(c(cell_vars, "design_row"))),
    by = join_by(arm, true_b1, n_subjects, n_trials, true_sd_slope, design_row)) |>
  select(-design_row) |>
  arrange(arm, true_sd_slope, n_trials, n_subjects, true_b1)
save_out(decomp, "bf_calibration_rs_sd_decomp")
save_out(decomp_reps, "bf_calibration_rs_sd_decomp_reps")
message("Decomposition: no |log BF| cutoff applied (bridge estimates; see header). ",
        "Flagged fits: pct_error_h1 > 5: ", sum(decomp_reps$flag_pct_error_gt5),
        "; |log10 bridge - logspline| > 0.5: ", sum(decomp_reps$flag_logspline_dlog10_gt0.5))

# ---- Console summary -------------------------------------------------------------

options(width = 200, pillar.sigfig = 3)
message("\n== Rates (input: ", input_source, ") ==")
print(select(rates, all_of(cell_vars), prior_label, n_reps, n_fail, n_inf,
             P_bf10_gt3, P_bf10_gt3_lo, P_bf10_gt3_hi, P_bf01_gt3, median_log10_bf10), n = 30)
message("\n== Paired differences (first rows) ==")
print(select(paired, arm, true_b1, n_subjects, comparison, label_a, label_b, n_pairs,
             diff_P_bf10_gt3, diff_P_bf10_gt3_lo, diff_P_bf10_gt3_hi, dlog10_median), n = 30)
message("\n== AUC (first rows) ==")
print(select(roc_auc, all_of(cell_vars), prior_label, n_signal, n_noise, auc, auc_lo, auc_hi), n = 20)
message("\n== Savage-Dickey decomposition ==")
print(select(decomp, all_of(cell_vars), n_reps, n_reps_all_finite, n_priors,
             iqr_log_post_0_median, iqr_log_bf_median, spearman_log_bf_median), n = 30)
message("DONE: script 15 (input: ", input_source, ", output suffix: '", out_tag, "')")
