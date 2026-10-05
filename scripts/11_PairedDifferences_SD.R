# ==============================================================================
# Script 11: Paired family differences and Savage-Dickey decomposition (T1)
# ==============================================================================
# Purpose (Reviewer 1, revision R1):
#   1. Count BF10 values at the 1e30 cap per cell and prior in the JZS v3 run.
#   2. Paired (same-dataset) differences between Logistic(0, s) and Cauchy(0, s)
#      in P(BF10 > 3) and P(BF01 > 3), with bootstrap CIs and McNemar counts,
#      plus the paired difference in log10 BF10.
#   3. Savage-Dickey decomposition: log10 BF10 = log10 pi(0) - log10 pi(0 | y).
#      Prior ordinates are analytic; the posterior ordinate is recovered as
#      pi(0) / BF10 for each prior and replication.
#
# Pairing note: in scripts/03d the design factor is sd_b1; Logistic and Cauchy
# are fitted to the SAME dataset only within a given sd_b1. Logistic(0, 0.25)
# and Cauchy(0, 0.707) live in different design rows and therefore on different
# datasets. That comparison is reported as unpaired (independent bootstrap).
#
# Input : output/Simulation_BFCalibration_JZS_v3/BFCalib_JZS_v3_Cond-<i> (qs2)
# Output: output/t1_cap_counts.rds
#         output/t1_paired_differences.rds
#         output/t1_unpaired_scale_comparison.rds
#         output/t1_sd_decomposition.rds
#         output/t1_prior_ordinates.rds
# ==============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(here)

cap_value <- 1e30
n_boot    <- 2000
set.seed(20261005)

# --- Load replication-level results ------------------------------------------

results_dir <- here("output", "Simulation_BFCalibration_JZS_v3")
cond_files  <- list.files(results_dir, pattern = "^BFCalib_JZS_v3_Cond-\\d+$",
                          full.names = TRUE)
stopifnot(length(cond_files) == 36)

read_cond <- function(path) {
  cond <- qs2::qd_read(path)
  as_tibble(cond$condition) |>
    select(sd_b1, true_b1, n_subjects, n_trials) |>
    slice(rep(1L, nrow(cond$results))) |>
    bind_cols(as_tibble(cond$results)) |>
    mutate(rep = row_number(), cond_file = basename(path))
}

reps <- map(cond_files, read_cond) |> list_rbind()

# Cross-check against the post-processed long file used by the manuscript.
# 04d writes, per condition i = 1..36 in order, one logistic and one cauchy row
# per replication. Compare element-wise after ordering by condition number.
raw_long <- readRDS(here("output", "bf_calibration_jzs_raw.rds"))
reps_ordered <- reps |>
  mutate(cond_id = as.integer(sub(".*Cond-", "", cond_file))) |>
  arrange(cond_id, rep)
stopifnot(
  "raw file row count does not match v3 condition files" =
    nrow(raw_long) == 2L * nrow(reps),
  "raw file logistic BF10 differs from v3 condition files" =
    isTRUE(all.equal(raw_long$BF10[raw_long$prior == "logistic"], reps_ordered$BF10_logistic)),
  "raw file cauchy BF10 differs from v3 condition files" =
    isTRUE(all.equal(raw_long$BF10[raw_long$prior == "cauchy"], reps_ordered$BF10_cauchy)),
  "raw file design columns differ from v3 condition files" =
    isTRUE(all.equal(raw_long$true_b1[raw_long$prior == "logistic"], reps_ordered$true_b1)) &&
    isTRUE(all.equal(raw_long$sd_b1[raw_long$prior == "logistic"], reps_ordered$sd_b1))
)
message("Cross-check vs bf_calibration_jzs_raw.rds passed (element-wise, ",
        nrow(raw_long), " rows)")

cell_vars <- c("sd_b1", "true_b1", "n_subjects", "n_trials")

# --- Part 1: values at the 1e30 cap ------------------------------------------
# 03d (line 84) maps BOTH fit failures (NA) and Inf (KDE underflow) to 1e30, so
# the two are not separable here and NA counts are structurally 0.

cap_counts <- reps |>
  pivot_longer(c(BF10_logistic, BF10_cauchy), names_to = "prior",
               names_prefix = "BF10_", values_to = "BF10") |>
  summarise(
    n_reps = n(),
    n_na_structurally_zero = sum(is.na(BF10)),
    n_at_cap_1e30_underflow_or_fit_failure = sum(BF10 >= cap_value, na.rm = TRUE),
    n_gt_1e6_below_cap = sum(BF10 > 1e6 & BF10 < cap_value, na.rm = TRUE),
    max_below_cap_log10 = max(log10(BF10[BF10 < cap_value])),
    .by = all_of(c(cell_vars, "prior"))
  )
stopifnot(sum(cap_counts$n_na_structurally_zero) == 0)
saveRDS(cap_counts, here("output", "t1_cap_counts.rds"))
message("BF10 capped at 1e30 (underflow or fit failure, not separable): ",
        sum(cap_counts$n_at_cap_1e30_underflow_or_fit_failure),
        " values in ", sum(cap_counts$n_at_cap_1e30_underflow_or_fit_failure > 0),
        " of ", nrow(cap_counts), " cell x prior combinations ",
        "(NA is structurally 0: 03d maps NA to 1e30)")

# --- Part 2: paired family differences ----------------------------------------

boot_ci <- function(x, n_boot) {
  # x: per-replication paired difference (numeric, no NA)
  stopifnot(is.numeric(x), !anyNA(x))
  b <- replicate(n_boot, mean(x[sample.int(length(x), replace = TRUE)]))
  unname(quantile(b, c(0.025, 0.975)))
}

paired_summary <- function(df) {
  h1_l <- df$BF10_logistic > 3
  h1_c <- df$BF10_cauchy   > 3
  h0_l <- 1 / df$BF10_logistic > 3
  h0_c <- 1 / df$BF10_cauchy   > 3
  d_h1 <- as.numeric(h1_l) - as.numeric(h1_c)
  d_h0 <- as.numeric(h0_l) - as.numeric(h0_c)
  ci_h1 <- boot_ci(d_h1, n_boot)
  ci_h0 <- boot_ci(d_h0, n_boot)

  mcnemar_p <- function(b, c) {
    if (b + c == 0) return(NA_real_)
    binom.test(b, b + c, 0.5)$p.value
  }
  b_h1 <- sum(h1_l & !h1_c); c_h1 <- sum(!h1_l & h1_c)
  b_h0 <- sum(h0_l & !h0_c); c_h0 <- sum(!h0_l & h0_c)

  ok    <- df$BF10_logistic < cap_value & df$BF10_cauchy < cap_value
  d_log <- log10(df$BF10_logistic[ok]) - log10(df$BF10_cauchy[ok])

  tibble(
    n_reps = nrow(df),
    P_BF10_gt3_logistic = mean(h1_l), P_BF10_gt3_cauchy = mean(h1_c),
    diff_BF10_gt3 = mean(d_h1), diff_BF10_gt3_lo = ci_h1[1], diff_BF10_gt3_hi = ci_h1[2],
    disc_BF10_L_only = b_h1, disc_BF10_C_only = c_h1,
    mcnemar_p_BF10 = mcnemar_p(b_h1, c_h1),
    P_BF01_gt3_logistic = mean(h0_l), P_BF01_gt3_cauchy = mean(h0_c),
    diff_BF01_gt3 = mean(d_h0), diff_BF01_gt3_lo = ci_h0[1], diff_BF01_gt3_hi = ci_h0[2],
    disc_BF01_L_only = b_h0, disc_BF01_C_only = c_h0,
    mcnemar_p_BF01 = mcnemar_p(b_h0, c_h0),
    n_excluded_cap = sum(!ok),
    dlog10_median = median(d_log),
    dlog10_q25 = unname(quantile(d_log, 0.25)),
    dlog10_q75 = unname(quantile(d_log, 0.75)),
    dlog10_iqr = IQR(d_log)
  )
}

paired <- reps |>
  nest(data = -all_of(cell_vars)) |>
  mutate(res = map(data, paired_summary)) |>
  select(-data) |>
  unnest(res) |>
  arrange(sd_b1, true_b1, n_subjects, n_trials)
saveRDS(paired, here("output", "t1_paired_differences.rds"))

# --- Part 2b: unpaired Logistic(0, .25) vs Cauchy(0, .707) ---------------------
# Different datasets (different design rows): independent bootstrap per arm.

unpaired_scale <- function(log_df, cau_df) {
  h1_l <- log_df$BF10_logistic > 3; h1_c <- cau_df$BF10_cauchy > 3
  h0_l <- 1 / log_df$BF10_logistic > 3; h0_c <- 1 / cau_df$BF10_cauchy > 3
  boot_diff <- function(a, b) {
    d <- replicate(n_boot, mean(sample(a, replace = TRUE)) - mean(sample(b, replace = TRUE)))
    unname(quantile(d, c(0.025, 0.975)))
  }
  ci1 <- boot_diff(h1_l, h1_c); ci0 <- boot_diff(h0_l, h0_c)
  tibble(
    P_BF10_gt3_log025 = mean(h1_l), P_BF10_gt3_cau0707 = mean(h1_c),
    diff_BF10_gt3 = mean(h1_l) - mean(h1_c), diff_BF10_gt3_lo = ci1[1], diff_BF10_gt3_hi = ci1[2],
    P_BF01_gt3_log025 = mean(h0_l), P_BF01_gt3_cau0707 = mean(h0_c),
    diff_BF01_gt3 = mean(h0_l) - mean(h0_c), diff_BF01_gt3_lo = ci0[1], diff_BF01_gt3_hi = ci0[2]
  )
}

unpaired <- reps |>
  filter(sd_b1 == 0.25) |>
  nest(log_df = -c(true_b1, n_subjects, n_trials)) |>
  inner_join(
    reps |> filter(sd_b1 > 0.7) |> nest(cau_df = -c(sd_b1, true_b1, n_subjects, n_trials)) |>
      select(-sd_b1),
    by = join_by(true_b1, n_subjects, n_trials)
  ) |>
  mutate(res = map2(log_df, cau_df, unpaired_scale)) |>
  select(true_b1, n_subjects, n_trials, res) |>
  unnest(res) |>
  arrange(true_b1, n_subjects, n_trials)
saveRDS(unpaired, here("output", "t1_unpaired_scale_comparison.rds"))

# --- Part 3: Savage-Dickey decomposition --------------------------------------

prior_ord <- tribble(
  ~prior,     ~scale,
  "logistic", 0.25,  "cauchy", 0.25,
  "logistic", 0.707, "cauchy", 0.707
) |>
  mutate(ordinate = if_else(prior == "logistic", dlogis(0, 0, scale), dcauchy(0, 0, scale)),
         log10_ordinate = log10(ordinate))
saveRDS(prior_ord, here("output", "t1_prior_ordinates.rds"))

ord <- function(p, s) log10(if (p == "logistic") dlogis(0, 0, s) else dcauchy(0, 0, s))
r1_checks <- c(
  family_shift_025   = ord("cauchy", 0.25)    - ord("logistic", 0.25),
  family_shift_0707  = ord("cauchy", 0.707)   - ord("logistic", 0.707),
  scale_shift_logist = ord("logistic", 0.25)  - ord("logistic", 0.707),
  scale_shift_cauchy = ord("cauchy", 0.25)    - ord("cauchy", 0.707)
)
print(round(r1_checks, 4))

decomp <- reps |>
  mutate(
    log10_prior_L = log10(dlogis(0, 0, sd_b1)),
    log10_prior_C = log10(dcauchy(0, 0, sd_b1)),
    ok = BF10_logistic < cap_value & BF10_cauchy < cap_value,
    log10_post_L = log10_prior_L - log10(BF10_logistic),
    log10_post_C = log10_prior_C - log10(BF10_cauchy),
    d_post = log10_post_L - log10_post_C,          # data-dependent term
    d_prior = log10_prior_L - log10_prior_C,       # constant within scale
    d_bf = log10(BF10_logistic) - log10(BF10_cauchy)
  )

# Savage-Dickey BFs here are KDE-based (4,000 draws). A bridge-sampling benchmark
# showed KDE errors of several log10 units when |log10 BF| >= 1, so the primary
# decomposition uses only replications with |log10 BF10| < 1 under BOTH priors.
# The full-sample version (excluding only 1e30 caps) is kept as a sensitivity
# column labelled "kde_tail_unreliable".

decomp <- decomp |>
  mutate(in_range = ok & abs(log10(BF10_logistic)) < 1 & abs(log10(BF10_cauchy)) < 1)

summarise_decomp <- function(df) {
  df |>
    summarise(
      n_used = n(),
      d_prior = first(d_prior),
      d_post_median = median(d_post),
      d_post_q25 = unname(quantile(d_post, 0.25)),
      d_post_q75 = unname(quantile(d_post, 0.75)),
      d_post_iqr = IQR(d_post),
      d_post_sd = sd(d_post),
      d_bf_median = median(d_bf),
      d_bf_iqr = IQR(d_bf),
      log10_bf_L_iqr = IQR(log10(BF10_logistic)),
      spearman_bf = if (n() >= 3) cor(BF10_logistic, BF10_cauchy, method = "spearman") else NA_real_,
      discordant_pairs = if (n() >= 3) {
        (1 - cor(log10(BF10_logistic), log10(BF10_cauchy), method = "kendall")) / 2
      } else NA_real_,
      .by = all_of(cell_vars)
    )
}

decomp_primary <- decomp |>
  filter(in_range) |>
  summarise_decomp()

decomp_full <- decomp |>
  filter(ok) |>
  summarise_decomp() |>
  select(all_of(cell_vars), n_used, d_post_median, d_post_iqr, spearman_bf) |>
  rename_with(\(x) paste0(x, "_kde_tail_unreliable"), -all_of(cell_vars))

decomp_summary <- decomp |>
  summarise(n_reps = n(),
            n_excluded_cap = sum(!ok),
            n_excluded_abs_log10_ge1 = sum(ok & !in_range),
            n_kept = sum(in_range),
            .by = all_of(cell_vars)) |>
  left_join(decomp_primary, by = join_by(sd_b1, true_b1, n_subjects, n_trials)) |>
  left_join(decomp_full, by = join_by(sd_b1, true_b1, n_subjects, n_trials)) |>
  arrange(sd_b1, true_b1, n_subjects, n_trials)

message("Decomposition restricted to |log10 BF10| < 1 under both priors: kept ",
        sum(decomp_summary$n_kept), " of ", sum(decomp_summary$n_reps),
        " replications (excluded: ", sum(decomp_summary$n_excluded_cap), " at cap, ",
        sum(decomp_summary$n_excluded_abs_log10_ge1), " with |log10 BF10| >= 1)")

saveRDS(list(summary = decomp_summary, r1_checks = r1_checks,
             replication_level = select(decomp, all_of(cell_vars), rep, ok, in_range,
                                        log10_post_L, log10_post_C, d_post, d_prior, d_bf)),
        here("output", "t1_sd_decomposition.rds"))

# --- Console summaries ----------------------------------------------------------

options(width = 200, pillar.sigfig = 3)
message("\n== Cap counts (nonzero only) ==")
print(filter(cap_counts, n_at_cap_1e30_underflow_or_fit_failure > 0), n = Inf)
message("\n== Paired family differences ==")
print(select(paired, all_of(cell_vars), P_BF10_gt3_logistic, P_BF10_gt3_cauchy,
             diff_BF10_gt3, diff_BF10_gt3_lo, diff_BF10_gt3_hi, disc_BF10_L_only, disc_BF10_C_only,
             mcnemar_p_BF10), n = Inf)
print(select(paired, all_of(cell_vars), P_BF01_gt3_logistic, P_BF01_gt3_cauchy,
             diff_BF01_gt3, diff_BF01_gt3_lo, diff_BF01_gt3_hi, disc_BF01_L_only, disc_BF01_C_only,
             mcnemar_p_BF01), n = Inf)
print(select(paired, all_of(cell_vars), n_excluded_cap, dlog10_median, dlog10_q25, dlog10_q75), n = Inf)
message("\n== Unpaired Logistic(0,.25) vs Cauchy(0,.707) ==")
print(unpaired, n = Inf)
message("\n== Decomposition (primary: |log10 BF10| < 1 under both priors) ==")
print(select(decomp_summary, all_of(cell_vars), n_kept, n_excluded_cap, n_excluded_abs_log10_ge1,
             d_post_median, d_post_q25, d_post_q75, d_post_iqr, d_bf_median, spearman_bf,
             discordant_pairs), n = Inf)
message("\n== Sensitivity: full sample excl. caps (KDE, unreliable in tail) ==")
print(select(decomp_summary, all_of(cell_vars), ends_with("_kde_tail_unreliable")), n = Inf)
message("\n== Pooled over cells, by scale (primary restriction) ==")
print(decomp |> filter(in_range) |>
        summarise(n = n(), d_post_median = median(d_post), d_post_iqr = IQR(d_post),
                  d_post_min = min(d_post), d_post_max = max(d_post),
                  d_bf_median = median(d_bf), .by = sd_b1))
message("DONE: script 11")
