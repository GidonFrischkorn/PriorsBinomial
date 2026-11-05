# rs_postprocess.R
# Helpers for post-processing the random-slopes BF calibration (Scripts 15/15b).
# All BF-based rates treat NA (failed fit) separately from Inf (extreme
# evidence): NA is counted in n_fail and never mapped to a number.
#
# The second prior parameter is a SCALE (Logistic SD = scale * pi / sqrt(3)).

#' Reshape the raw per-replication results of Script 13/14 to long format
#'
#' The raw object has one row per replication and one column per
#' `<stat>_<family>_<scale>` (e.g. `bf10_bridge_logistic_0.25`). Arms have
#' different prior grids, so SimResults() pads other arms' columns with NA;
#' only the priors of each arm's own grid are kept (`make_prior_grids_rs()`).
#'
#' @param raw Data frame as written by Script 14 (`res_bf_calibration_rs_raw.rds`)
#'   or `SimResults()` of a Script 13 run with the bridge engine.
#' @param prior_grids Named list of per-arm prior grids (`make_prior_grids_rs()`).
#' @param stats Character vector of per-prior statistic prefixes to keep.
#' @return Tibble with one row per replication x prior: `arm`, `true_b1`,
#'   `n_subjects`, `n_trials`, `true_sd_slope`, `design_row`, `rep` (1-based
#'   index within the design row, in the order stored by SimDesign),
#'   `family`, `scale`, `prior_label`, one column per element of `stats`,
#'   `log10_bf10` (bridge; NA for failed fits, Inf/-Inf kept) and the
#'   replication-level columns `n_fail`, `true_p0`, `data_sum_y`. Priors of
#'   the grid without a `bf10_bridge_*` column in `raw` (not fitted in that
#'   run) are dropped.
rs_raw_to_long <- function(raw, prior_grids = make_prior_grids_rs(),
                           stats = c("bf10_bridge", "bf10_logspline",
                                     "pct_error_h1", "post_mean", "post_sd",
                                     "n_div", "rhat_max", "ess_bulk_min",
                                     "ess_tail_min", "wall_s")) {
  raw <- as.data.frame(raw)
  design_cols <- c("arm", "true_b1", "n_subjects", "n_trials", "true_sd_slope",
                   "design_row")
  stopifnot(all(design_cols %in% names(raw)))
  raw <- raw[order(raw$design_row), ]
  raw$rep <- stats::ave(seq_len(nrow(raw)), raw$design_row, FUN = seq_along)

  grid_long <- purrr::imap(prior_grids, function(g, a) {
    data.frame(arm = a, family = as.character(g$dist_b1), scale = g$scale_b1,
               scale_chr = format_scale_rs(g$scale_b1), stringsAsFactors = FALSE)
  }) |> purrr::list_rbind()

  rep_cols <- intersect(c("n_fail", "true_p0", "data_sum_y", "rep_wall_s",
                          "wall_s_h0"), names(raw))
  per_prior <- purrr::pmap(grid_long, function(arm, family, scale, scale_chr) {
    rows <- raw[raw$arm == arm, , drop = FALSE]
    # A prior that was not fitted in this run (e.g. the 3-prior smoke test)
    # has no column, or an all-NA column padded in by SimResults(); it is
    # dropped, not counted as failed. A failed fit is NA in every stat column
    # of one replication only, and is counted in n_fail (checked below).
    bf_col <- paste0("bf10_bridge_", family, "_", scale_chr)
    if (nrow(rows) == 0 || !bf_col %in% names(rows)) return(NULL)
    all_na <- all(is.na(rows[[bf_col]])) &&
      all(is.na(rows[[paste0("wall_s_", family, "_", scale_chr)]]))
    if (all_na) return(NULL)
    out <- rows[, c(design_cols, "rep", rep_cols)]
    out$family <- family
    out$scale  <- scale
    for (s in stats) {
      col <- paste0(s, "_", family, "_", scale_chr)
      out[[s]] <- if (col %in% names(rows)) rows[[col]] else NA_real_
    }
    out
  }) |> purrr::list_rbind() |> tibble::as_tibble()

  fitted <- unique(per_prior[, c("arm", "family", "scale")])
  message("Priors with results per arm: ",
          paste(names(table(fitted$arm)), table(fitted$arm), sep = " = ", collapse = ", "))
  # Integrity: the NA BFs among kept priors are exactly the failures Analyse counted
  n_na_bf <- sum(is.na(per_prior$bf10_bridge))
  n_fail_total <- if ("n_fail" %in% names(raw)) sum(raw$n_fail) else NA
  if (!is.na(n_fail_total) && n_na_bf != n_fail_total) {
    stop("NA bf10 values among fitted priors (", n_na_bf, ") differ from total n_fail (",
         n_fail_total, "): a prior was dropped or padded wrongly")
  }
  per_prior$prior_label <- prior_label_rs(per_prior$family, per_prior$scale)
  per_prior$log10_bf10  <- log10(per_prior$bf10_bridge)
  per_prior[order(per_prior$design_row, per_prior$rep, per_prior$family,
                  per_prior$scale), ]
}

#' Prior label "Family(0, scale)"
#'
#' @param family Character vector of family names (lower case).
#' @param scale Numeric vector of scales.
#' @return Character vector such as "Logistic(0, 0.25)".
prior_label_rs <- function(family, scale) {
  fam <- c(normal = "Normal", logistic = "Logistic", cauchy = "Cauchy")[family]
  sprintf("%s(0, %s)", fam, format_scale_rs(scale))
}

#' Scale formatted as in the column names of Script 13 (element-wise)
#'
#' `fit_prior_grid_rs()` formats each scale on its own with
#' `format(scale, nsmall = 2)`, giving "0.15" and "0.707"; formatting a vector
#' at once would pad to a common width ("0.150"), so this is done element-wise.
#'
#' @param scale Numeric vector.
#' @return Character vector.
format_scale_rs <- function(scale) {
  vapply(scale, function(s) format(s, nsmall = 2), character(1))
}

#' Wilson score interval for a binomial proportion
#'
#' @param k Integer successes.
#' @param n Integer trials.
#' @param level Confidence level.
#' @return Numeric vector `c(lo, hi)`; `NA` when `n == 0`.
wilson_ci <- function(k, n, level = 0.95) {
  if (is.na(n) || n == 0) return(c(NA_real_, NA_real_))
  z  <- stats::qnorm(1 - (1 - level) / 2)
  p  <- k / n
  dn <- 1 + z^2 / n
  ce <- (p + z^2 / (2 * n)) / dn
  hw <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / dn
  c(max(0, ce - hw), min(1, ce + hw))
}

#' Area under the ROC curve by the Wilcoxon (rank-sum) statistic
#'
#' `AUC = P(score_signal > score_noise) + 0.5 * P(tie)`. Infinite scores are
#' allowed (they rank above every finite score); NA scores must be removed
#' by the caller.
#'
#' @param signal Numeric scores of the signal class (true effect).
#' @param noise Numeric scores of the noise class (null).
#' @return Numeric AUC in [0, 1]; `NA` if either vector is empty.
auc_wilcoxon <- function(signal, noise) {
  n_s <- length(signal); n_n <- length(noise)
  if (n_s == 0 || n_n == 0) return(NA_real_)
  stopifnot(!anyNA(signal), !anyNA(noise))
  r <- rank(c(signal, noise))               # average ranks for ties
  (sum(r[seq_len(n_s)]) - n_s * (n_s + 1) / 2) / (n_s * n_n)
}

#' Empirical ROC curve interpolated on a false-positive-rate grid
#'
#' For every threshold t (all observed scores), FPR = P(noise > t) and
#' TPR = P(signal > t); the step curve is linearly interpolated at `grid`.
#' Same construction as the earlier JZS comparison, vectorised with `ecdf()`.
#'
#' @param noise,signal Numeric scores without NA (Inf allowed).
#' @param grid Numeric vector of FPR values in [0, 1].
#' @return Numeric TPR at each grid point.
roc_at_grid <- function(noise, signal, grid) {
  thr <- sort(unique(c(noise, signal)), decreasing = TRUE)
  thr <- thr[is.finite(thr)]
  f_n <- stats::ecdf(noise); f_s <- stats::ecdf(signal)
  # P(x > t) = 1 - F(t); ecdf is right-continuous so this is exact at observed
  # values. Inf scores exceed every finite threshold.
  fpr <- c(1, 1 - f_n(thr), 0)
  tpr <- c(1, 1 - f_s(thr), 0)
  ord <- order(fpr, tpr)
  fpr <- fpr[ord]; tpr <- tpr[ord]
  keep <- !duplicated(fpr, fromLast = TRUE)   # at equal FPR keep the max TPR
  stats::approx(fpr[keep], tpr[keep], xout = grid, method = "linear",
                yleft = 0, yright = 1, ties = "ordered")$y
}

#' Scale of another family with the same inter-quartile range as Logistic(0, s)
#'
#' IQR(Logistic(0, s)) = 2 s log 3; IQR(Normal(0, sd)) = 2 sd qnorm(0.75);
#' IQR(Cauchy(0, g)) = 2 g. Checked numerically against `qlogis()`,
#' `qnorm()` and `qcauchy()` by `check_equal_iqr_formulas()`.
#'
#' @param s Numeric Logistic scale(s).
#' @param family Character, "normal" or "cauchy".
#' @return Numeric scale(s) of `family` with the same IQR.
equal_iqr_scale <- function(s, family = c("normal", "cauchy")) {
  family <- match.arg(family)
  switch(family,
         normal = s * 2 * log(3) / (2 * stats::qnorm(0.75)),
         cauchy = s * log(3))
}

#' Numerical check of the equal-IQR formulas and prior ordinates
#'
#' @param s Numeric vector of Logistic scales to check.
#' @return Invisibly `TRUE`; stops if any identity fails.
check_equal_iqr_formulas <- function(s = c(0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.707)) {
  iqr <- function(q, ...) q(0.75, ...) - q(0.25, ...)
  iqr_l <- iqr(stats::qlogis, location = 0, scale = s)
  iqr_n <- iqr(stats::qnorm, mean = 0, sd = equal_iqr_scale(s, "normal"))
  iqr_c <- iqr(stats::qcauchy, location = 0, scale = equal_iqr_scale(s, "cauchy"))
  stopifnot(isTRUE(all.equal(iqr_l, iqr_n)), isTRUE(all.equal(iqr_l, iqr_c)))
  stopifnot(isTRUE(all.equal(stats::dlogis(0, 0, s), 1 / (4 * s))),
            isTRUE(all.equal(stats::dcauchy(0, 0, s), 1 / (pi * s))),
            isTRUE(all.equal(stats::dnorm(0, 0, s), 1 / (s * sqrt(2 * pi)))))
  invisible(TRUE)
}

#' Self-check of the ROC/AUC helpers on known inputs
#'
#' @return Invisibly `TRUE`; stops on failure.
check_roc_helpers <- function() {
  stopifnot(auc_wilcoxon(c(4, 5, 6), c(1, 2, 3)) == 1,
            auc_wilcoxon(c(1, 2, 3), c(4, 5, 6)) == 0,
            auc_wilcoxon(c(1, 2), c(1, 2)) == 0.5,
            auc_wilcoxon(c(Inf, 5), c(1, Inf)) == 0.625)
  g <- c(0, 0.5, 1)
  stopifnot(isTRUE(all.equal(roc_at_grid(c(1, 2, 3), c(4, 5, 6), g), c(1, 1, 1))),
            isTRUE(all.equal(roc_at_grid(c(1, 2, 3), c(1, 2, 3), g), g)))
  w <- wilson_ci(5, 10)
  stopifnot(w[1] > 0.2, w[1] < 0.5, w[2] > 0.5, w[2] < 0.8,
            all(is.na(wilson_ci(0, 0))))
  invisible(TRUE)
}
