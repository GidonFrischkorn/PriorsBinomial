# ==============================================================================
# 12_ConditionalTransforms.R
# ------------------------------------------------------------------------------
# Reviewer 3: prior on delta_p as a TRANSFORM at fixed baseline eta0 = b0 + u
# (Jacobian, no marginalisation over b0 or u).
#   1. Verify analytical density / CDF / quantile functions against Monte Carlo
#   2. Conditional summaries on a grid (nominal and IQR-matched scales)
#   3. Participants at u = -sigma_u, 0, +sigma_u (+ draft figure)
#   4. Criterion sensitivity (intercept floor/ceiling mass; effect percentiles)
#   5. Probit check: Phi(b0 + b1) uniform only if s0^2 + s1^2 = 1
# Note: `scale` is the family's scale parameter. For Logistic it is NOT an SD
# (SD = scale * pi / sqrt(3)).
# ==============================================================================

library(here)
library(ggplot2)

source(here("R", "link_functions.R"))
source(here("R", "prior_sampling.R"))

mc_seed <- 20261005
out_dir <- here("output")
fig_dir <- here("figures")
dir.create(fig_dir, showWarnings = FALSE)

families <- c("normal", "logistic", "cauchy")
links    <- c("logit", "probit")

# ------------------------------------------------------------------------------
# 1. Verification against Monte Carlo (1e6 draws)
# ------------------------------------------------------------------------------
n_mc <- 1e6
verify_combos <- data.frame(
  link   = c("logit", "logit", "logit", "logit", "probit", "probit", "probit", "probit"),
  family = c("logistic", "cauchy", "normal", "cauchy", "normal", "cauchy", "logistic", "normal"),
  scale  = c(0.25, 0.707, 0.50, 0.10, 0.25, 0.10, 0.30, 0.707),
  p0     = c(0.70, 0.90, 0.60, 0.50, 0.50, 0.80, 0.70, 0.90),
  stringsAsFactors = FALSE
)

verify_one <- function(link, family, scale, p0, sat = 1 - 1e-8) {
  eta0 <- apply_link(p0, link)
  b1   <- sample_prior(n_mc, family, location = 0, scale = scale)
  dl   <- delta_p_from_b1(b1, eta0, link)
  s_dl  <- sort(dl)
  s_abs <- sort(abs(dl))
  n <- length(dl)
  # right-continuous ECDF and its left limit (tie-aware: heavy Cauchy tails
  # produce draws with delta_p exactly +-1 in double precision)
  F_r <- function(x, s) findInterval(x, s) / n
  F_l <- function(x, s) findInterval(x, s, left.open = TRUE) / n
  # distance of a target probability from the ECDF jump interval at x
  tie_dist <- function(x, target, s) pmax(0, F_l(x, s) - target, target - F_r(x, s))

  probs_grid <- seq(0.0001, 0.9999, length.out = 5001)
  pq <- seq(0.001, 0.999, by = 0.001)

  # CDF of delta_p on a dense grid of empirical quantiles
  d_grid <- unique(quantile(dl, probs_grid, names = FALSE))
  cdf_err <- abs(F_r(d_grid, s_dl) - p_delta_p(d_grid, eta0, link, family, scale))
  in_d <- abs(d_grid) < sat

  # density: CDF built by integrating d_delta_p from 0 (F(0) = 0.5 by symmetry)
  d_int <- sort(unique(c(0, quantile(dl, seq(0.001, 0.999, length.out = 801),
                                     names = FALSE))))
  # drop points where delta_p is numerically indistinguishable from +-1;
  # b1(d) is not recoverable there in double precision
  d_int <- d_int[abs(d_int) < sat]
  i0 <- which(d_int == 0)
  seg <- vapply(seq_len(length(d_int) - 1), function(i) {
    integrate(d_delta_p, d_int[i], d_int[i + 1], eta0 = eta0, link = link,
              family = family, scale = scale, rel.tol = 1e-10)$value
  }, numeric(1))
  cdf_int <- numeric(length(d_int))
  cdf_int[i0] <- 0.5
  if (i0 < length(d_int)) cdf_int[(i0 + 1):length(d_int)] <- 0.5 + cumsum(seg[i0:length(seg)])
  if (i0 > 1) cdf_int[(i0 - 1):1] <- 0.5 - cumsum(rev(seg[1:(i0 - 1)]))
  ks_density <- max(abs(F_r(d_int, s_dl) - cdf_int))

  # quantile function (tie-aware)
  qd <- q_delta_p(pq, eta0, link, family, scale)
  q_err <- tie_dist(qd, pq, s_dl)
  in_q <- abs(qd) < sat

  # |delta_p| CDF and quantile function
  a_grid <- unique(quantile(abs(dl), probs_grid, names = FALSE))
  abs_cdf_err <- abs(F_r(a_grid, s_abs) - p_abs_delta_p(a_grid, eta0, link, family, scale))
  in_a <- a_grid < sat
  qa <- q_abs_delta_p(pq, eta0, link, family, scale)
  abs_q_err <- tie_dist(qa, pq, s_abs)
  in_qa <- qa < sat

  data.frame(
    link = link, family = family, scale = scale, p0 = p0,
    ks_cdf = max(cdf_err[in_d]), ks_density = ks_density,
    ks_quantile = max(q_err[in_q]),
    ks_abs_cdf = max(abs_cdf_err[in_a]), ks_abs_quantile = max(abs_q_err[in_qa]),
    sat_mass_mc = mean(abs(dl) >= sat),
    sat_mass_analytic = 1 - p_abs_delta_p(sat, eta0, link, family, scale),
    ks_unrestricted_max = max(cdf_err, q_err, abs_cdf_err, abs_q_err)
  )
}

verification <- withr::with_seed(mc_seed, do.call(rbind, Map(verify_one, verify_combos$link, verify_combos$family,
                                   verify_combos$scale, verify_combos$p0)))
rownames(verification) <- NULL
verification$max_ks <- apply(verification[, c("ks_cdf", "ks_density", "ks_quantile",
                                              "ks_abs_cdf", "ks_abs_quantile")], 1, max)
verification$sat_mass_diff <- abs(verification$sat_mass_mc - verification$sat_mass_analytic)
# verified: KS < 0.005 on |delta_p| < 1 - 1e-8 AND saturated mass agrees within 0.005
verification$verified <- verification$max_ks < 0.005 & verification$sat_mass_diff < 0.005
saveRDS(verification, file.path(out_dir, "t2_verification.rds"))
cat("\n== 1. Verification (KS distances vs 1e6 MC draws) ==\n")
print(verification, digits = 3)

# ------------------------------------------------------------------------------
# 2. Conditional summaries on a grid
# ------------------------------------------------------------------------------
summarise_cond <- function(link, family, scale, eta0) {
  data.frame(
    p_gt_0.1 = 1 - p_abs_delta_p(0.1, eta0, link, family, scale),
    p_gt_0.2 = 1 - p_abs_delta_p(0.2, eta0, link, family, scale),
    p_gt_0.3 = 1 - p_abs_delta_p(0.3, eta0, link, family, scale),
    q50 = q_abs_delta_p(0.50, eta0, link, family, scale),
    q90 = q_abs_delta_p(0.90, eta0, link, family, scale),
    q95 = q_abs_delta_p(0.95, eta0, link, family, scale)
  )
}

scale_grid <- c(0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.707)
p0_grid    <- c(0.5, 0.6, 0.7, 0.8, 0.9)

grid_nominal <- expand.grid(link = links, family = families, scale = scale_grid,
                            p0 = p0_grid, stringsAsFactors = FALSE)
grid_nominal$scale_type <- "nominal"
grid_nominal$logistic_ref_scale <- NA_real_

grid_matched <- expand.grid(link = links, family = families,
                            logistic_ref_scale = scale_grid, p0 = p0_grid,
                            stringsAsFactors = FALSE)
grid_matched$scale <- mapply(matched_iqr_scale, grid_matched$logistic_ref_scale,
                             grid_matched$family)
grid_matched$scale_type <- "iqr_matched_to_logistic"

grid_all <- rbind(grid_nominal[, c("scale_type", "link", "family", "scale", "logistic_ref_scale", "p0")],
                  grid_matched[, c("scale_type", "link", "family", "scale", "logistic_ref_scale", "p0")])
grid_summ <- do.call(rbind, lapply(seq_len(nrow(grid_all)), function(i) {
  r <- grid_all[i, ]
  summarise_cond(r$link, r$family, r$scale, apply_link(r$p0, r$link))
}))
conditional_grid <- cbind(grid_all, grid_summ)
saveRDS(conditional_grid, file.path(out_dir, "t2_conditional_grid.rds"))

equiv_scales <- data.frame(
  logistic_scale = scale_grid,
  normal_scale   = matched_iqr_scale(scale_grid, "normal"),
  cauchy_scale   = matched_iqr_scale(scale_grid, "cauchy"),
  b1_iqr         = 2 * scale_grid * log(3)
)
saveRDS(equiv_scales, file.path(out_dir, "t2_matched_iqr_scales.rds"))
cat("\n== 2a. IQR-matched equivalent scales ==\n")
print(equiv_scales, digits = 4)
cat("\n== 2b. Conditional grid (logit, nominal, scale 0.25/0.30) ==\n")
print(subset(conditional_grid, scale_type == "nominal" & link == "logit" &
               scale %in% c(0.25, 0.30)), digits = 3)
cat("\n== 2c. Conditional grid (logit, IQR-matched to Logistic 0.25) ==\n")
print(subset(conditional_grid, scale_type != "nominal" & link == "logit" &
               logistic_ref_scale == 0.25), digits = 3)

# ------------------------------------------------------------------------------
# 3. Participants at u = -sigma_u, 0, +sigma_u
# ------------------------------------------------------------------------------
part_grid <- expand.grid(link = links, family = families, scale = c(0.25, 0.30),
                         p0 = c(0.5, 0.7, 0.9), sigma_u = c(0.25, 0.5, 1.0),
                         u_sd = c(-1, 0, 1), stringsAsFactors = FALSE)
part_grid$u    <- part_grid$u_sd * part_grid$sigma_u
part_grid$eta0 <- mapply(apply_link, part_grid$p0, part_grid$link) + part_grid$u
part_grid$p_baseline <- mapply(apply_inverse_link, part_grid$eta0, part_grid$link)
part_summ <- do.call(rbind, lapply(seq_len(nrow(part_grid)), function(i) {
  r <- part_grid[i, ]
  summarise_cond(r$link, r$family, r$scale, r$eta0)
}))
participants <- cbind(part_grid, part_summ)
saveRDS(participants, file.path(out_dir, "t2_participants.rds"))
cat("\n== 3. Participants (logit, Logistic 0.25) ==\n")
print(subset(participants, link == "logit" & family == "logistic" & scale == 0.25)[
  , c("p0", "sigma_u", "u", "p_baseline", "p_gt_0.1", "p_gt_0.2", "q50", "q90", "q95")],
  digits = 3, row.names = FALSE)

# Draft figure: Logistic(0, 0.25), sigma_u = 0.5, p0 = 0.7, logit
fig_sigma <- 0.5
fig_p0    <- 0.7
d_seq <- seq(-0.6, 0.6, length.out = 1201)
fig_df <- do.call(rbind, lapply(c(-1, 0, 1), function(k) {
  eta0 <- logit(fig_p0) + k * fig_sigma
  data.frame(
    delta_p = d_seq,
    density = d_delta_p(d_seq, eta0, "logit", "logistic", 0.25),
    participant = sprintf("u = %+.0f SD (baseline p = %.3f)", k, ilogit(eta0))
  )
}))
fig_df$participant <- factor(fig_df$participant, levels = unique(fig_df$participant))
p_fig <- ggplot(fig_df, aes(delta_p, density, colour = participant, linetype = participant)) +
  geom_line(linewidth = 0.8) +
  scale_colour_viridis_d(end = 0.85) +
  labs(x = expression(Delta * p ~ "at fixed baseline"), y = "Density",
       colour = NULL, linetype = NULL,
       title = "Prior on delta_p: Logistic(0, scale 0.25) effect, logit link",
       subtitle = "Baseline p0 = 0.7 for u = 0; sigma_u = 0.5") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "vertical")
ggsave(file.path(fig_dir, "t2_delta_p_by_participant.png"), p_fig,
       width = 6.5, height = 4.5, dpi = 300, bg = "white")

# ------------------------------------------------------------------------------
# 4. Criterion sensitivity
# ------------------------------------------------------------------------------
# 4a. Intercept prior: exact floor/ceiling mass F(g(lo)) + 1 - F(g(hi))
b0_scales <- round(seq(0.50, 1.25, by = 0.05), 2)
bounds <- data.frame(lower = c(0.01, 0.025, 0.05, 0.10),
                     upper = c(0.99, 0.975, 0.95, 0.90))
ic_grid <- expand.grid(link = links, family = c("logistic", "normal"),
                       scale = b0_scales, bound_id = seq_len(nrow(bounds)),
                       stringsAsFactors = FALSE)
ic_grid$lower <- bounds$lower[ic_grid$bound_id]
ic_grid$upper <- bounds$upper[ic_grid$bound_id]
ic_grid$mass <- mapply(floor_ceiling_mass, ic_grid$scale, ic_grid$family,
                       ic_grid$link, ic_grid$lower, ic_grid$upper)
ic_grid$bound_id <- NULL
saveRDS(ic_grid, file.path(out_dir, "t2_intercept_floor_ceiling.rds"))

# exact largest scale meeting mass <= threshold (symmetric bounds):
# 2 F0(g(lo) / s) <= t  <=>  s <= g(lo) / F0^{-1}(t / 2)
thr_grid <- expand.grid(link = links, family = c("logistic", "normal"),
                        lower = bounds$lower, threshold = c(0.05, 0.10, 0.15),
                        stringsAsFactors = FALSE)
thr_grid$upper <- 1 - thr_grid$lower
thr_grid$max_scale_exact <- mapply(function(link, family, lower, t) {
  apply_link(lower, link) / effect_prior_fun(t / 2, family, 1, "q")
}, thr_grid$link, thr_grid$family, thr_grid$lower, thr_grid$threshold)
thr_grid$max_scale_on_grid <- vapply(seq_len(nrow(thr_grid)), function(i) {
  r <- thr_grid[i, ]
  ok <- b0_scales[floor_ceiling_mass(b0_scales, r$family, r$link, r$lower, r$upper) <= r$threshold]
  if (length(ok) == 0) NA_real_ else max(ok)
}, numeric(1))
thr_grid$mass_at_0.75 <- mapply(floor_ceiling_mass, 0.75, thr_grid$family, thr_grid$link,
                                thr_grid$lower, thr_grid$upper)
thr_grid$meets_at_0.75 <- thr_grid$mass_at_0.75 <= thr_grid$threshold
saveRDS(thr_grid, file.path(out_dir, "t2_intercept_thresholds.rds"))
cat("\n== 4a. Intercept: largest scale meeting floor/ceiling mass threshold ==\n")
print(thr_grid[order(thr_grid$family, thr_grid$link, thr_grid$lower, thr_grid$threshold), ],
      digits = 3, row.names = FALSE)

# 4b. Effect prior: scale that puts the q-th percentile of |delta_p| at target
# Q_|dp|(q) = delta_p(s * Q0((1 + q) / 2))  =>  s = b1(target) / Q0((1 + q) / 2)
eff_grid <- expand.grid(link = links, family = families, p0 = c(0.5, 0.7),
                        percentile = c(0.80, 0.90, 0.95), target = c(0.20, 0.25, 0.30),
                        stringsAsFactors = FALSE)
eff_grid$scale_required <- mapply(function(link, family, p0, q, target) {
  eta0 <- apply_link(p0, link)
  b1_from_delta_p(target, eta0, link) / effect_prior_fun((1 + q) / 2, family, 1, "q")
}, eff_grid$link, eff_grid$family, eff_grid$p0, eff_grid$percentile, eff_grid$target)
# round-trip check
eff_grid$check <- mapply(function(link, family, p0, q, s) {
  q_abs_delta_p(q, apply_link(p0, link), link, family, s)
}, eff_grid$link, eff_grid$family, eff_grid$p0, eff_grid$percentile, eff_grid$scale_required)
stopifnot(max(abs(eff_grid$check - eff_grid$target)) < 1e-8)
eff_grid$check <- NULL
saveRDS(eff_grid, file.path(out_dir, "t2_effect_percentile_scales.rds"))
cat("\n== 4b. Effect: scale giving percentile of |delta_p| = target ==\n")
print(eff_grid[order(eff_grid$link, eff_grid$family, eff_grid$p0, eff_grid$percentile), ],
      digits = 3, row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. Probit check (Reviewer 1)
# ------------------------------------------------------------------------------
probit_check <- function(s0, s1, n = 1e6) {
  b0 <- rnorm(n, 0, s0)
  b1 <- rnorm(n, 0, s1)
  p  <- pnorm(b0 + b1)
  s_tot <- sqrt(s0^2 + s1^2)
  u <- seq(0.0005, 0.9995, by = 0.0005)
  data.frame(
    s0 = s0, s1 = s1,
    sd_total_analytic = s_tot,
    sd_total_mc = sd(b0 + b1),
    ks_uniform_mc = suppressWarnings(ks.test(p, "punif")$statistic[[1]]),
    ks_uniform_exact = max(abs(pnorm(qnorm(u) / s_tot) - u)),
    mass_out_0.05_0.95_exact = 2 * pnorm(qnorm(0.05) / s_tot),
    mass_out_0.05_0.95_mc = mean(p < 0.05 | p > 0.95)
  )
}
probit_res <- withr::with_seed(mc_seed + 1, rbind(
  probit_check(0.75, 0.25),
  probit_check(sqrt(1 - 0.25^2), 0.25),
  probit_check(0.75, 0),
  probit_check(1, 0)
))
saveRDS(probit_res, file.path(out_dir, "t2_probit_check.rds"))
cat("\n== 5. Probit check ==\n")
print(probit_res, digits = 4)

cat("\nT2_DONE\n")
