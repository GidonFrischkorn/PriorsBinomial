# ==============================================================================
# Dual ROC Comparison: H1 detection and H0 detection
# ------------------------------------------------------------------------------
# H1 detection: score = BF10, signal = H1 data, noise = H0 data.
#   TPR = P(BF10 > t | H1),  FPR = P(BF10 > t | H0)
# H0 detection: score = BF01 = 1/BF10, signal = H0 data, noise = H1 data.
#   Derived analytically: (fpr_H0, tpr_H0) = (1 - tpr_H1, 1 - fpr_H1)
# Both ROC curves share the same AUC — discriminability is invariant to
# detection direction as well as to prior scale and family.
# ==============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(here)

df        <- readRDS(here("output", "bf_calibration_jzs_raw.rds"))
set.seed(42)
n_boot    <- 500
fpr_grid  <- seq(0, 1, length.out = 100)
jzs_scale <- round(sqrt(2) / 2, 3)

make_prior_label <- function(dist, sc) {
  sprintf("%s\n(0, %s)", ifelse(dist == "logistic", "Logistic", "Cauchy"), sc)
}

prior_levels <- c(
  "Logistic\n(0, 0.25)", "Cauchy\n(0, 0.25)",
  "Logistic\n(0, 0.707)", "Cauchy\n(0, 0.707)"
)

# Viridis (option D, end = 0.85) to match the in-document figures:
# scale 0.25 -> dark half (purple/blue), scale 0.707 -> light half (teal/green).
prior_colors <- c(
  "Logistic\n(0, 0.25)"  = "#440154",
  "Cauchy\n(0, 0.25)"    = "#375B8D",
  "Logistic\n(0, 0.707)" = "#1FA188",
  "Cauchy\n(0, 0.707)"   = "#9AD93C"
)

roc_at_grid <- function(null_vals, alt_vals, grid) {
  all_t  <- sort(unique(c(null_vals, alt_vals)), decreasing = TRUE)
  fpr_pt <- c(1, sapply(all_t, function(t) mean(null_vals > t)), 0)
  tpr_pt <- c(1, sapply(all_t, function(t) mean(alt_vals  > t)), 0)
  ord    <- order(fpr_pt)
  fpr_s  <- fpr_pt[ord]; tpr_s <- tpr_pt[ord]
  keep   <- !duplicated(fpr_s)
  approx(fpr_s[keep], tpr_s[keep], xout = grid,
         method = "linear", yleft = 0, yright = 1)$y
}

# --- Bootstrap H1 ROC curves --------------------------------------------------

results <- list()
idx     <- 1

for (sc in c(0.25, jzs_scale)) {
  for (dist in c("logistic", "cauchy")) {
    for (ns in c(20, 30, 50)) {
      for (nt in c(20, 50)) {
        nv <- df$BF10[df$prior == dist & df$sd_b1 == sc &
                      df$n_subjects == ns & df$n_trials == nt & df$true_b1 == 0]
        for (tb1 in c(0.1, 0.2)) {
          av <- df$BF10[df$prior == dist & df$sd_b1 == sc &
                        df$n_subjects == ns & df$n_trials == nt & df$true_b1 == tb1]
          if (length(nv) == 0 || length(av) == 0) next

          tpr_orig <- roc_at_grid(nv, av, fpr_grid)
          boot_mat <- matrix(NA_real_, n_boot, length(fpr_grid))
          for (b in seq_len(n_boot)) {
            boot_mat[b, ] <- roc_at_grid(sample(nv, replace = TRUE),
                                         sample(av, replace = TRUE),
                                         fpr_grid)
          }

          results[[idx]] <- data.frame(
            prior = dist, sd_b1 = sc, n_subjects = ns, n_trials = nt, true_b1 = tb1,
            fpr   = fpr_grid,
            tpr   = tpr_orig,
            ci_lo = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE),
            ci_hi = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE)
          )
          idx <- idx + 1
        }
      }
    }
  }
}

h1_roc <- bind_rows(results) |>
  mutate(
    prior_label = factor(make_prior_label(prior, sd_b1), levels = prior_levels),
    b1_label    = factor(paste0("b₁ = ", true_b1),
                         levels = c("b₁ = 0.1", "b₁ = 0.2"))
  )

# --- Bootstrap H0 ROC curves directly -----------------------------------------
# Score = BF01 = 1/BF10; signal = H0 data (b1=0); noise = H1 data (b1>0).
# CI bands are on the TPR axis at each fixed FPR grid point, same as H1 ROC.

h0_results <- list()
idx_h0     <- 1

for (sc in c(0.25, jzs_scale)) {
  for (dist in c("logistic", "cauchy")) {
    for (ns in c(20, 30, 50)) {
      for (nt in c(20, 50)) {
        h0_vals <- 1 / df$BF10[df$prior == dist & df$sd_b1 == sc &
                                df$n_subjects == ns & df$n_trials == nt &
                                df$true_b1 == 0]
        for (tb1 in c(0.1, 0.2)) {
          h1_vals <- 1 / df$BF10[df$prior == dist & df$sd_b1 == sc &
                                  df$n_subjects == ns & df$n_trials == nt &
                                  df$true_b1 == tb1]
          if (length(h0_vals) == 0 || length(h1_vals) == 0) next

          tpr_orig <- roc_at_grid(h1_vals, h0_vals, fpr_grid)
          boot_mat <- matrix(NA_real_, n_boot, length(fpr_grid))
          for (b in seq_len(n_boot)) {
            boot_mat[b, ] <- roc_at_grid(sample(h1_vals, replace = TRUE),
                                         sample(h0_vals, replace = TRUE),
                                         fpr_grid)
          }

          h0_results[[idx_h0]] <- data.frame(
            prior = dist, sd_b1 = sc, n_subjects = ns, n_trials = nt, true_b1 = tb1,
            fpr   = fpr_grid,
            tpr   = tpr_orig,
            ci_lo = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE),
            ci_hi = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE)
          )
          idx_h0 <- idx_h0 + 1
        }
      }
    }
  }
}

h0_roc <- bind_rows(h0_results) |>
  mutate(
    prior_label = factor(make_prior_label(prior, sd_b1), levels = prior_levels),
    b1_label    = factor(paste0("b₁ = ", true_b1),
                         levels = c("b₁ = 0.1", "b₁ = 0.2"))
  )

# --- BF criterion points ------------------------------------------------------

# H1: BF10 > 3
bf3_h1 <- df |>
  filter(true_b1 %in% c(0.1, 0.2)) |>
  group_by(prior, sd_b1, n_subjects, n_trials, true_b1) |>
  summarise(tpr_3 = mean(BF10 > 3), .groups = "drop") |>
  left_join(
    df |>
      filter(true_b1 == 0) |>
      group_by(prior, sd_b1, n_subjects, n_trials) |>
      summarise(fpr_3 = mean(BF10 > 3), .groups = "drop"),
    by = c("prior", "sd_b1", "n_subjects", "n_trials")
  ) |>
  mutate(
    prior_label = factor(make_prior_label(prior, sd_b1), levels = prior_levels),
    b1_label    = factor(paste0("b₁ = ", true_b1), levels = c("b₁ = 0.1", "b₁ = 0.2"))
  )

# H0: BF01 > 3, i.e., BF10 < 1/3
bf3_h0 <- df |>
  filter(true_b1 %in% c(0.1, 0.2)) |>
  group_by(prior, sd_b1, n_subjects, n_trials, true_b1) |>
  summarise(fpr_3 = mean(BF10 < 1/3), .groups = "drop") |>  # FPR_H0: H1 data flagged as H0
  left_join(
    df |>
      filter(true_b1 == 0) |>
      group_by(prior, sd_b1, n_subjects, n_trials) |>
      summarise(tpr_3 = mean(BF10 < 1/3), .groups = "drop"),  # TPR_H0: H0 data correctly flagged
    by = c("prior", "sd_b1", "n_subjects", "n_trials")
  ) |>
  mutate(
    prior_label = factor(make_prior_label(prior, sd_b1), levels = prior_levels),
    b1_label    = factor(paste0("b₁ = ", true_b1), levels = c("b₁ = 0.1", "b₁ = 0.2"))
  )

# --- Shared plot elements -----------------------------------------------------

roc_panel <- function(roc_df, bf3_df, title) {
  ggplot(roc_df,
         aes(x = fpr, y = tpr,
             color = prior_label, fill = prior_label,
             linetype = b1_label)) +
    geom_abline(slope = 1, intercept = 0, color = "grey70", linetype = "dotted") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
                alpha = 0.08, color = NA, show.legend = FALSE) +
    geom_line(linewidth = 0.7) +
    geom_segment(data = bf3_df,
                 aes(x = fpr_3, xend = fpr_3, y = 0, yend = tpr_3),
                 linetype = "dotted", linewidth = 0.4, color = "grey40",
                 inherit.aes = FALSE) +
    geom_segment(data = bf3_df,
                 aes(x = 0, xend = fpr_3, y = tpr_3, yend = tpr_3),
                 linetype = "dotted", linewidth = 0.4, color = "grey40",
                 inherit.aes = FALSE) +
    geom_point(data = bf3_df,
               aes(x = fpr_3, y = tpr_3, color = prior_label),
               size = 2.5, shape = 21, fill = "white", stroke = 1.2,
               show.legend = FALSE, inherit.aes = FALSE) +
    facet_grid(n_subjects ~ n_trials,
               labeller = labeller(
                 n_subjects = as_labeller(function(x) paste0("N = ", x)),
                 n_trials   = as_labeller(function(x) paste0("T = ", x)))) +
    scale_color_manual(values = prior_colors, name = "Prior") +
    scale_fill_manual( values = prior_colors, name = "Prior") +
    scale_linetype_manual(
      values = c("b₁ = 0.1" = "solid", "b₁ = 0.2" = "dashed"),
      name   = "Effect size"
    ) +
    labs(x = "False Positive Rate", y = "True Positive Rate",
         title = title) +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "bottom",
      strip.background = element_rect(fill = "grey92"),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(size = 12, face = "bold", hjust = 0.5)
    )
}

p_h1 <- roc_panel(h1_roc, bf3_h1, "H₁ Detection  (BF₁₀ > 3)")
p_h0 <- roc_panel(h0_roc, bf3_h0, "H₀ Detection  (BF₀₁ > 3)")

combined <- (p_h1 | p_h0) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

dir.create(here("figures"), showWarnings = FALSE)
ggsave(here("figures", "roc_dual_comparison.png"),
       combined, width = 14, height = 9, dpi = 150)
message("Saved: figures/roc_dual_comparison.png")
