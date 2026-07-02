# ==============================================================================
# AUC Comparison: All four prior combinations
# Logistic/Cauchy × scale 0.25/0.707
# ------------------------------------------------------------------------------
# Computes AUC (via Wilcoxon statistic) with 95% bootstrap CI (n = 500) for the
# discriminability of the Bayes Factor between H1 (true effect) and H0 (null).
#
# A single AUC characterises the test for BOTH detection directions: because
# BF01 = 1 / BF10 is a strictly decreasing transform of BF10, ranking datasets
# by BF01 is the exact reverse of ranking them by BF10, so
#   AUC(H0 detection) = P(BF01 | H0 > BF01 | H1) = P(BF10 | H1 > BF10 | H0)
#                     = AUC(H1 detection).
# We therefore report one AUC panel rather than two analytically identical ones.
#
# Key finding (verified): AUC is invariant to prior scale and family; the small
# scale-to-scale differences are non-systematic and lie within bootstrap CIs.
# Discriminability is driven by sample size, trial count, and true effect size.
# ==============================================================================

library(dplyr)
library(ggplot2)
library(here)

df        <- readRDS(here("output", "bf_calibration_jzs_raw.rds"))
set.seed(42)
n_boot    <- 500
jzs_scale <- round(sqrt(2) / 2, 3)

make_prior_label <- function(dist, sc) {
  sprintf("%s(0, %s)", ifelse(dist == "logistic", "Logistic", "Cauchy"), sc)
}

prior_levels <- c(
  "Logistic(0, 0.25)", "Cauchy(0, 0.25)",
  "Logistic(0, 0.707)", "Cauchy(0, 0.707)"
)

compute_auc <- function(signal_vals, noise_vals) {
  W <- wilcox.test(signal_vals, noise_vals, alternative = "greater")$statistic
  as.numeric(W) / (length(signal_vals) * length(noise_vals))
}

# --- Compute AUC (H1 vs H0 discriminability) ----------------------------------

auc_rows <- list()
idx      <- 1

for (sc in c(0.25, jzs_scale)) {
  for (dist in c("logistic", "cauchy")) {
    for (ns in c(20, 30, 50)) {
      for (nt in c(20, 50)) {
        h0_bf10 <- df$BF10[df$prior == dist & df$sd_b1 == sc &
                            df$n_subjects == ns & df$n_trials == nt & df$true_b1 == 0]
        for (tb1 in c(0.1, 0.2)) {
          h1_bf10 <- df$BF10[df$prior == dist & df$sd_b1 == sc &
                              df$n_subjects == ns & df$n_trials == nt & df$true_b1 == tb1]
          if (length(h0_bf10) == 0 || length(h1_bf10) == 0) next

          sv <- h1_bf10; nv <- h0_bf10          # signal = H1, noise = H0
          auc_obs  <- compute_auc(sv, nv)
          auc_boot <- vapply(seq_len(n_boot), function(b) {
            compute_auc(sample(sv, replace = TRUE), sample(nv, replace = TRUE))
          }, numeric(1))

          auc_rows[[idx]] <- data.frame(
            prior      = dist, sd_b1 = sc,
            n_subjects = ns,  n_trials = nt,
            true_b1    = tb1,
            auc        = auc_obs,
            auc_lo     = quantile(auc_boot, 0.025),
            auc_hi     = quantile(auc_boot, 0.975)
          )
          idx <- idx + 1
        }
      }
    }
  }
}

auc_df <- bind_rows(auc_rows) |>
  mutate(
    prior_label = factor(make_prior_label(prior, sd_b1), levels = prior_levels),
    scale_label = factor(paste0("Scale = ", sd_b1),
                         levels = c("Scale = 0.25", "Scale = 0.707")),
    family_label = factor(ifelse(prior == "logistic", "Logistic", "Cauchy"),
                          levels = c("Logistic", "Cauchy")),
    n_subjects     = factor(n_subjects, levels = c(20, 30, 50)),
    b1_label       = factor(paste0("b₁ = ", true_b1),
                            levels = c("b₁ = 0.1", "b₁ = 0.2")),
    n_trials_label = factor(paste0("T = ", n_trials),
                            levels = c("T = 20", "T = 50"))
  )

# Viridis (option D, end = 0.85) to match the in-document figures:
# scale 0.25 -> dark end, scale 0.707 -> light end.
scale_colors <- c("Scale = 0.25" = "#440154", "Scale = 0.707" = "#9AD93C")
dodge        <- position_dodge(width = 0.35)

p_auc <- ggplot(auc_df,
                aes(x = n_subjects, y = auc,
                    color = scale_label, linetype = family_label,
                    group = interaction(scale_label, family_label))) +
  facet_grid(b1_label ~ n_trials_label) +
  geom_hline(yintercept = 0.5, color = "grey60", linetype = "dashed",
             linewidth = 0.4) +
  geom_linerange(aes(ymin = auc_lo, ymax = auc_hi),
                 position = dodge, linewidth = 0.6) +
  geom_line(position = dodge, linewidth = 0.75) +
  geom_point(position = dodge, size = 2.5) +
  scale_color_manual(values = scale_colors, name = "Prior scale") +
  scale_linetype_manual(
    values = c("Logistic" = "solid", "Cauchy" = "dashed"),
    name   = "Prior family"
  ) +
  scale_y_continuous("AUC", limits = c(0.45, 1.0),
                     breaks = seq(0.5, 1.0, by = 0.1)) +
  scale_x_discrete("N subjects") +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey92"),
    panel.grid.minor = element_blank()
  )

dir.create(here("figures"), showWarnings = FALSE)
ggsave(here("figures", "auc_jzs_comparison.png"),
       p_auc, width = 8, height = 6, dpi = 150)
message("Saved: figures/auc_jzs_comparison.png")
