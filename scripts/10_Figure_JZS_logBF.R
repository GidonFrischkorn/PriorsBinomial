# ==============================================================================
# log10 BF10 Distributions:
# Scale vs family effect — ridgeline plot across N and T conditions
# ------------------------------------------------------------------------------
# Shows the log10 BF10 density for all four prior combinations
# (Logistic/Cauchy x scale 0.25/0.707) under H0, b1=0.10, and b1=0.20.
# Facets: T (rows) x N (columns); scale on ridge y-axis; family as linetype.
# Colored tick marks at the mean of each distribution.
# Key message: within each ridge, solid (Logistic) and dashed (Cauchy) nearly
# overlap — family is negligible. Ridge rows differ substantially — scale drives
# the criterion difference.
# ==============================================================================

library(dplyr)
library(ggplot2)
library(ggridges)
library(here)

df <- readRDS(here("output", "bf_calibration_jzs_raw.rds"))

thr <- log10(3)

plot_df <- df |>
  filter(true_b1 %in% c(0.00, 0.10, 0.20)) |>
  mutate(
    # 0.25 at bottom, 0.707 at top — top ridge is more conservative
    scale_label = factor(
      paste0("Scale = ", sd_b1),
      levels = c("Scale = 0.25", "Scale = 0.707")
    ),
    b1_label = factor(
      case_when(
        true_b1 == 0.00 ~ "b₁ = 0",
        true_b1 == 0.10 ~ "b₁ = 0.10",
        true_b1 == 0.20 ~ "b₁ = 0.20"
      ),
      levels = c("b₁ = 0", "b₁ = 0.10", "b₁ = 0.20")
    ),
    prior_family = factor(
      prior,
      levels = c("logistic", "cauchy"),
      labels = c("Logistic", "Cauchy")
    ),
    n_subjects_label = factor(
      paste0("N = ", n_subjects),
      levels = paste0("N = ", sort(unique(df$n_subjects)))
    ),
    n_trials_label = factor(
      paste0("T = ", n_trials),
      levels = paste0("T = ", sort(unique(df$n_trials)))
    )
  )

# Viridis (option D, end = 0.85) to match the in-document figures:
# sequential by effect size (0 -> dark, larger -> light).
effect_colors <- c(
  "b₁ = 0"    = "#440154",
  "b₁ = 0.10" = "#277E8E",
  "b₁ = 0.20" = "#9AD93C"
)

p <- ggplot(plot_df,
            aes(x = log10_BF10, y = scale_label,
                group = interaction(scale_label, b1_label, prior_family),
                color = b1_label, linetype = prior_family)) +
  facet_grid(n_subjects_label ~ n_trials_label) +
  # Decision region shading (viridis endpoints, low alpha): evidence for H0
  # (dark/purple end) vs evidence for H1 (light/green end); grey = inconclusive.
  annotate("rect", xmin = -Inf, xmax = -thr,
           ymin = -Inf, ymax = Inf, fill = "#440154", alpha = 0.06) +
  annotate("rect", xmin = -thr, xmax = thr,
           ymin = -Inf, ymax = Inf, fill = "grey60",  alpha = 0.06) +
  annotate("rect", xmin = thr,  xmax = Inf,
           ymin = -Inf, ymax = Inf, fill = "#9AD93C", alpha = 0.10) +
  # Threshold lines
  geom_vline(xintercept =  thr, color = "grey45", linetype = "dashed",
             linewidth = 0.4) +
  geom_vline(xintercept = -thr, color = "grey45", linetype = "dashed",
             linewidth = 0.4) +
  # Ridge densities — alpha = 0 gives outlines only for clean overlap comparison
  geom_density_ridges(alpha = 0, scale = 0.80, linewidth = 0.75,
                      panel_scaling = FALSE) +
  coord_cartesian(xlim = c(-2.5, 6)) +
  scale_color_manual(values = effect_colors, name = "True effect") +
  scale_linetype_manual(
    values = c("Logistic" = "solid", "Cauchy" = "dashed"),
    name   = "Prior family"
  ) +
  labs(
    x = expression(log[10]~BF[10]),
    y = "Prior scale"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92"),
    plot.caption     = element_text(size = 9, color = "grey40", hjust = 0)
  )

dir.create(here("figures"), showWarnings = FALSE)
out <- here("figures", "logbf_distributions.png")
ggsave(out, p, width = 8, height = 12, dpi = 150)
message("Saved: figures/logbf_distributions.png")
