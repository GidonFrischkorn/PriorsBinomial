# plotting.R
# ggplot2 figure functions for the prior predictive analysis.
# All figures use theme_minimal() + viridis color palette.
#
# Required packages: ggplot2, ggridges, viridis, dplyr, tidyr

# ==============================================================================
# Shared theme and helpers
# ==============================================================================

#' Base ggplot2 theme for all figures
prior_theme <- function() {
  ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      panel.grid.minor  = ggplot2::element_blank(),
      strip.text        = ggplot2::element_text(face = "bold"),
      legend.position   = "bottom",
      legend.title      = ggplot2::element_text(face = "bold")
    )
}

# labels for facets and legends
link_labels    <- c(logit = "Logit link", probit = "Probit link")
dist_labels    <- c(cauchy = "Cauchy", normal = "Normal", logistic = "Logistic")


# ==============================================================================
# Figure 1: Density ridges of p_intercept by sd_b0
# ==============================================================================

plot_p_intercept_ridges <- function(draws_long, sd_b1_focus = 0.25) {
  ridge_scale <- 1.8
  dist_b0_levels <- c("Cauchy", "Normal", "Logistic")

  draws_long <- draws_long |>
    dplyr::mutate(
      sd_b0_label = paste0("sd[b[0]] == ", sd_b0),
      dist_b0     = factor(dist_b0, levels = c("cauchy", "normal", "logistic"),
                           labels = dist_b0_levels)
    )

  bin_breaks <- seq(0, 1, length.out = 81)  # 80 bins
  panel_max <- draws_long |>
    dplyr::group_by(sd_b0_label, link, dist_b0) |>
    dplyr::summarise(
      ridge_max = max(hist(p_intercept, breaks = bin_breaks, plot = FALSE)$density),
      .groups = "drop"
    ) |>
    dplyr::group_by(sd_b0_label, link) |>
    dplyr::summarise(
      hist_max = max(ridge_max),
      .groups = "drop"
    )

  ref_data <- expand.grid(
    p       = seq(0.001, 0.999, length.out = 500),
    sd_b0   = unique(draws_long$sd_b0),
    link    = c("logit", "probit"),
    dist_b0 = dist_b0_levels,
    stringsAsFactors = FALSE
  )
  ref_data$density <- mapply(function(p, sd, lnk) {
    if (lnk == "logit")  d_logistic_on_p(p, scale = sd)
    else                  d_normal_on_p(p, sigma = sd)
  }, ref_data$p, ref_data$sd_b0, ref_data$link)
  ref_data$sd_b0_label <- paste0("sd[b[0]] == ", ref_data$sd_b0)
  ref_data$dist_b0 <- factor(ref_data$dist_b0, levels = dist_b0_levels)

  ref_data <- ref_data |>
    dplyr::left_join(panel_max, by = c("sd_b0_label", "link")) |>
    dplyr::mutate(
      y_pos = as.numeric(dist_b0) + (density / hist_max) * ridge_scale
    )

  ggplot2::ggplot(
    draws_long,
    ggplot2::aes(x = p_intercept, y = dist_b0, fill = dist_b0, colour = dist_b0)
  ) +
    ggridges::geom_density_ridges(stat = "binline", bins = 80,
                                   alpha = 0.4, scale = ridge_scale) +
    ggplot2::geom_line(
      data = ref_data,
      ggplot2::aes(x = p, y = y_pos,
                   group = interaction(dist_b0, sd_b0_label, link)),
      linetype = "dashed", linewidth = 0.6, colour = "black",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_vline(xintercept = c(0.05, 0.95), linetype = "dashed",
                        linewidth = 0.4, colour = "grey40") +
    ggplot2::scale_x_continuous(
      expression("Success probability" ~ italic(p)),
      limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0", "0.25", "0.50", "0.75", "1"),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::scale_y_discrete("Prior family",
                               expand = ggplot2::expansion(add = c(0.3, 1.5))) +
    ggplot2::scale_fill_viridis_d("Prior family", option = "D", end = 0.85) +
    ggplot2::scale_colour_viridis_d("Prior family", option = "D", end = 0.85) +
    ggplot2::facet_grid(
      sd_b0_label ~ link,
      labeller = ggplot2::labeller(
        link = link_labels,
        sd_b0_label = ggplot2::label_parsed
      )
    ) +
    prior_theme() +
    ggplot2::coord_cartesian(clip = "off")
}


# ==============================================================================
# Figure 2: Density ridges of abs_delta_p by sd_b1
# ==============================================================================

plot_delta_p_ridges <- function(draws_long, sd_b0_focus = 0.75) {
  draws_long$dist_b1 <- factor(draws_long$dist_b1,
    levels = c("cauchy", "normal", "logistic"),
    labels = c("Cauchy", "Normal", "Logistic"))
  draws_long$sd_b1_label <- paste0("sd[b[1]] == ", draws_long$sd_b1)
  # Order facets by numeric sd_b1 value
  sd_b1_order <- sort(unique(as.numeric(as.character(
    factor(draws_long$sd_b1)
  ))))
  draws_long$sd_b1_label <- factor(
    draws_long$sd_b1_label,
    levels = paste0("sd[b[1]] == ", sd_b1_order)
  )

  ggplot2::ggplot(
    draws_long,
    ggplot2::aes(x = delta_p, y = dist_b1, fill = dist_b1, color = dist_b1)
  ) +
    ggridges::geom_density_ridges(stat = "binline", bins = 80,
                                   alpha = 0.4, scale = 0.9) +
    ggplot2::geom_vline(xintercept = 0,
                        linetype = "solid", linewidth = 0.5, color = "grey20") +
    ggplot2::geom_vline(xintercept = c(-0.30, -0.20, -0.10, 0.10, 0.20, 0.30),
                        linetype = "dotted", linewidth = 0.5, color = "grey40") +
    ggplot2::scale_x_continuous(
      expression(delta * italic(p)),
      limits = c(-1, 1),
      breaks = seq(-1, 1, 0.25)
    ) +
    ggplot2::scale_y_discrete("Prior family") +
    ggplot2::scale_fill_viridis_d("Prior family", option = "D", end = 0.85) +
    ggplot2::scale_color_viridis_d("Prior family", option = "D", end = 0.85) +
    ggplot2::facet_wrap(~ sd_b1_label, ncol = 3,
                        labeller = ggplot2::label_parsed) +
    prior_theme() +
    ggplot2::theme(
      legend.position = "inside",
      legend.position.inside = c(0.83, 0.08),
      legend.justification = c(0.5, 0.5),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key.size = ggplot2::unit(0.8, "cm"),
      axis.title.x = ggplot2::element_text(hjust = 0.17)
    )
}


# ==============================================================================
# Figure 3: Floor/ceiling heatmap
# ==============================================================================

plot_floor_ceiling_heatmap <- function(summaries, threshold = 0.10) {
  dat <- summaries |>
    dplyr::summarise(
      prob_floor_ceiling = mean(prob_floor_ceiling),
      .by = c(link, dist_b0, sd_b0)
    )
  dat$sd_b0   <- factor(dat$sd_b0)
  dat$dist_b0 <- factor(dat$dist_b0, levels = c("cauchy", "normal", "logistic"),
                        labels = c("Cauchy", "Normal", "Logistic"))
  fill_mid <- max(dat$prob_floor_ceiling, 0.20) / 2
  dat$text_color <- ifelse(dat$prob_floor_ceiling > fill_mid, "white", "black")

  ggplot2::ggplot(dat, ggplot2::aes(x = sd_b0, y = dist_b0, fill = prob_floor_ceiling)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", prob_floor_ceiling), color = text_color),
      size = 3.5, fontface = "bold"
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_viridis_c(
      expression(P(italic(p) < 0.05 ~ "or" ~ italic(p) > 0.95)),
      option = "C", direction = -1,
      limits = c(0, max(dat$prob_floor_ceiling, 0.20))
    ) +
    ggplot2::scale_x_discrete(expression(sd[b[0]])) +
    ggplot2::scale_y_discrete("Intercept prior family") +
    ggplot2::facet_wrap(~link, labeller = ggplot2::labeller(link = link_labels)) +
    prior_theme() +
    ggplot2::theme(legend.position = "right")
}


# ==============================================================================
# Figure 4: Quantile profile of abs_delta_p vs sd_b1
# ==============================================================================

plot_quantile_profile <- function(summaries, sd_b0_focus = 0.75) {
  dat <- summaries |>
    dplyr::filter(
      sd_b0 == sd_b0_focus,
      (link == "logit"  & dist_b0 == "logistic") |
      (link == "probit" & dist_b0 == "normal")
    ) |>
    tidyr::pivot_longer(
      cols      = c(adp_q50, adp_q90, adp_q95),
      names_to  = "quantile",
      values_to = "value"
    ) |>
    dplyr::mutate(
      quantile = dplyr::recode(quantile,
        adp_q50 = "50th percentile",
        adp_q90 = "90th percentile",
        adp_q95 = "95th percentile"
      ),
      dist_b1 = factor(dist_b1, levels = c("cauchy", "normal", "logistic"),
                       labels = c("Cauchy", "Normal", "Logistic"))
    )

  ggplot2::ggplot(dat, ggplot2::aes(x = sd_b1, y = value, color = dist_b1,
                                     linetype = quantile)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = c(0.10, 0.20, 0.30),
                        linetype = "dotted", linewidth = 0.4, color = "grey60") +
    ggplot2::scale_x_continuous(expression(sd[b[1]])) +
    ggplot2::scale_y_continuous(
      expression("Quantile of" ~ group("|", delta * italic(p), "|")),
      limits = c(0, 1)
    ) +
    ggplot2::scale_color_viridis_d("Prior family", option = "D", end = 0.85) +
    ggplot2::scale_linetype_manual("Quantile",
      values = c("50th percentile" = "solid",
                 "90th percentile" = "dashed",
                 "95th percentile" = "dotted")) +
    ggplot2::facet_wrap(~link, labeller = ggplot2::labeller(link = link_labels)) +
    prior_theme()
}


# ==============================================================================
# Figure 5: Logit–probit equivalence
# ==============================================================================

plot_link_equivalence <- function(summaries, sd_b0_focus = 0.75) {
  dat <- summaries |>
    dplyr::filter(sd_b0 == sd_b0_focus, dist_b0 %in% c("logistic", "normal")) |>
    dplyr::mutate(
      matched = (link == "logit" & dist_b0 == "logistic") |
                (link == "probit" & dist_b0 == "normal"),
      dist_b1 = factor(dist_b1, levels = c("cauchy", "normal", "logistic"),
                       labels = c("Cauchy", "Normal", "Logistic"))
    ) |>
    dplyr::filter(matched)

  ggplot2::ggplot(dat, ggplot2::aes(x = sd_b1, y = adp_q90, color = link,
                                     linetype = dist_b1)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = c(0.10, 0.20, 0.30),
                        linetype = "dotted", linewidth = 0.4, color = "grey60") +
    ggplot2::scale_x_continuous(expression(sd[b[1]])) +
    ggplot2::scale_y_continuous(
      expression("90th percentile of" ~ group("|", delta * italic(p), "|")),
      limits = c(0, 1)
    ) +
    ggplot2::scale_color_viridis_d("Link function",
                                    labels = link_labels, option = "D", end = 0.85) +
    ggplot2::scale_linetype_discrete("Prior family") +
    prior_theme()
}


# ==============================================================================
# Figure 6: Misfit heatmap
# ==============================================================================

plot_misfit_heatmap <- function(summaries) {
  dat <- summaries |>
    dplyr::summarise(
      prob_floor_ceiling = mean(prob_floor_ceiling),
      .by = c(link, dist_b0, sd_b0)
    ) |>
    dplyr::mutate(
      matched = (link == "logit" & dist_b0 == "logistic") |
                (link == "probit" & dist_b0 == "normal"),
      label = paste0(
        dplyr::recode(link, logit = "Logit", probit = "Probit"),
        " + ",
        dplyr::recode(dist_b0, normal = "Normal", logistic = "Logistic", cauchy = "Cauchy")
      )
    )
  fill_mid <- max(dat$prob_floor_ceiling) / 2
  dat$text_color <- ifelse(dat$prob_floor_ceiling > fill_mid, "white", "black")

  ggplot2::ggplot(dat, ggplot2::aes(x = factor(sd_b0), y = label,
                                     fill = prob_floor_ceiling)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", prob_floor_ceiling), color = text_color),
      size = 3.5, fontface = "bold"
    ) +
    ggplot2::geom_tile(
      data = dplyr::filter(dat, matched),
      color = "#FFD700", linewidth = 1.5, fill = NA
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_viridis_c(
      expression(P(italic(p) < 0.05 ~ "or" ~ italic(p) > 0.95)),
      option = "C", direction = -1
    ) +
    ggplot2::scale_x_discrete(expression(sd[b[0]])) +
    ggplot2::scale_y_discrete("Link + prior family") +
    prior_theme() +
    ggplot2::theme(legend.position = "right")
}


# ==============================================================================
# Figure S1: BF Validation scatter plot (Savage-Dickey vs Bridge Sampling)
# ==============================================================================

#' Scatter plot comparing Savage-Dickey and Bridge Sampling Bayes Factors
#'
#' Each point is one simulation replication. The identity line marks perfect
#' agreement between methods. Per-panel Pearson correlations are annotated.
#'
#' @param bf_data Data frame with columns: log10_SD, log10_BS, dist_b1,
#'   true_b1_label, sd_b1_label (one row per replication).
#' @return A ggplot2 object.
plot_bf_validation_scatter <- function(bf_data) {

  # per-panel correlation annotation
  cor_df <- bf_data |>
    dplyr::summarise(
      r     = stats::cor(log10_SD, log10_BS, use = "complete.obs"),
      label = sprintf("italic(r) == %.2f", r),
      .by   = c(sd_b1_label, true_b1_label)
    )

  ggplot2::ggplot(bf_data, ggplot2::aes(x = log10_SD, y = log10_BS,
                                         colour = dist_b1)) +
    ggplot2::geom_abline(intercept = 0, slope = 1,
                         colour = "grey50", linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = c(log10(3), log10(1/3), 0),
                        linetype = "dashed", colour = "grey60",
                        linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = c(log10(3), log10(1/3), 0),
                        linetype = "dashed", colour = "grey60",
                        linewidth = 0.4) +
    ggplot2::geom_point(alpha = 0.5, size = 1.5) +
    ggplot2::geom_text(
      data    = cor_df,
      ggplot2::aes(x = -Inf, y = Inf, label = label),
      hjust = -0.1, vjust = 1.5, size = 4, colour = "black",
      parse = TRUE, inherit.aes = FALSE
    ) +
    ggplot2::facet_grid(
      sd_b1_label ~ true_b1_label,
      labeller = ggplot2::label_parsed,
      scales   = "free"
    ) +
    ggplot2::scale_colour_viridis_d("Prior family", option = "D", end = 0.85) +
    ggplot2::labs(
      x = expression(log[10](BF[10]^{SD})),
      y = expression(log[10](BF[10]^{BS}))
    ) +
    prior_theme()
}


# ==============================================================================
# Figure S2: Bland-Altman plot (method discrepancy vs evidence strength)
# ==============================================================================

#' Bland-Altman plot for Savage-Dickey vs Bridge Sampling Bayes Factors
#'
#' X-axis shows average evidence strength, y-axis shows the discrepancy
#' (SD minus BS on log10 scale). Per-panel mean bias is shown as a dashed
#' red line.
#'
#' @param bf_data Data frame with columns: mean_log10, diff_log10, dist_b1,
#'   sd_b1_factor, true_b1_label (one row per replication).
#' @return A ggplot2 object.
plot_bf_validation_bland_altman <- function(bf_data) {

  # per-panel mean bias
  bias_df <- bf_data |>
    dplyr::summarise(
      mean_bias = mean(diff_log10, na.rm = TRUE),
      .by = true_b1_label
    )

  ggplot2::ggplot(bf_data, ggplot2::aes(x = mean_log10, y = diff_log10,
                                         colour = dist_b1,
                                         shape  = sd_b1_factor)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.5) +
    ggplot2::geom_hline(
      data     = bias_df,
      ggplot2::aes(yintercept = mean_bias),
      linetype = "dashed", colour = "red", linewidth = 0.4,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(alpha = 0.5, size = 1.5) +
    ggplot2::facet_wrap(
      ~ true_b1_label, ncol = 3,
      labeller = ggplot2::label_parsed,
      scales   = "free"
    ) +
    ggplot2::scale_colour_viridis_d("Prior family", option = "D", end = 0.85) +
    ggplot2::scale_shape_manual(
      expression(sd[b[1]]),
      values = c("0.25" = 16, "0.5" = 17)
    ) +
    ggplot2::labs(
      x = expression((log[10](BF[10]^{SD}) + log[10](BF[10]^{BS})) / 2),
      y = expression(log[10](BF[10]^{SD}) - log[10](BF[10]^{BS}))
    ) +
    prior_theme()
}
