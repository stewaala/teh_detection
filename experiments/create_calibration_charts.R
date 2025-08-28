source("simulate_and_fit_forests.R")

suppressPackageStartupMessages({
  library(dplyr); library(purrr); library(ggplot2); library(tibble)
})

dir.create("charts", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Common settings
# ------------------------------------------------------------
ns   <- c(100, 500, 1000, 2500)
rho0 <- 0.5

# ----------------------------------------------------------------------
# Robust loader: works whether or not there is a Z-suffix bundle on disk
# ----------------------------------------------------------------------
load_bundle_match <- function(n, rho, link_type, p_base = NULL, Z = NULL,
                              data_dir = "data") {
  # Try the standard loader first (uses your make_base_filename internally)
  b <- try(load_bundle(n = n, rho = rho, link_type = link_type, p_base = p_base),
           silent = TRUE)
  if (!inherits(b, "try-error")) return(b)
  
  # Fallback: construct the filename ourselves and optionally append _Z###
  baseline_tag <- if (is.null(p_base)) "defaultBaseline" else {
    bstr <- formatC(p_base, format = "f", digits = 3)
    bstr <- sub("0+$", "", sub("\\.$", "", bstr))
    paste0("baseline", bstr)
  }
  z_tag <- if (is.null(Z)) "" else sprintf("_Z%03d", Z)
  fname <- sprintf("bundle_n%05d_rho%02d_%s_%s%s.rds",
                   n, round(100 * rho), link_type, baseline_tag, z_tag)
  path <- file.path(data_dir, fname)
  if (!file.exists(path))
    stop("Cannot find bundle file via either route: ", path)
  readRDS(path)
}

# ----------------------------------------------------------------------
# Internal: extract (trial_id, pred, truth) as a tibble
# pred_col ∈ {crte_hat_grf, crte_hat_glm, cate_hat_grf, cate_hat_glm}
# truth_col ∈ {CRTE, CATE}; transform() may be log, identity, ...
# ----------------------------------------------------------------------
.extract_pred_truth <- function(bundle, pred_col, truth_col, transform = identity) {
  stopifnot(all(c("trial_id", pred_col, truth_col) %in% names(bundle)))
  d <- as.data.frame(bundle[, c("trial_id", pred_col, truth_col), with = FALSE])
  names(d) <- c("trial_id", "pred", "truth")
  d$pred  <- transform(d$pred)
  d$truth <- transform(d$truth)
  d <- d[is.finite(d$pred) & is.finite(d$truth), , drop = FALSE]
  tibble::as_tibble(d)
}

# ----------------------------------------------------------------------
# One bundle -> per-trial decile reliability + trialwise slope/intercept + ECE
# ECE here is weighted L1 gap across decile bins (equal-count bins per trial)
# ----------------------------------------------------------------------
calibration_from_bundle <- function(bundle, pred_col, truth_col,
                                    transform = identity, n_bins = 10) {
  d <- .extract_pred_truth(bundle, pred_col, truth_col, transform)
  
  # Bin by predicted value within each trial (equal-count deciles)
  d <- d %>%
    dplyr::group_by(trial_id) %>%
    dplyr::mutate(bin = dplyr::ntile(pred, n_bins)) %>%
    dplyr::ungroup()
  
  # Per-trial, per-bin means (reliability curve points)
  by_bin <- d %>%
    dplyr::group_by(trial_id, bin) %>%
    dplyr::summarise(pred_mean = mean(pred),
                     obs_mean  = mean(truth),
                     n         = dplyr::n(),
                     .groups = "drop")
  
  # Per-trial slope/intercept of (truth ~ pred)
  by_trial_lm <- d %>%
    dplyr::group_by(trial_id) %>%
    dplyr::group_modify(~{
      if (var(.x$pred) == 0) return(tibble::tibble(intercept = NA_real_, slope = NA_real_))
      cf <- coef(stats::lm(truth ~ pred, data = .x))
      tibble::tibble(intercept = unname(cf[1]), slope = unname(cf[2]))
    }) %>%
    dplyr::ungroup()
  
  # Per-trial ECE (weighted absolute gap between obs and pred means over bins)
  by_trial_ece <- by_bin %>%
    dplyr::group_by(trial_id) %>%
    dplyr::summarise(ece = sum(abs(obs_mean - pred_mean) * n / sum(n)),
                     .groups = "drop")
  
  list(
    by_bin   = by_bin,
    by_trial = dplyr::left_join(by_trial_lm, by_trial_ece, by = "trial_id")
  )
}

# ----------------------------------------------------------------------
# Aggregate across trials -> mean curve + SE, and summary metrics + SE
# ----------------------------------------------------------------------
aggregate_calibration <- function(calib) {
  curve <- calib$by_bin %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(pred = mean(pred_mean),
                     obs  = mean(obs_mean),
                     se   = stats::sd(obs_mean) / sqrt(dplyr::n()),
                     .groups = "drop")
  
  mets <- calib$by_trial
  metrics <- tibble::tibble(
    slope      = mean(mets$slope, na.rm = TRUE),
    slope_se   = stats::sd(mets$slope, na.rm = TRUE) / sqrt(sum(is.finite(mets$slope))),
    intercept  = mean(mets$intercept, na.rm = TRUE),
    int_se     = stats::sd(mets$intercept, na.rm = TRUE) / sqrt(sum(is.finite(mets$intercept))),
    ece        = mean(mets$ece, na.rm = TRUE),
    ece_se     = stats::sd(mets$ece, na.rm = TRUE) / sqrt(sum(is.finite(mets$ece)))
  )
  list(curve = curve, metrics = metrics)
}

# ----------------------------------------------------------------------
# Convenience: build an aggregated curve tibble for (bundle, method, target)
# target ∈ {"crte","cate"}; scale controls transform (e.g., "log" or "identity")
# ----------------------------------------------------------------------
curve_for <- function(bundle, method = c("GRF","RRCF"),
                      target = c("crte","cate"),
                      scale  = c("log","identity"),
                      n_bins = 10) {
  method <- match.arg(method)
  target <- match.arg(target)
  scale  <- match.arg(scale)
  
  if (target == "crte") {
    pred_col <- if (method == "GRF") "crte_hat_grf" else "crte_hat_glm"
    truth_col <- "CRTE"
    transform <- if (scale == "log") log else identity
  } else {
    pred_col <- if (method == "GRF") "cate_hat_grf" else "cate_hat_glm"
    truth_col <- "CATE"
    transform <- identity   # CATE stays on additive scale
  }
  
  agg <- aggregate_calibration(
    calibration_from_bundle(bundle, pred_col, truth_col,
                            transform = transform, n_bins = n_bins)
  )
  dplyr::mutate(agg$curve, method = method)
}

# ----------------------------------------------------------------------
# Plot observed vs predicted with 45° line; facet by any column name
# df_curve must have columns: pred, obs, se, method, and facet columns
# ----------------------------------------------------------------------
plot_calibration <- function(df_curve, facet_by = NULL, xlab, title = "",
                             save = NULL) {
  p <- ggplot2::ggplot(df_curve,
                       ggplot2::aes(x = pred, y = obs, color = method, linetype = method)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.4,
                         linetype = "dotted", color = "grey40") +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = obs - se, ymax = obs + se),
                           width = 0.0, linewidth = 0.35, alpha = 0.7) +
    ggplot2::labs(x = xlab, y = "Empirical mean (truth)", color = "Method",
                  linetype = "Method", title = title) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::scale_color_manual(values = c("RRCF" = "steelblue",
                                           "GRF"  = "forestgreen"))
  
  if (!is.null(facet_by)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)))
  }
  
  if (is.list(save)) {
    dir.create(save$dir %||% "charts", recursive = TRUE, showWarnings = FALSE)
    ts_fmt <- save$timestamp_fmt %||% "%Y%m%d%H%M"
    stamp  <- format(Sys.time(), ts_fmt)
    base   <- if (is.null(save$name_components)) "calibration" else
      paste(save$name_components, collapse = "_")
    outfile <- file.path(save$dir %||% "charts", paste0(base, "_", stamp, ".pdf"))
    ggplot2::ggsave(outfile, plot = p,
                    width  = save$width_in  %||% 10.5,
                    height = save$height_in %||% 5.5,
                    units  = "in")
    message("Saved plot to: ", outfile)
  }
  invisible(p)
}

rhos      <- c(0.25, 0.5, 0.75)
n_default <- 2500

# ======================================================================
# (1) MAPE — single chart log link
# ======================================================================
df_cal_log <- purrr::map_dfr(rhos, function(rho) {
  b <- load_bundle_match(n = n_default, rho = rho, link_type = "log")
  dplyr::bind_rows(
    dplyr::mutate(curve_for(b, "RRCF", target = "crte", scale = "log"), rho = rho),
    dplyr::mutate(curve_for(b, "GRF",  target = "crte", scale = "log"), rho = rho)
  )
})

plot_calibration(
  df_cal_log, facet_by = "rho", xlab = "Predicted log CRTE",
  title = "",
  save = list(dir = "charts",
              name_components = c("calib", "log", "crte_by_rho"),
              timestamp_fmt = "%Y%m%d%H%M",
              width_in = 11, height_in = 5.5)
)
# ======================================================================
# (2) MAE — single chart log link
# ======================================================================
# --- CATE (identity link) calibration by rho, n = 2500, default baseline, Z = 0 ---
df_cal_cate <- purrr::map_dfr(rhos, function(rho) {
  b <- load_bundle_match(n = n_default, rho = rho, link_type = "identity")
  dplyr::bind_rows(
    dplyr::mutate(curve_for(b, method = "RRCF", target = "cate", scale = "identity"), rho = rho),
    dplyr::mutate(curve_for(b, method = "GRF",  target = "cate", scale = "identity"), rho = rho)
  )
})

plot_calibration(
  df_curve = df_cal_cate %>% dplyr::mutate(rho = factor(rho, levels = rhos)),
  facet_by = "rho",
  xlab     = "Predicted CATE",
  title    = "",
  save     = list(
    dir            = "charts",
    name_components= c("calib","identity","cate_by_rho"),
    timestamp_fmt  = "%Y%m%d%H%M",
    width_in       = 11,
    height_in      = 5.5
  )
)

# ------------------------------------------------------------
# A) LOG-LINK DGP: Calibration of log-CRTE by sample size (rho = 0.5)
# ------------------------------------------------------------
df_small_log <- purrr::map_dfr(ns, function(n) {
  b  <- load_bundle_match(n = n, rho = rho0, link_type = "log")
  nb <- if (n == 100) 5 else 10  # fewer bins for very small n
  bind_rows(
    mutate(curve_for(b, method = "RRCF", target = "crte", scale = "log",      n_bins = nb),
           panel = paste0("n = ", n)),
    mutate(curve_for(b, method = "GRF",  target = "crte", scale = "log",      n_bins = nb),
           panel = paste0("n = ", n))
  )
})

plot_calibration(
  df_curve = df_small_log %>% mutate(panel = factor(panel, levels = paste0("n = ", ns))),
  facet_by = "panel",
  xlab     = "Predicted log CRTE",
  title    = "",
  save     = list(dir = "charts",
                  name_components = c("calib","log","crte_by_n_rho50"),
                  timestamp_fmt   = "%Y%m%d%H%M",
                  width_in        = 12, height_in = 6)
)

# ------------------------------------------------------------
# B) IDENTITY-LINK DGP: Calibration of CATE by sample size (rho = 0.5)
# ------------------------------------------------------------
df_small_id <- purrr::map_dfr(ns, function(n) {
  b  <- load_bundle_match(n = n, rho = rho0, link_type = "identity")
  nb <- if (n == 100) 5 else 10
  bind_rows(
    mutate(curve_for(b, method = "RRCF", target = "cate", scale = "identity", n_bins = nb),
           panel = paste0("n = ", n)),
    mutate(curve_for(b, method = "GRF",  target = "cate", scale = "identity", n_bins = nb),
           panel = paste0("n = ", n))
  )
})

plot_calibration(
  df_curve = df_small_id %>% mutate(panel = factor(panel, levels = paste0("n = ", ns))),
  facet_by = "panel",
  xlab     = "Predicted CATE",
  title    = "",
  save     = list(dir = "charts",
                  name_components = c("calib","identity","cate_by_n_rho50"),
                  timestamp_fmt   = "%Y%m%d%H%M",
                  width_in        = 12, height_in = 6)
)

# ------------------------------------------------------------
# C/D) Summary metrics vs n (slope and ECE) for each DGP
# ------------------------------------------------------------

# Helper: compute metrics (slope, intercept, ECE) for a bundle/method/target
metrics_for <- function(bundle, method, target, scale, n_bins = 10) {
  if (target == "crte") {
    pred_col <- if (method == "RRCF") "crte_hat_glm" else "crte_hat_grf"
    truth_col <- "CRTE"
    transform <- log
  } else {
    pred_col <- if (method == "RRCF") "cate_hat_glm" else "cate_hat_grf"
    truth_col <- "CATE"
    transform <- identity
  }
  agg <- aggregate_calibration(
    calibration_from_bundle(bundle, pred_col, truth_col,
                            transform = transform, n_bins = n_bins)
  )
  cbind(agg$metrics, method = method)
}

# ---- Log-link metrics vs n (rho = 0.5) ----
df_mets_log <- purrr::map_dfr(ns, function(n) {
  b  <- load_bundle_match(n = n, rho = rho0, link_type = "log")
  nb <- if (n == 100) 5 else 10
  bind_rows(
    mutate(metrics_for(b, "RRCF", "crte", "log", n_bins = nb), n = n),
    mutate(metrics_for(b, "GRF",  "crte", "log", n_bins = nb), n = n)
  )
})

df_mets_log_long <- bind_rows(
  transmute(df_mets_log, n, method, metric = "Slope", value = slope, se = slope_se),
  transmute(df_mets_log, n, method, metric = "ECE",   value = ece,   se = ece_se)
)

p_log_mets <- ggplot(df_mets_log_long,
                     aes(x = n, y = value, color = method, group = method)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 30, alpha = 0.7) +
  scale_x_continuous(breaks = ns) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = "Sample size n", y = NULL, color = "Method",
       title = "") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank()) +
  scale_color_manual(values = c("RRCF" = "steelblue", "GRF" = "forestgreen"))

ggsave(file.path("charts",
                 paste0("calib_metrics_log_crte_by_n_rho50_",
                        format(Sys.time(), "%Y%m%d%H%M"), ".pdf")),
       p_log_mets, width = 10, height = 5.5, units = "in")

# ---- Identity-link metrics vs n (rho = 0.5) ----
df_mets_id <- purrr::map_dfr(ns, function(n) {
  b  <- load_bundle_match(n = n, rho = rho0, link_type = "identity")
  nb <- if (n == 100) 5 else 10
  bind_rows(
    mutate(metrics_for(b, "RRCF", "cate", "identity", n_bins = nb), n = n),
    mutate(metrics_for(b, "GRF",  "cate", "identity", n_bins = nb), n = n)
  )
})

df_mets_id_long <- bind_rows(
  transmute(df_mets_id, n, method, metric = "Slope", value = slope, se = slope_se),
  transmute(df_mets_id, n, method, metric = "ECE",   value = ece,   se = ece_se)
)

p_id_mets <- ggplot(df_mets_id_long,
                    aes(x = n, y = value, color = method, group = method)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 30, alpha = 0.7) +
  scale_x_continuous(breaks = ns) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = "Sample size n", y = NULL, color = "Method",
       title = "") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank()) +
  scale_color_manual(values = c("RRCF" = "steelblue", "GRF" = "forestgreen"))

ggsave(file.path("charts",
                 paste0("calib_metrics_identity_cate_by_n_rho50_",
                        format(Sys.time(), "%Y%m%d%H%M"), ".pdf")),
       p_id_mets, width = 10, height = 5.5, units = "in")

