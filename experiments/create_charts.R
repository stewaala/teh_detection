############################################################
## Master script: MAPE/MAE charts with 5–95% CIs
## Uses a single plotting function, called by six scenarios:
##  1) MAPE — two panels: log (left), identity (right)
##  2) MAE  — two panels: log (left), identity (right)
##  3) MAPE by sample size (log link) — panels: n in {100, 500, 1000, 2500}
##  4) MAE  by sample size (log link) — panels: n in {100, 500, 1000, 2500}
##  5) MAPE by baseline (log link)    — panels: default, 0.05, 0.01
##  6) MAE  by baseline (log link)    — panels: default, 0.05, 0.01
##
## Requires in scope:
##   - load_bundle()
##   - mape_from_bundle_ci(bundle, forest_type=c("grf","rrcf"))
##   - mae_from_bundle_ci(bundle,  forest_type=c("grf","rrcf"))
## If needed:
# source("simulate_and_fit_forests.R")
############################################################

source("simulate_and_fit_forests.R")

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)

`%||%` <- function(a, b) if (!is.null(a)) a else b

plot_with_ci <- function(data,
                         facet_by,
                         metric = c("mape","mae"),
                         title = NULL,
                         subtitle = NULL,
                         share_y = FALSE,
                         save = NULL,   # list(dir="charts", name_components=c(...), timestamp_fmt="%Y%m%d%H%M", width_in=11, height_in=6)
                         style = list(),
                         facet_nrow = 1,          # <— NEW
                         facet_ncol = NULL) {
  
  metric <- match.arg(metric)
  stopifnot(is.data.frame(data),
            all(c("rho","method","estimate","ci_lower","ci_upper") %in% names(data)))
  
  defaults <- list(
    colors       = c("RRCF" = "steelblue", "GRF" = "forestgreen"),
    linetypes    = c("RRCF" = "solid",     "GRF" = "dashed"),  # <— NEW
    line_width   = 1.1,
    point_size   = 2.6,
    ci_linewidth = 0.35,
    ci_cap_width = 0.018,
    ci_alpha     = 0.7,
    theme_base   = 12,
    pad_x        = 0.02,
    breaks_rho   = c(0.25, 0.5, 0.75),
    dodge_width  = 0.02
  )
  style <- utils::modifyList(defaults, style)
  
  # Preprocess for consistent ordering/labels
  data_proc <- data
  data_proc$method <- factor(as.character(data_proc$method), levels = c("RRCF","GRF"))
  if ("link_type" %in% names(data_proc))
    data_proc$link_type <- factor(as.character(data_proc$link_type), levels = c("log","identity"))
  if ("n" %in% names(data_proc))
    if (is.numeric(data_proc$n)) data_proc$n <- factor(data_proc$n, levels = sort(unique(data_proc$n)))
  if ("p_base" %in% names(data_proc)) {
    pb <- data_proc$p_base
    num_lev <- sort(unique(pb[!is.na(pb)]))
    lev <- c(if (any(is.na(pb))) "default" else character(0),
             if (length(num_lev)) format(num_lev, trim = TRUE))
    lab <- ifelse(is.na(pb), "default", format(pb, trim = TRUE))
    data_proc$p_base <- factor(lab, levels = lev)
  }
  
  y_label <- if (metric == "mape") "MAPE" else "MAE"
  #if (is.null(title))  title <- paste0(y_label, " vs Heterogeneity")
  if (is.null(title))  title <- ""
  
  # Percent axis for MAPE; raw for MAE
  y_scale <- if (metric == "mape") ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1))
  else ggplot2::scale_y_continuous()
  
  dodge <- ggplot2::position_dodge(width = style$dodge_width)
  xmin <- min(data_proc$rho, na.rm = TRUE) - style$pad_x
  xmax <- max(data_proc$rho, na.rm = TRUE) + style$pad_x
  
  p <- ggplot2::ggplot(data_proc, ggplot2::aes(x = rho, y = estimate, color = method, group = method)) +
    ggplot2::geom_line(ggplot2::aes(linetype = method), linewidth = style$line_width, position = dodge) +  # <— NEW
    ggplot2::geom_point(size = style$point_size, position = dodge) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      position = dodge, width = style$ci_cap_width, linewidth = style$ci_linewidth, alpha = style$ci_alpha
    ) +
    ggplot2::scale_x_continuous(breaks = style$breaks_rho, limits = c(xmin, xmax)) +
    y_scale +
    ggplot2::scale_color_manual(values = style$colors) +
    ggplot2::scale_linetype_manual(values = style$linetypes) +  # <— NEW
    ggplot2::labs(
      x = expression(rho ~ "(heterogeneity)"),
      y = y_label,
      color = "Method",
      linetype = "Method",
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal(base_size = style$theme_base) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  
  if (!missing(facet_by) && length(facet_by) > 0) {
    facet_formula <- stats::as.formula(paste("~", paste(facet_by, collapse = " + ")))
    p <- p + ggplot2::facet_wrap(
      facet_formula,
      nrow   = facet_nrow,                      # <— use new arg
      ncol   = facet_ncol,                      # <— use new arg
      scales = if (share_y) "fixed" else "free_y"
    )
  }
  
  outfile <- NULL
  if (is.list(save)) {
    dir.create(save$dir %||% "charts", recursive = TRUE, showWarnings = FALSE)
    ts_fmt <- save$timestamp_fmt %||% "%Y%m%d%H%M"
    stamp  <- format(Sys.time(), ts_fmt)
    base   <- if (is.null(save$name_components)) "plot" else paste(save$name_components, collapse = "_")
    outfile <- file.path(save$dir %||% "charts", paste0(base, "_", stamp, ".pdf"))
    ggplot2::ggsave(outfile, plot = p,
                    width = save$width_in  %||% 11,
                    height= save$height_in %||% 6,
                    units = "in")
    message("Saved plot to: ", outfile)
  }
  
  invisible(list(plot = p, file = outfile))
}

# ------------------------------
# Common parameters
# ------------------------------
rhos      <- c(0.25, 0.5, 0.75)
links     <- c("log","identity")
ns        <- c(100, 500, 1000, 2500)
p_bases   <- c(NA, 0.05, 0.01)  # NA => default baseline
n_default <- 2500
pbase_def <- NULL

# Small helpers to compute MAPE/MAE + CIs per (bundle, method)
get_mape_ci <- function(b) {
  d_grf  <- mape_from_bundle_ci(b, "grf")
  d_rrcf <- mape_from_bundle_ci(b, "rrcf")
  tibble::tibble(
    method    = c("GRF","RRCF"),
    estimate  = c(as.numeric(d_grf$mape),  as.numeric(d_rrcf$mape)),
    ci_lower  = c(as.numeric(d_grf$lower_se),  as.numeric(d_rrcf$lower_se)),
    ci_upper  = c(as.numeric(d_grf$upper_se), as.numeric(d_rrcf$upper_se))
  )
}

get_mae_ci <- function(b) {
  d_grf  <- mae_from_bundle_ci(b, "grf")
  d_rrcf <- mae_from_bundle_ci(b, "rrcf")
  tibble::tibble(
    method    = c("GRF","RRCF"),
    estimate  = c(as.numeric(d_grf$mae),  as.numeric(d_rrcf$mae)),
    ci_lower  = c(as.numeric(d_grf$lower_se),  as.numeric(d_rrcf$lower_se)),
    ci_upper  = c(as.numeric(d_grf$upper_se), as.numeric(d_rrcf$upper_se))
  )
}

# ======================================================================
# (1) MAPE — single chart log link
# ======================================================================
# Build only the log-link data (no identity), same rho grid
grid_log <- tidyr::expand_grid(link_type = "log", rho = rhos)

df_mape_log <- purrr::pmap_dfr(grid_log, function(link_type, rho) {
  b <- load_bundle(n = n_default, rho = rho, link_type = link_type, p_base = pbase_def)
  get_mape_ci(b) |>
    dplyr::mutate(link_type = link_type, rho = rho, metric = "mape", n = n_default, p_base = NA_real_)
}) |>
  dplyr::mutate(method = factor(method, levels = c("RRCF","GRF"))) |>
  dplyr::select(-link_type)  # drop since we're plotting a single panel

plot_with_ci(
  data     = df_mape_log,
  facet_by = NULL,                 # single panel
  metric   = "mape",
  #title    = "MAPE vs Heterogeneity (log link, CRTE scale)",
  #subtitle = "Blue: RRCF; Green: GRF. Error bars: mean ± 1 SE across replicates.",
  share_y  = FALSE,
  save     = list(dir = "charts",
                  name_components = c("log", "mape_single"),
                  timestamp_fmt = "%Y%m%d%H%M",
                  width_in = 7, height_in = 4.5)
)

# ======================================================================
# (2) MAE — single panel: identity link
# ======================================================================
# Build only the identity-link data (no log), same rho grid
grid_id <- tidyr::expand_grid(link_type = "identity", rho = rhos)

df_mae_id <- purrr::pmap_dfr(grid_id, function(link_type, rho) {
  b <- load_bundle(n = n_default, rho = rho, link_type = link_type, p_base = pbase_def)
  get_mae_ci(b) |>
    dplyr::mutate(link_type = link_type, rho = rho, metric = "mae", n = n_default, p_base = NA_real_)
}) |>
  dplyr::mutate(method = factor(method, levels = c("RRCF","GRF"))) |>
  dplyr::select(-link_type)  # drop since we're plotting a single panel

plot_with_ci(
  data     = df_mae_id,
  facet_by = NULL,                 # single panel
  metric   = "mae",
  #title    = "MAE vs Heterogeneity (identity link, CATE scale)",
  #subtitle = "Blue: RRCF; Green: GRF. Error bars: mean ± 1 SE across replicates.",
  share_y  = FALSE,
  save     = list(dir = "charts",
                  name_components = c("identity", "mae_single"),
                  timestamp_fmt = "%Y%m%d%H%M",
                  width_in = 7, height_in = 4.5)
)

# ======================================================================
# (3) MAPE by sample size (log link) — panels: n in {100, 500, 1000, 2500}
# ======================================================================

grid_n <- tidyr::expand_grid(n = ns, rho = rhos)

df_mape_n <- purrr::pmap_dfr(grid_n, function(n, rho) {
  b <- load_bundle(n = n, rho = rho, link_type = "log", p_base = pbase_def)
  get_mape_ci(b) |>
    mutate(n = n, rho = rho, link_type = "log", metric = "mape", p_base = NA_real_)
}) |>
  mutate(method = factor(method, levels = c("RRCF","GRF")),
         n = factor(n, levels = ns))

plot_with_ci(
  data         = df_mape_n,
  facet_by     = "n",
  metric       = "mape",
  #title        = "MAPE vs Heterogeneity by Sample Size (log link, CRTE scale)",
  #subtitle     = "Blue: RRCF; Green: GRF. Error bars: mean ± 1 SE across replicates.",
  share_y      = FALSE,
  facet_nrow   = 2,                 # <— 2 rows × 2 columns
  save         = list(dir = "charts",
                      name_components = c("log","mape_ci_by_n"),
                      timestamp_fmt   = "%Y%m%d%H%M",
                      width_in        = 12,      # slightly taller/wider for 2×2
                      height_in       = 8)
)

# ======================================================================
# (4) MAE by sample size (identity link) — panels: n in {100, 500, 1000, 2500}
# ======================================================================
# If not already defined earlier:
# grid_n <- tidyr::expand_grid(n = ns, rho = rhos)

# data build stays the same
df_mae_n_id <- purrr::pmap_dfr(grid_n, function(n, rho) {
  b <- load_bundle(n = n, rho = rho, link_type = "identity", p_base = pbase_def)
  get_mae_ci(b) |>
    dplyr::mutate(n = n, rho = rho, link_type = "identity", metric = "mae", p_base = NA_real_)
}) |>
  dplyr::mutate(method = factor(method, levels = c("RRCF","GRF")),
                n = factor(n, levels = ns))

# 2×2 facet layout, same filename pattern
plot_with_ci(
  data         = df_mae_n_id,
  facet_by     = "n",
  metric       = "mae",
  #title        = "MAE vs Heterogeneity by Sample Size (identity link, CATE scale)",
  #subtitle     = "Blue: RRCF; Green: GRF. Error bars: mean ± 1 SE across replicates.",
  share_y      = FALSE,
  facet_nrow   = 2,                    # <-- new: 2 rows × 2 cols
  save         = list(dir = "charts",
                      name_components = c("identity","mae_ci_by_n"),
                      timestamp_fmt   = "%Y%m%d%H%M",
                      width_in        = 12,       # adjust size for 2×2
                      height_in       = 8)
)

# ======================================================================
# (5) MAPE by baseline (log link) — panels: default, 0.05, 0.01
# ======================================================================
grid_pb <- tidyr::expand_grid(p_base = p_bases, rho = rhos)

df_mape_pb <- purrr::pmap_dfr(grid_pb, function(p_base, rho) {
  pb <- if (is.na(p_base)) NULL else p_base
  b  <- load_bundle(n = n_default, rho = rho, link_type = "log", p_base = pb)
  get_mape_ci(b) |>
    dplyr::mutate(p_base = p_base, rho = rho, link_type = "log", metric = "mape", n = n_default)
}) |>
  dplyr::mutate(
    method     = factor(method, levels = c("RRCF","GRF")),
    # create a facet label with desired ordering
    p_base_lab = factor(
      ifelse(is.na(p_base), "default", format(p_base, trim = TRUE)),
      levels = c("default", "0.05", "0.01")
    )
  )

plot_with_ci(
  data         = df_mape_pb,
  facet_by     = "p_base_lab",  # <- facet on the ordered label
  metric       = "mape",
  # title        = "MAPE vs Heterogeneity by Baseline (log link, CRTE scale)",
  # subtitle     = "Blue: RRCF; Green: GRF. Error bars: mean ± 1 SE across replicates.",
  share_y      = FALSE,
  save         = list(dir="charts",
                      name_components=c("log","mape_ci_by_baseline"),
                      timestamp_fmt="%Y%m%d%H%M", width_in=13, height_in=5.5)
)

# ======================================================================
# (6) MAE by baseline (log link) — panels: default, 0.05, 0.01
# ======================================================================
df_mae_pb <- purrr::pmap_dfr(grid_pb, function(p_base, rho) {
  pb <- if (is.na(p_base)) NULL else p_base
  b  <- load_bundle(n = n_default, rho = rho, link_type = "log", p_base = pb)
  get_mae_ci(b) |>
    dplyr::mutate(p_base = p_base, rho = rho, link_type = "log", metric = "mae", n = n_default)
}) |>
  dplyr::mutate(
    method     = factor(method, levels = c("RRCF","GRF")),
    # ordered facet labels: default -> 0.05 -> 0.01
    p_base_lab = factor(
      ifelse(is.na(p_base), "default", format(p_base, trim = TRUE)),
      levels = c("default", "0.05", "0.01")
    )
  )

plot_with_ci(
  data         = df_mae_pb,
  facet_by     = "p_base_lab",  # facet on the ordered label
  metric       = "mae",
  # title        = "MAE vs Heterogeneity by Baseline (log link, CATE scale)",
  # subtitle     = "Blue: RRCF; Green: GRF. Error bars: mean ± 1 SE across replicates.",
  share_y      = FALSE,
  save         = list(dir = "charts",
                      name_components = c("log","mae_ci_by_baseline"),
                      timestamp_fmt   = "%Y%m%d%H%M",
                      width_in        = 13, height_in = 5.5)
)


# Done — six PDFs saved under ./charts/
# --- Side-by-side MAPE (log link, n = 2500): Z = 0 vs Z = 50 -----------------

rhos        <- c(0.25, 0.5, 0.75)
n_default   <- 2500
pbase_def   <- NULL  # NULL => default baseline

build_mape_df <- function(n_extra) {
  grid <- tidyr::expand_grid(rho = rhos)
  purrr::pmap_dfr(grid, function(rho) {
    b <- load_bundle(n = n_default, rho = rho, link_type = "log",
                     p_base = pbase_def, n_extra = n_extra)
    get_mape_ci(b) |>
      dplyr::mutate(rho = rho, metric = "mape")
  }) |>
    dplyr::mutate(method = factor(method, levels = c("RRCF","GRF")))
}

df0  <- build_mape_df(0)  |> dplyr::mutate(feature_set = "Z = 0")
df50 <- build_mape_df(50) |> dplyr::mutate(feature_set = "Z = 50")

df_mape_side <- dplyr::bind_rows(df0, df50) |>
  dplyr::mutate(feature_set = factor(feature_set, levels = c("Z = 0","Z = 50")))

plot_with_ci(
  data       = df_mape_side,
  facet_by   = "feature_set",   # two panels: left = Z=0, right = Z=50
  metric     = "mape",
  share_y    = TRUE,            # same y-axis so lines are directly comparable
  facet_nrow = 1,
  facet_ncol = 2,
  save       = list(dir = "charts",
                    name_components = c("log","mape_ci_side_by_side","n2500","Z050"),
                    timestamp_fmt   = "%Y%m%d%H%M",
                    width_in        = 12,
                    height_in       = 4.5)
)