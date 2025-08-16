# Assumes:
#   load_bundle()
#   poisson_power_from_bundle_detail()  # returns estimate + CI

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)

source("simulate_and_fit_forests.R")

# --- Parameters ---
link_type <- "log"
n         <- 2500
rhos      <- c(0.25, 0.5, 0.75)
p_bases   <- c(NA, 0.05, 0.01)   # NA => default baseline (p_base = NULL)

# --- Collect power + CI across p_base & rho ---
grid <- tidyr::expand_grid(p_base = p_bases, rho = rhos)

get_power_row <- function(p_base, rho) {
  pb <- if (is.na(p_base)) NULL else p_base
  b  <- load_bundle(n = n, rho = rho, link_type = link_type, p_base = pb)
  
  # small_sample=TRUE for p_base 0.05 and 0.01; FALSE for default
  ss <- !is.na(p_base) && p_base %in% c(0.05, 0.01)
  
  d_grf  <- poisson_power_from_bundle_detail(b, forest_type = "grf",  small_sample = ss)
  d_rrcf <- poisson_power_from_bundle_detail(b, forest_type = "rrcf", small_sample = ss)
  
  tibble(
    p_base     = p_base,
    rho        = rho,
    method     = c("GRF", "RRCF"),
    value      = c(d_grf$estimate,  d_rrcf$estimate),
    ci_lower   = c(d_grf$ci_lower,  d_rrcf$ci_lower),
    ci_upper   = c(d_grf$ci_upper,  d_rrcf$ci_upper),
    small_sample = ss
  )
}

df_power_pb <- purrr::pmap_dfr(grid, get_power_row) |>
  mutate(
    method       = factor(method, levels = c("RRCF", "GRF")),  # blue RRCF, green GRF
    p_base_label = case_when(
      is.na(p_base) ~ "default",
      TRUE          ~ format(p_base, trim = TRUE)
    ),
    p_base_label = factor(p_base_label, levels = c("default", "0.05", "0.01"))
  )

# --- Plot: three panels (default, 0.05, 0.01), percent y, dodged/thinner CI bars ---
pd   <- position_dodge(width = 0.02)
xmin <- min(df_power_pb$rho) - 0.02
xmax <- max(df_power_pb$rho) + 0.02

p_power_by_pbase <- ggplot(df_power_pb, aes(x = rho, y = value, color = method, group = method)) +
  geom_line(linewidth = 1.1, position = pd) +
  geom_point(size = 2.6, position = pd) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    position = pd,
    width = 0.018,     # shorter caps for readability
    linewidth = 0.35,  # thinner lines
    alpha = 0.7
  ) +
  scale_x_continuous(breaks = rhos, limits = c(xmin, xmax)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_manual(values = c("RRCF" = "steelblue", "GRF" = "forestgreen")) +
  labs(
    x = expression(rho ~ "(heterogeneity)"),
    y = "Power",
    color = "Method",
    title = "Power vs Heterogeneity by Baseline (log link)",
    subtitle = "Panels: default, p_base=0.05 (small_sample), p_base=0.01 (small_sample)"
  ) +
  facet_wrap(~ p_base_label, nrow = 1, scales = "free_y") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

print(p_power_by_pbase)

# --- Save PDF to ./charts with timestamp yyyymmddhhmm ---
if (!dir.exists("charts")) dir.create("charts", recursive = TRUE)
timestamp <- format(Sys.time(), "%Y%m%d%H%M")
outfile   <- file.path("charts", paste0("log_power_by_baseline_", timestamp, ".pdf"))
ggsave(outfile, plot = p_power_by_pbase, width = 13, height = 5.5, units = "in")
cat("Saved plot to:", outfile, "\n")

# Optional: open in browser
# browseURL(normalizePath(outfile))
