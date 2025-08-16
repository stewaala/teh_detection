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
p_base    <- NULL
rhos      <- c(0.25, 0.5, 0.75)
ns        <- c(100, 500, 1000, 2500)

# --- Collect power + CI across n & rho (small_sample for n <= 500) ---
grid <- tidyr::expand_grid(n = ns, rho = rhos)

get_power_row <- function(n, rho) {
  b  <- load_bundle(n = n, rho = rho, link_type = link_type, p_base = p_base)
  ss <- n <= 500
  
  d_grf  <- poisson_power_from_bundle_detail(b, forest_type = "grf",  small_sample = ss)
  d_rrcf <- poisson_power_from_bundle_detail(b, forest_type = "rrcf", small_sample = ss)
  
  tibble(
    n        = n,
    rho      = rho,
    method   = c("GRF", "RRCF"),
    value    = c(d_grf$estimate,  d_rrcf$estimate),
    ci_lower = c(d_grf$ci_lower,  d_rrcf$ci_lower),
    ci_upper = c(d_grf$ci_upper,  d_rrcf$ci_upper),
    small_sample = ss
  )
}

df_power_n <- purrr::pmap_dfr(grid, get_power_row) |>
  mutate(
    method = factor(method, levels = c("RRCF", "GRF")),  # blue RRCF, green GRF
    n      = factor(n, levels = ns)
  )

# --- Plot (4 panels in one row), independent y-axes, percent y, dodged & lighter CIs ---
pd   <- position_dodge(width = 0.02)
xmin <- min(df_power_n$rho) - 0.02
xmax <- max(df_power_n$rho) + 0.02

p_power_by_n <- ggplot(df_power_n, aes(x = rho, y = value, color = method, group = method)) +
  geom_line(linewidth = 1.1, position = pd) +
  geom_point(size = 2.6, position = pd) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    position = pd,
    width = 0.014,     # narrower caps
    linewidth = 0.25,  # lighter/thinner CI lines
    alpha = 0.6        # lighter appearance
  ) +
  scale_x_continuous(breaks = rhos, limits = c(xmin, xmax)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = c("RRCF" = "steelblue", "GRF" = "forestgreen")) +
  labs(
    x = expression(rho ~ "(heterogeneity)"),
    y = "Power",
    color = "Method",
    title = "Power vs Heterogeneity by Sample Size (log link)",
    subtitle = "95% CIs (Wilson). small_sample=TRUE for n ∈ {100, 500}"
  ) +
  facet_wrap(~ n, nrow = 1, scales = "free_y") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

print(p_power_by_n)

# --- Save PDF to ./charts with timestamp yyyymmddhhmm ---
if (!dir.exists("charts")) dir.create("charts", recursive = TRUE)
timestamp <- format(Sys.time(), "%Y%m%d%H%M")
outfile   <- file.path("charts", paste0("log_power_by_n_", timestamp, ".pdf"))
ggsave(outfile, plot = p_power_by_n, width = 18, height = 5.5, units = "in")
cat("Saved plot to:", outfile, "\n")

# Optional: open in browser
# browseURL(normalizePath(outfile))
