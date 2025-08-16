# Assumes:
#   load_bundle()
#   poisson_power_from_bundle_detail()  # returns estimate + CI

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)

source("simulate_and_fit_forests.R")

# --- Parameters for this set ---
n <- 2500
p_base <- NULL
rhos <- c(0.25, 0.5, 0.75)
link_types <- c("log", "identity")

# --- Collect power + CI for both methods across rho & link types ---
grid <- tidyr::expand_grid(link_type = link_types, rho = rhos)

get_power_row <- function(link_type, rho) {
  b <- load_bundle(n = n, rho = rho, link_type = link_type, p_base = p_base)
  
  d_grf  <- poisson_power_from_bundle_detail(b, forest_type = "grf")
  d_rrcf <- poisson_power_from_bundle_detail(b, forest_type = "rrcf")
  
  tibble(
    link_type = link_type,
    rho       = rho,
    method    = c("GRF", "RRCF"),
    value     = c(d_grf$estimate,  d_rrcf$estimate),
    ci_lower  = c(d_grf$ci_lower,  d_rrcf$ci_lower),
    ci_upper  = c(d_grf$ci_upper,  d_rrcf$ci_upper),
    ci_note   = d_grf$ci_method
  )
}

df_power <- purrr::pmap_dfr(grid, get_power_row) |>
  dplyr::mutate(
    link_type = factor(link_type, levels = c("log", "identity")),
    method    = factor(method, levels = c("RRCF", "GRF")) # blue RRCF, green GRF
  )

# --- Plot: two panels in one row, independent y-axes, thin CI lines ---
pd   <- position_dodge(width = 0.02)   # separate CI caps
xmin <- min(df_power$rho) - 0.02
xmax <- max(df_power$rho) + 0.02

p_power <- ggplot(df_power, aes(x = rho, y = value, color = method, group = method)) +
  geom_line(linewidth = 1.1, position = pd) +
  geom_point(size = 2.6, position = pd) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    position = pd,
    width = 0.014,     # narrow caps
    linewidth = 0.25,  # thin CI lines
    alpha = 0.6
  ) +
  scale_x_continuous(breaks = rhos, limits = c(xmin, xmax)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +  # power as %
  scale_color_manual(values = c("RRCF" = "steelblue", "GRF" = "forestgreen")) +
  labs(
    x = expression(rho ~ "(heterogeneity)"),
    y = "Power",
    color = "Method",
    title = "Power vs Heterogeneity (n = 2500, default baseline)",
    subtitle = "Poisson omnibus test; 95% CIs (Wilson); slight horizontal dodge to avoid overlap"
  ) +
  facet_wrap(~ link_type, nrow = 1, scales = "free_y") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

print(p_power)

# --- Save PDF to ./charts with timestamp yyyymmddhhmm ---
if (!dir.exists("charts")) dir.create("charts", recursive = TRUE)
timestamp   <- format(Sys.time(), "%Y%m%d%H%M")
link_labels <- paste(sort(unique(as.character(df_power$link_type))), collapse = "_")
outfile     <- file.path("charts", paste0(link_labels, "_power_", timestamp, ".pdf"))
ggsave(outfile, plot = p_power, width = 11, height = 6, units = "in")
cat("Saved plot to:", outfile, "\n")

# Optional: open in browser
# browseURL(normalizePath(outfile))
