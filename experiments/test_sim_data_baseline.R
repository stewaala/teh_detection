# Load required libraries
library(data.table)
library(ggplot2)
library(patchwork)

# Simulation settings
n_sim   <- 1000           # default sample size
rho_sim <- 0.5            # heterogeneity scale
p_bases <- c(0.1, 0.5, 0.9)

# Helper to simulate and extract p0 for a given link type
get_p0_dt <- function(link_type) {
  lapply(p_bases, function(pb) {
    dat <- copula_rct_sim(
      n         = n_sim,
      rho       = rho_sim,
      p_base    = pb,
      link_type = link_type,
      seed      = 123,
      clip_C4   = TRUE
    )
    data.table(
      p0     = dat$p0,
      p_base = factor(pb, levels = p_bases)
    )
  }) |> rbindlist()
}

# Prepare data
df_log <- get_p0_dt("log")
df_id  <- get_p0_dt("identity")

# Compute shared y-axis limits
y_limits <- range(c(df_log$p0, df_id$p0))

# Build log-link panel with internal title
p_log <- ggplot(df_log, aes(x = p_base, y = p0)) +
  geom_boxplot(outlier.size = 0.6) +
  ggtitle("Log link") +
  labs(x = "Target untreated risk") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    plot.title   = element_text(hjust = 0.5)
  ) +
  coord_cartesian(ylim = y_limits)

# Build identity-link panel with internal title
p_id <- ggplot(df_id, aes(x = p_base, y = p0)) +
  geom_boxplot(outlier.size = 0.6) +
  ggtitle("Identity link") +
  labs(x = "Target untreated risk") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    plot.title   = element_text(hjust = 0.5)
  ) +
  coord_cartesian(ylim = y_limits)

# Combine side by side
combined_plot <- p_log + p_id + plot_layout(ncol = 2)

# Display the plot
print(combined_plot)

# Save to PDF with timestamped filename
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
chart_path = '/home/rstudio/teh_detection/experiments/charts/'
filename  <- paste0(chart_path, "p0_boxplots_", timestamp, ".pdf")

ggsave(
  filename,
  plot   = combined_plot,
  width  = 10,        # in inches
  height = 5,         # in inches
  device = cairo_pdf  # vector output for LaTeX
)

message("Saved plot to: ", filename)
