# ================================================
# Reproduce Figure: Average MAPE vs. Sample Size
# ================================================

# --- 0) PARAMETERS (modify as needed) ---
out_dir <- "TrialEnvs"         # folder where your .Rdata files live
output_plot <- "mape_plot.png" # name of the saved chart
# (The code will auto‐discover whatever n and rho you actually ran.)

# --- 1) PACKAGES ---
library(ggplot2)
library(dplyr)
library(tidyr)

# --- 2) FIND ALL RESULT FILES ---
files <- list.files(
  path       = out_dir,
  pattern    = "\\.Rdata$",
  full.names = TRUE
)

# --- 3) READ & EXTRACT METRICS INTO A DATA FRAME ---
res_list <- lapply(files, function(fpath) {
  # ---- 3a) Parse trial, n, rho from the filename ----
  # expects names like "T001_n00100_rho25.Rdata"
  fname   <- basename(fpath)
  parts   <- sub("^T(\\d+)_n(\\d+)_rho(\\d+)\\.Rdata$",
                 "\\1;\\2;\\3",
                 fname)
  parts   <- strsplit(parts, ";")[[1]]
  trial   <- as.integer(parts[1])
  n_val   <- as.integer(parts[2])
  rho_val <- as.integer(parts[3]) / 100
  
  # ---- 3b) Load into a fresh environment ----
  env <- new.env()
  load(fpath, envir = env)
  
  # ---- 3c) Grab the two MAPE metrics ----
  data.frame(
    trial = trial,
    n     = n_val,
    rho   = rho_val,
    GRF   = env$mape_grf,  # dotted line in your figure
    RRCF  = env$mape_glm   # solid line in your figure
  )
})

res_df <- bind_rows(res_list)

# --- 4) RESHAPE AND AVERAGE ACROSS TRIALS ---
plot_df <- res_df %>%
  pivot_longer(
    cols      = c(GRF, RRCF),
    names_to  = "Method",
    values_to = "MAPE"
  ) %>%
  group_by(n, rho, Method) %>%
  summarize(Avg_MAPE = mean(MAPE), .groups = "drop")

# --- 5) MAKE THE CHART ---
p <- ggplot(plot_df, aes(x = n, y = Avg_MAPE,
                         color     = factor(rho),
                         linetype  = Method)) +
  #geom_line(size = 1) +
  geom_line(linewidth = 1) +

  scale_color_manual(
    name   = "Heterogeneity (ρ)",
    # note: levels of factor(rho) are "0.25", "0.5", "0.75"
      values = c("0.25" = "red",
                                 "0.5"  = "green",
                                 "0.75" = "blue"),
      breaks = c("0.25","0.5","0.75")
    ) +
  
  scale_linetype_manual(
    name   = "Method",
    values = c("GRF"  = "dotted",
               "RRCF" = "solid")
  ) +
  scale_x_continuous(breaks = unique(plot_df$n)) +
  labs(x = "Sample Size (n)",
       y = "Average MAPE") +
  theme_minimal()

# --- 6) DISPLAY & SAVE ---
print(p)
ggsave(output_plot, plot = p, width = 6, height = 4, dpi = 300)

