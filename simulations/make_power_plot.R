# ================================================
# Reproduce Figure: Power vs. Sample Size
# ================================================

# --- 0) PARAMETERS (modify as needed) ---
out_dir     <- "TrialEnvs"         # where your .Rdata files live
output_plot <- "power_plot.png"    # where the chart will be saved
alpha       <- 0.05                # significance threshold for “power”
# (This will pick up whatever n and rho you actually ran.)

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

# --- 3) READ & EXTRACT P-VALUES INTO A DATA FRAME ---
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
  
  # ---- 3c) Grab the two ANOVA p-values ----
  data.frame(
    trial    = trial,
    n        = n_val,
    rho      = rho_val,
    GRF_p    = env$pval_grf,  # dotted baseline
    RRCF_p   = env$pval_glm   # solid RR‐forest
  )
})

res_df <- bind_rows(res_list)

# --- 4) RESHAPE & COMPUTE POWER (%) ACROSS TRIALS ---
power_df <- res_df %>%
  pivot_longer(
    cols      = c(GRF_p, RRCF_p),
    names_to  = "Method",
    values_to = "p_value"
  ) %>%
  mutate(
    Method = ifelse(Method == "GRF_p",  "GRF", 
                    ifelse(Method == "RRCF_p", "RRCF", NA))
  ) %>%
  group_by(n, rho, Method) %>%
  summarize(
    Power = mean(p_value < alpha) * 100,
    .groups = "drop"
  )

# --- 5) PLOT POWER vs. SAMPLE SIZE ---
p <- ggplot(power_df, 
            aes(x = n, y = Power,
                color     = factor(rho),
                linetype  = Method)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    name   = "Heterogeneity (ρ)",
    # factor(rho) levels: "0.25","0.5","0.75"
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
  scale_x_continuous(breaks = unique(power_df$n)) +
  labs(x = "Sample Size (n)",
       y = "Power (%)") +
  theme_minimal()

# --- 6) DISPLAY & SAVE ---
print(p)
ggsave(output_plot, plot = p, width = 6, height = 4, dpi = 300)
