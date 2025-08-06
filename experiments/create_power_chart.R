library(ggplot2)

# ------------------------------------------------------------
# Dummy power calculators (replace with real ones later)
# ------------------------------------------------------------
power_rr_dummy <- function(n, rho, ...) {
  val_glm <- pmin(1, 0.10 + 0.30*rho + 0.00009 * n)   # RRCF
  val_grf <- pmin(1, 0.05 + 0.20*rho + 0.00006 * n)   # GRF
  c(power_grf = val_grf, power_glm = val_glm)
}

power_rd_dummy <- function(n, rho, ...) {
  val_grf <- pmin(1, 0.08 + 0.28*rho + 0.00008 * n)
  val_glm <- pmin(1, 0.04 + 0.18*rho + 0.00005 * n)
  c(power_grf = val_grf, power_glm = val_glm)
}

# ------------------------------------------------------------
# Figure‑1 style power plot
# ------------------------------------------------------------
create_power_plot <- function(link_type  = c("log", "identity"),
                              p_base     = NULL,
                              test_type  = c("RR", "RD"),   # RR = Shirvaikar, RD = Athey
                              n_values   = c(2500, 5000, 7500, 10000),
                              rho_values = c(0.25, 0.50, 0.75)) {
  
  link_type <- match.arg(link_type)
  test_type <- match.arg(test_type)
  get_power <- if (test_type == "RR") power_rr_dummy else power_rd_dummy
  
  df <- data.frame()
  for (n in n_values) {
    for (rho in rho_values) {
      pow <- get_power(n = n, rho = rho, link_type = link_type, p_base = p_base)
      df <- rbind(df,
                  data.frame(n = n, rho = factor(rho), Method = "GRF",  Power = 100*pow["power_grf"]),
                  data.frame(n = n, rho = factor(rho), Method = "RRCF", Power = 100*pow["power_glm"]))
    }
  }
  
  plt <- ggplot(df, aes(x = n, y = Power,
                        colour = rho,
                        linetype = Method,
                        group = interaction(Method, rho))) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_colour_manual(values = c("0.25" = "#D55E00",
                                   "0.5"  = "#009E73",
                                   "0.75" = "#0072B2")) +
    scale_linetype_manual(values = c("GRF" = "dotted",
                                     "RRCF" = "solid")) +
    labs(x = "Sample Size (n)",
         y = "Power (%)",
         colour = "Heterogeneity (rho)",
         linetype = "Method") +
    ggtitle(sprintf("Power curves (%s DGP, %s test)",
                    link_type,
                    ifelse(test_type == "RR", "RR‑scale", "RD‑scale"))) +
    theme_minimal(base_size = 12)
  
  print(plt)
  
  caption <- sprintf(
    "Power (proportion of trials where %s omnibus test detected heterogeneity) \
across 100 simulated RCTs. Curves vary by sample size n and heterogeneity level rho; \
dotted lines denote GRF, solid lines RRCF. DGP link: %s%s.",
    ifelse(test_type == "RR", "relative‑risk", "risk‑difference"),
    link_type,
    ifelse(is.null(p_base), "", paste0(", baseline risk = ", p_base))
  )
  cat("\nFigure caption:\n", caption, "\n")
  
  invisible(plt)
}

# --- example call ---------------------------------------------------
create_power_plot(link_type = "log",
                  p_base    = 0.13,
                  test_type = "RR")
