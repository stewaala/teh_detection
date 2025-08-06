library(ggplot2)

# ------------------------------------------------------------
# Dummy MAPE calculator (swap in mape_from_bundle() later)
# ------------------------------------------------------------
mape_dummy <- function(n, rho, ...) {
  ## MAPE declines with n, rises with rho, and GRF is ~10–15 % worse
  base_rrcf <- 0.50 - 0.00003 * n + 0.10 * (rho - 0.25)
  base_grf  <- base_rrcf * 1.12                 # ~12 % higher error
  c(mean_mape_grf = base_grf,
    mean_mape_glm = base_rrcf)
}

# ------------------------------------------------------------
# Figure‑2 style MAPE plot
# ------------------------------------------------------------
create_mape_plot <- function(link_type  = c("log", "identity"),
                             p_base     = NULL,
                             n_values   = c(2500, 5000, 7500, 10000),
                             rho_values = c(0.25, 0.50, 0.75)) {
  
  link_type <- match.arg(link_type)
  
  df <- data.frame()
  for (n in n_values) {
    for (rho in rho_values) {
      mpe <- mape_dummy(n = n, rho = rho, link_type = link_type, p_base = p_base)
      df  <- rbind(df,
                   data.frame(n = n, rho = factor(rho), Method = "GRF",  MAPE = mpe["mean_mape_grf"]),
                   data.frame(n = n, rho = factor(rho), Method = "RRCF", MAPE = mpe["mean_mape_glm"]))
    }
  }
  
  plt <- ggplot(df, aes(x = n, y = MAPE,
                        colour   = rho,
                        linetype = Method,
                        group    = interaction(Method, rho))) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_colour_manual(values = c("0.25" = "#D55E00",
                                   "0.5"  = "#009E73",
                                   "0.75" = "#0072B2")) +
    scale_linetype_manual(values = c("GRF" = "dotted",
                                     "RRCF" = "solid")) +
    labs(x = "Sample Size (n)",
         y = "Average MAPE",
         colour = "Heterogeneity (rho)",
         linetype = "Method") +
    ggtitle(sprintf("Average MAPE (%s DGP)", link_type)) +
    theme_minimal(base_size = 12)
  
  print(plt)
  
  caption <- sprintf(
    "Average mean absolute percentage error (MAPE) on predicted CRTE across \
100 simulated RCTs.  Curves vary by sample size n and heterogeneity level rho; \
dotted lines denote GRF, solid lines RRCF.  DGP link: %s%s.",
    link_type,
    ifelse(is.null(p_base), "", paste0(", baseline risk = ", p_base))
  )
  cat("\nFigure caption:\n", caption, "\n")
  
  invisible(plt)
}

## --- example usage ----------------------------------------------
create_mape_plot(link_type = "log", p_base = NULL)
