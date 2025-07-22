library(data.table)
library(ggplot2)

# Simulation settings
n_sim   <- 1000
p_base  <- 0.1
rho_vals <- c(0.25, 0.50, 0.75)

# Simulate for each rho and collect p0 and CRTE
dt_list <- lapply(rho_vals, function(rho) {
  dat <- copula_rct_sim(
    n         = n_sim,
    rho       = rho,
    p_base    = p_base,
    link_type = "log",
    seed      = 123,
    clip_C4   = TRUE
  )
  data.table(
    p0   = dat$p0,
    CRTE = dat$CRTE,
    rho  = factor(rho, levels = rho_vals)
  )
})
df <- rbindlist(dt_list)

# Scatter plot: p0 vs CRTE, colored by rho
ggplot(df, aes(x = p0, y = CRTE, color = rho)) +
  geom_point(alpha = 0.4, size = 1.2) +
  scale_color_manual(values = c("red","blue","green"),
                     name = expression(rho)) +
  labs(
    title = "CRTE vs Untreated Risk (log link)",
    x     = expression(p[0]),
    y     = "CRTE"
  ) +
  theme_minimal()
