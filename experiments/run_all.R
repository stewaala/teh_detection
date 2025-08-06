# ------------------------------------------------------------------
# 0.  create bundle files
# ------------------------------------------------------------------
n_values = c(2500, 5000, 7500, 10000)
rho_values = c(0.25, 0.50, 0.75)
run_trials_grid(n_values, rho_values, link_type="log", p_base= NULL, trials=100, num_trees=200)
run_trials_grid(n_values, rho_values, link_type="log", p_base= 0.01, trials=100, num_trees=200)
run_trials_grid(n_values, rho_values, link_type="identity", p_base=NULL, trials=100, num_trees=200)

# ------------------------------------------------------------------
# 1.  POWER (Experiment 1)
# ------------------------------------------------------------------
create_power_plot(link_type="log", p_base=NULL, test_type="RR", n_values=n_values, rho_values=rho_values)
create_power_plot(link_type="log",  p_base=0.01, test_type="RR", n_values=n_values, rho_values=rho_values)
create_power_plot(link_type="log", p_base=NULL, test_type = "RD", n_values=n_values, rho_values=rho_values)
create_power_plot(link_type="identity", p_base=NULL, test_type="RD", n_values=n_values, rho_values=rho_values)

# ------------------------------------------------------------------
# 2.  MAPE (Experiment 2)
# ------------------------------------------------------------------
create_mape_plot(link_type="log", p_base=NULL, n_values=n_values, rho_values=rho_values)
create_mape_plot(link_type="log", p_base=0.01, n_values=n_values, rho_values=rho_values)
create_mape_plot(link_type="identity", p_base=NULL, n_values=n_values, rho_values=rho_values)

