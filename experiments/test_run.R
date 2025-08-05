source("simulate_and_fit_forests.R")

n = 2500
rho = 0.75
trials = 3
link_type = "identity"
p_base = NULL

run_trials(n, rho, trials, link_type, p_base, num_trees=200)

# power <- power_from_rds(n, rho, link_type, p_base)
# print(power)
# power <- power_from_bundle(n, rho, link_type, p_base)
# print(power)

# mape <- mape_from_rds(n, rho, link_type, p_base)
# print(mape)
# mape <- mape_from_bundle(n, rho, link_type, p_base)
# print(mape)

