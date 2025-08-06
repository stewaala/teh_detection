source("simulate_and_fit_forests.R")

n = 2500
rho = 0.75
trials = 100
link_type = "log"
p_base = NULL

run_trials(n, rho, trials, link_type, p_base, num_trees=200)

print('Shirvaikar power test')
power <- power_from_rds(n, rho, link_type, p_base)
print(power)
power <- power_from_bundle(n, rho, link_type, p_base)
print(power)

print('MAPE')
mape <- mape_from_rds(n, rho, link_type, p_base)
print(mape)
mape <- mape_from_bundle(n, rho, link_type, p_base)
print(mape)

print('Athey power test')
power <- power_rd_from_bundle(n, rho, link_type, p_base)
print(power)