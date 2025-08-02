source('copula_rct.R')
source('alan_run_trials.R')

alan_run_trials(n=5000, rho=0.5, trials=100)


# pvals <- read.csv("data/pvals_n00200_rho50.csv", stringsAsFactors = FALSE)
# View(pvals)
