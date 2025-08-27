source("simulate_and_fit_forests.R")

# log link - different sample sizes
make_and_save_bundle(n=50, rho=0.25, link_type="log", skip_if_file_exists=TRUE)
make_and_save_bundle(n=50, rho=0.5, link_type="log", skip_if_file_exists=TRUE)
make_and_save_bundle(n=50, rho=0.75, link_type="log", skip_if_file_exists=TRUE)
make_and_save_bundle(n=100, rho=0.25, link_type="log", skip_if_file_exists=TRUE)
make_and_save_bundle(n=100, rho=0.5, link_type="log", skip_if_file_exists=TRUE)
make_and_save_bundle(n=100, rho=0.75, link_type="log", skip_if_file_exists=TRUE)
make_and_save_bundle(n=500, rho=0.25, link_type="log", skip_if_file_exists=TRUE)
make_and_save_bundle(n=500, rho=0.5, link_type="log", skip_if_file_exists=TRUE)
make_and_save_bundle(n=500, rho=0.75, link_type="log", skip_if_file_exists=TRUE)
make_and_save_bundle(n=1000, rho=0.25, link_type="log", skip_if_file_exists=TRUE)
make_and_save_bundle(n=1000, rho=0.5, link_type="log", skip_if_file_exists=TRUE)
make_and_save_bundle(n=1000, rho=0.75, link_type="log", skip_if_file_exists=TRUE)
make_and_save_bundle(n=2500, rho=0.25, link_type="log", skip_if_file_exists=TRUE)
make_and_save_bundle(n=2500, rho=0.5, link_type="log", skip_if_file_exists=TRUE)
make_and_save_bundle(n=2500, rho=0.75, link_type="log", skip_if_file_exists=TRUE)

# identity link - different sample sizes
make_and_save_bundle(n=50, rho=0.25, link_type="identity", skip_if_file_exists=TRUE)
make_and_save_bundle(n=50, rho=0.5, link_type="identity", skip_if_file_exists=TRUE)
make_and_save_bundle(n=50, rho=0.75, link_type="identity", skip_if_file_exists=TRUE)
make_and_save_bundle(n=100, rho=0.25, link_type="identity", skip_if_file_exists=TRUE)
make_and_save_bundle(n=100, rho=0.5, link_type="identity", skip_if_file_exists=TRUE)
make_and_save_bundle(n=100, rho=0.75, link_type="identity", skip_if_file_exists=TRUE)
make_and_save_bundle(n=500, rho=0.25, link_type="identity", skip_if_file_exists=TRUE)
make_and_save_bundle(n=500, rho=0.5, link_type="identity", skip_if_file_exists=TRUE)
make_and_save_bundle(n=500, rho=0.75, link_type="identity", skip_if_file_exists=TRUE)
make_and_save_bundle(n=1000, rho=0.25, link_type="identity", skip_if_file_exists=TRUE)
make_and_save_bundle(n=1000, rho=0.5, link_type="identity", skip_if_file_exists=TRUE)
make_and_save_bundle(n=1000, rho=0.75, link_type="identity", skip_if_file_exists=TRUE)
make_and_save_bundle(n=2500, rho=0.25, link_type="identity", skip_if_file_exists=TRUE)
make_and_save_bundle(n=2500, rho=0.5, link_type="identity", skip_if_file_exists=TRUE)
make_and_save_bundle(n=2500, rho=0.75, link_type="identity", skip_if_file_exists=TRUE)

#low baselines
make_and_save_bundle(n=2500, rho=0.25, link_type="log", p_base=0.05, skip_if_file_exists=TRUE)
make_and_save_bundle(n=2500, rho=0.5, link_type="log", p_base=0.05, skip_if_file_exists=TRUE)
make_and_save_bundle(n=2500, rho=0.75, link_type="log", p_base=0.05, skip_if_file_exists=TRUE)
make_and_save_bundle(n=2500, rho=0.25, link_type="log", p_base=0.01, skip_if_file_exists=TRUE)
make_and_save_bundle(n=2500, rho=0.5, link_type="log", p_base=0.01, skip_if_file_exists=TRUE)
make_and_save_bundle(n=2500, rho=0.75, link_type="log", p_base=0.01, skip_if_file_exists=TRUE)

#extra covariates
make_and_save_bundle(n=2500, rho=0.25, link_type="log", skip_if_file_exists=TRUE, n_extra=50)
make_and_save_bundle(n=2500, rho=0.5, link_type="log", skip_if_file_exists=TRUE, n_extra=50)
make_and_save_bundle(n=2500, rho=0.75, link_type="log", skip_if_file_exists=TRUE, n_extra=50)

