#!/usr/bin/env Rscript
#
# run_trials_parallel.R
#  - Parallelized trials with link_type & baseline args
#

# ----- 0) Load packages -----
suppressMessages({
  library(data.table)
  library(grf)
  library(rrcf)
  library(parallel)
})

# ----- 1) Parse command-line arguments -----
args <- commandArgs(trailingOnly = TRUE)
# Defaults
link_type <- "log"
#p_base    <- 0.10
p_base    <- 0.16253

# Helper to extract named args: --link_type and --baseline
if(length(args) > 0) {
  for(i in seq(1, length(args), by = 2)) {
    name  <- args[i]
    value <- args[i + 1]
    if(name == "--link_type") {
      link_type <- value
    } else if(name == "--baseline") {
      p_base <- as.numeric(value)
    }
  }
}

cat("Running with link_type =", link_type, "and p_base =", p_base, "\n")

# ----- 2) Source the simulator -----
source("copula_rct_sim.R")  # must define copula_rct_sim()

# ----- 3) Setup output directory -----
out_dir <- Sys.getenv("OUTPUT_DIR", unset = "TrialEnvs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----- 4) Define parameter grid -----
trials     <- 100
ns         <- c(2500, 5000, 7500, 10000)
rhos       <- c(0.25, 0.50, 0.75)
param_grid <- expand.grid(
  trial     = seq_len(trials),
  n         = ns,
  rho       = rhos,
  stringsAsFactors = FALSE
)

# ----- 5) Worker function -----
run_one <- function(row) {
  trial <- row$trial
  n     <- row$n
  rho   <- row$rho
  
  # filename includes all key parameters
  fname <- sprintf(
    "T%03d_n%05d_rho%02d_link-%s_base-%.5f.Rdata",
    trial, n, as.integer(rho * 100),
    link_type, p_base
  )
  outfile <- file.path(out_dir, fname)
  if (file.exists(outfile)) {
    message("[SKIP] ", fname)
    return(NULL)
  }
  
  # --- Data generation using new simulator ---
  dat <- copula_rct_sim(
    n         = n,
    rho       = rho,
    p_base    = p_base,
    link_type = link_type,
    seed      = trial,
    clip_C4   = TRUE
  )
  
  # Extract true effects depending on link_type
  if (link_type == "log") {
    true_eff <- dat$CRTE
  } else {
    true_eff <- dat$CATE
  }
  
  # --- Train/test split ---
  idx       <- sample.int(n, size = floor(0.8 * n))
  x_train   <- dat[idx,   !c("Y","CRTE","CATE"), with = FALSE]
  x_test    <- dat[-idx,  !c("Y","CRTE","CATE"), with = FALSE]
  y_train   <- dat$Y[idx]
  y_test    <- dat$Y[-idx]
  t_train   <- dat$A[idx]
  t_test    <- dat$A[-idx]
  eff_test  <- true_eff[-idx]
  
  # --- Fit forests + predictions ---
  f_grf  <- causal_forest(x_train, y_train, t_train,
                          W.hat = 0.5, num.trees = 500, seed = 1234)
  f_rrcf <- rr_causal_forest(x_train, y_train, t_train,
                             rct = TRUE, num.trees = 500, seed = 1234)
  
  pred_grf  <- rr_predict(f_grf,  x_test)
  pred_rrcf <- rr_predict(f_rrcf, x_test)
  
  # --- Compute MAPE ---
  mape_grf  <- mean(abs((eff_test - pred_grf ) / eff_test))
  mape_rrcf <- mean(abs((eff_test - pred_rrcf) / eff_test))
  
  # --- Heterogeneity power test via Poisson deviance ---
  anova_df  <- data.frame(y_test, x_test)
  base_glm  <- glm(y_test ~ ., family = poisson, data = anova_df)
  
  # safe clamp to avoid log(0)
  eps <- .Machine$double.xmin
  p1_grf  <- pmax(pred_grf,  eps)
  p1_rrcf <- pmax(pred_rrcf, eps)
  
  cov_grf  <- t_test * log(p1_grf)
  cov_rrcf <- t_test * log(p1_rrcf)
  
  glm_grf  <- glm(y_test ~ . + cov_grf,
                  family = poisson,
                  data   = anova_df)
  glm_rrcf <- glm(y_test ~ . + cov_rrcf,
                  family = poisson,
                  data   = anova_df)
  
  an_grf   <- anova(base_glm, glm_grf)
  an_rrcf  <- anova(base_glm, glm_rrcf)
  pval_grf <- 1 - pchisq(an_grf$Deviance[2], df = 1)
  pval_rrcf<- 1 - pchisq(an_rrcf$Deviance[2], df = 1)
  
  # --- Save results ---
  save(
    mape_grf, mape_rrcf,
    pval_grf, pval_rrcf,
    file      = outfile,
    compress  = "xz"
  )
  
  message(sprintf("[DONE] %s | MAPE=(%.4f,%.4f) | p=(%.4g,%.4g)",
                  fname, mape_grf, mape_rrcf, pval_grf, pval_rrcf))
  return(NULL)
}

# ----- 6) Dispatch in parallel -----
n_cores <- max(detectCores(logical = FALSE) - 1, 1)
mclapply(split(param_grid, seq_len(nrow(param_grid))),
         run_one,
         mc.cores = n_cores)

message("All trials complete.")
