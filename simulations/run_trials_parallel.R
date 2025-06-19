#!/usr/bin/env Rscript
#
# run_trials.R
#  - Parallelized trials with logging and output files
#

# ----- 0) Load packages -----
suppressMessages({
  library(causl)
  library(data.table)
  library(grf)
  library(rrcf)
  library(ggplot2)
  library(dplyr)
  library(parallel)
})

source("copula_rct.R")

# ----- 1) Setup output directory -----
out_dir <- Sys.getenv("OUTPUT_DIR", unset = "TrialEnvs")
cat("Initializing output directory:", out_dir, "\n")
unlink(out_dir, recursive = TRUE, force = TRUE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----- 2) Define parameters -----
trials <- 11
ns    <- c(200)
rhos  <- c(0.5)

param_grid <- expand.grid(
  trial = seq_len(trials),
  n     = ns,
  rho   = rhos,
  stringsAsFactors = FALSE
)

# ----- 3) Worker function -----
run_one <- function(row) {
  trial <- row$trial
  n     <- row$n
  rho   <- row$rho
  
  start_time <- Sys.time()
  cat(sprintf("[START] trial=%03d | n=%05d | rho=%.2f at %s\n",
              trial, n, rho, start_time))
  
  # --- data generation & split ---
  data <- copula_rct(n, rho, trial)
  idx  <- sample.int(n, size = floor(0.8 * n))
  x_train <- data[idx, !c("A","Y","CRTE"), with = FALSE]
  x_test  <- data[-idx, !c("A","Y","CRTE"), with = FALSE]
  y_train <- data$Y[idx]
  y_test  <- data$Y[-idx]
  t_train <- data$A[idx]
  t_test  <- data$A[-idx]
  crte_test <- data$CRTE[-idx]
  
  # --- fit forests + predictions ---
  f_grf <- causal_forest(x_train, y_train, t_train,
                         W.hat = 0.5, num.trees = 500, seed = 1234)
  f_glm <- rr_causal_forest(x_train, y_train, t_train,
                            rct = TRUE, num.trees = 500, seed = 1234)
  p_grf <- rr_predict(f_grf, x_test)
  p_glm <- rr_predict(f_glm, x_test)
  
  # --- compute metrics ---
  mape_grf <- mean(abs((crte_test - p_grf) / crte_test))
  mape_glm <- mean(abs((crte_test - p_glm) / crte_test))
  
  # --- optional calibration tests / p-values ---
  anova_data <- data.frame(y_test, t_test, x_test)
  base_mod   <- glm(y_test ~ ., family = poisson, data = anova_data)
  mod_grf    <- glm(y_test ~ ., family = poisson,
                    data = cbind(anova_data, t_test * log(p_grf)))
  mod_glm    <- glm(y_test ~ ., family = poisson,
                    data = cbind(anova_data, t_test * log(p_glm)))
  an_grf     <- anova(base_mod, mod_grf, test = "Chisq")
  an_glm     <- anova(base_mod, mod_glm, test = "Chisq")
  pval_grf   <- 1 - pchisq(an_grf$Deviance[2], df = 1)
  pval_glm   <- 1 - pchisq(an_glm$Deviance[2], df = 1)
  
  # --- save results ---
  fname <- sprintf("T%03d_n%05d_rho%02d.Rdata",
                   trial, n, as.integer(rho * 100))
  save.image(file.path(out_dir, fname))
  
  end_time <- Sys.time()
  cat(sprintf("[DONE ] trial=%03d | duration=%s | mape=(%.4f,%.4f) | p=(%.4g,%.4g) at %s\n",
              trial,
              round(difftime(end_time, start_time, units = "secs"),1),
              mape_grf, mape_glm,
              pval_grf, pval_glm,
              end_time))
  
  return(NULL)
}

# ----- 4) Dispatch in parallel -----
n_cores <- max(detectCores(logical = FALSE) - 1, 1)
cat("Running on", n_cores, "cores\n")
mclapply(
  split(param_grid, seq_len(nrow(param_grid))),
  run_one,
  mc.cores = n_cores
)

cat("All trials finished at", Sys.time(), "\n")
