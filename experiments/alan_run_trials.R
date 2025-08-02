alan_run_trials <- function(n, rho, trials=100) {
  
  ## --------------------  package setup  -------------------- ##
  library(causl)
  library(data.table)
  library(grf)
  library(rrcf)
  library(parallel)
  
  out_dir <- "data"
  num_trees <- 200
  
  ## --------------------  helper : one trial  ---------------- ##
  # idx runs 1:trials; we pass idx as RNG seed to copula_rct()
  one_trial <- function(idx) {
    
    ## 1.  draw data  -------------------------------------------------
    dat <- copula_rct(n, rho, idx)
    
    split <- sample(seq_len(n), size = 0.8 * n)
    x          <- dat[, !c("A", "Y", "CRTE"), with = FALSE]
    x.train    <- x[ split, ]
    x.test     <- x[-split, ]
    y.train    <- dat$Y[ split]
    y.test     <- dat$Y[-split]
    t.train    <- dat$A[ split]
    t.test     <- dat$A[-split]
    crte.test  <- dat$CRTE[-split]
    
    ## 2.  fit forests  -----------------------------------------------
    forest.grf <- causal_forest(x.train, y.train, t.train, W.hat = 0.5, num.trees = num_trees, seed = 1234, num.threads  = 1)
    forest.glm <- rr_causal_forest(x.train, y.train, t.train, rct = TRUE, num.trees = num_trees, seed = 1234, num.threads  = 1)
    
    ## 3.  predictions & MAPE  ----------------------------------------
    pred.grf <- rr_predict(forest.grf, x.test)
    pred.glm <- rr_predict(forest.glm, x.test)
    
    mape.grf <- mean(abs((crte.test - pred.grf) / crte.test))
    mape.glm <- mean(abs((crte.test - pred.glm) / crte.test))
    
    ## 4.  calibration p‑values (LR test) ------------------------------
    anova.data <- data.frame(cbind(y.test, t.test, x.test))
    base       <- glm(y.test ~ ., family = poisson, data = anova.data)
    
    add.grf <- glm(y.test ~ ., family = poisson, data = cbind(anova.data, t.test * log(pred.grf)))
    add.glm <- glm(y.test ~ ., family = poisson, data = cbind(anova.data, t.test * log(pred.glm)))
    
    p.grf <- 1 - pchisq(anova(base, add.grf)$Deviance[2], 1)
    p.glm <- 1 - pchisq(anova(base, add.glm)$Deviance[2], 1)
    
    ## 5.  return summary rows  ---------------------------------------
    list(pvals = c(n, rho, p.grf,  p.glm), mapes = c(n, rho, mape.grf, mape.glm))
  }
  ## --------------------  parallel execution  ----------------------- ##
  res_lst <- mclapply(seq_len(trials), FUN=one_trial, mc.cores=parallel::detectCores()-1, mc.preschedule = FALSE)
  ## --------------------  bind results  ----------------------------- ##
  pvals <- do.call(rbind, lapply(res_lst, `[[`, "pvals"))
  mapes <- do.call(rbind, lapply(res_lst, `[[`, "mapes"))
  colnames(pvals) <- c("n", "rho", "p_grf", "p_glm")
  colnames(mapes) <- c("n", "rho", "mape_grf", "mape_glm")
  
  ## --------------------  write CSVs  ------------------------------- ##
  fname_base <- sprintf("n%05d_rho%02d", n, round(rho * 100))
  write.csv(pvals, file = file.path(out_dir, paste0("pvals_", fname_base, ".csv")), row.names = FALSE)
  write.csv(mapes, file = file.path(out_dir, paste0("mapes_", fname_base, ".csv")), row.names = FALSE)
  
  invisible(list(pvals = pvals, mapes = mapes))
}
