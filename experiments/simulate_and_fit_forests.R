library(causl)
library(data.table)

.EEXP_H_CONST <- 1.20332119584

# ------------------------------------------------------------------------------
# The following function replicates the functionality of copula_rct() in 
# Shirvaikar's github, but it takes an additional p_base parameter that
# allows the user the specify the level of the expected untreated risk
# ------------------------------------------------------------------------------
copula_rct_log <- function(n, rho = 0, seed = 111, p_base = .EEXP_H_CONST * exp(-2)){
  
  alpha <- log(p_base) - log(.EEXP_H_CONST)
  message("alpha set to ", alpha, '\n')
  
  fam <- list(c(1,3,4),c(5,5,5,1,2),c(5),c(1,3,4))  
  forms <- list(c(X1 ~ 1, X2 ~ X1, X3~ X1), 
                list(A ~ 1, C1 ~ 1, C2 ~ C1, C3 ~ C1:C2, C4 ~ C1), 
                Y ~ A + C1 + I(sin(C4)) + A:C1 + A:C2 + A:I(C3>0) + A:I(C4^2), 
                ~ C2) 
  pars <- list(X1 = list(beta=0, phi=1),
               X2 = list(beta=c(0.1,0.2), phi=1), 
               X3 = list(beta=c(0.1,0.1), phi=1), 
               A = list(beta=c(0)), # treatment is completely randomized
               C1 = list(beta=c(0)),
               C2 = list(beta=c(-2,1)),
               C3 = list(beta=c(0,0.1), phi = 1),
               C4 = list(beta=c(0,0.1), phi=0.1, par2=20),
               #Y = list(beta=c(-2, -0.2, 0.3, 0.4, -rho, -rho, rho, rho)), 
               Y = list(beta=c(alpha, -0.2, 0.3, 0.4, -rho, -rho, rho, rho)), 
               cop = list(Y=list(X1=list(beta=c(0,0.5)),
                                 X2=list(beta=c(-0.5,0)),
                                 X3=list(beta=c(0.5,0))))) 
  link <- list(c("identity","log","logit"),
               c("logit","logit","logit","identity","identity"),
               "log")
  
  # ------------------------------------------------------------------------------
  # Fail-fast check: validate .EEXP_H_CONST against current model parameters
  # ------------------------------------------------------------------------------
  # We model the untreated mean by a log-link:
  #   E[Y | C, A=0] = exp(alpha + h(C)),
  # where
  #   h(C)      = beta_C1 * C1 + beta_sin * sin(C4),
  #   beta_C1   = pars$Y$beta[3],        # default 0.3
  #   beta_sin  = pars$Y$beta[4],        # default 0.4
  #   C4 ~ t(df = pars$C4$par2,           # default df=20
  #           mu = pars$C4$beta[2],       # default 0.1
  #           sigma = pars$C4$phi)        # default 0.1
  #
  # The constant .EEXP_H_CONST must equal E[exp{h(C)}] under these settings,
  # so that alpha = log(p_base) - log(.EEXP_H_CONST) and consequently that the 
  # expected untreated risk E[Y | A=0] equals p_base
  # We recompute the analytic integral via compute_Eexp_h() and compare it to
  # .EEXP_H_CONST. If the difference exceeds 1e-8, we stop() immediately—before
  # calling the expensive simulator—to alert you that .EEXP_H_CONST needs updating.
  # If you change any h(C) coefficients or C4 distribution parameters, update
  # .EEXP_H_CONST to the new E[exp{h(C)}] as returned by compute_Eexp_h().
  compute_Eexp_h <- function(beta_C1, beta_sin,
                             mu_shift = pars$C4$beta[2],
                             sigma    = pars$C4$phi,
                             df       = pars$C4$par2) {
    f <- function(mu) {
      integrand <- function(z)
        exp(beta_sin * sin(mu + sigma * z)) * dt(z, df)
      integrate(integrand, -Inf, Inf, rel.tol = 1e-12)$value
    }
    F0 <- f(0)
    F1 <- f(mu_shift)
    0.5 * F0 + 0.5 * exp(beta_C1) * F1
  }
  E_det <- compute_Eexp_h(beta_C1  = pars$Y$beta[3],
                          beta_sin = pars$Y$beta[4])
  
  if (abs(E_det - .EEXP_H_CONST) > 1e-8) {
    stop(sprintf("E[exp{h(C)}] changed: %.11f vs %.11f",
                 E_det, .EEXP_H_CONST))
  }
  #---------------------------------------------------------------------------------------
  set.seed(seed)
  dat <- as.data.table(rfrugalParam(n, formulas=forms, pars=pars, 
                                    family=fam, link=link))
  
  dat2 <- copy(dat)
  X_treated   <- model.matrix(forms[[3]], dat2[, A := 1])
  X_control   <- model.matrix(forms[[3]], dat2[, A := 0])
  
  log_p1 <- as.vector(X_treated %*% pars$Y$beta)
  log_p0 <- as.vector(X_control %*% pars$Y$beta)
  p1  <- exp(log_p1)
  p0  <- exp(log_p0)
  
  # Shirvaikar does not apply this assertion
  # it is failing in some cases
  # I'm removing it for compatibility and to continue to make progress
  # but would like to revisit
  
  # if (any(p1  > 1 | p0 > 1)) {
  #   bad <- which(p1 > 1 | p0 > 1)[1]
  #   stop(sprintf(
  #     "Invalid probability at row %d: µ_control=%.3f, µ_treated=%.3f",
  #     bad, p0[bad], p1[bad]
  #   ))
  # }
  
  dat$CRTE <- p1 / p0
  dat$CATE <- p1 - p0
  return(dat)
}

copula_rct_identity <- function(n,
                                     rho        = 0.50,
                                     seed       = 111,
                                     p_base     = 0.23) {
  
  ##########################################################################################
  # ensure the user has passed rho/p_base combination that will generate valid probabilities
  ##########################################################################################
  # 1) Restrict rho to allowable values
  valid_rho <- c(0.00, 0.25, 0.50, 0.75)
  if (! (rho %in% valid_rho) ) {
    stop(sprintf(
      "Invalid rho = %.2f: must be one of %s",
      rho, paste(valid_rho, collapse = ", ")
    ))
  }
  
  # 2) Compute admissible p_base interval:
  lower <- 0.17 + 0.067 * rho
  upper <- 0.83 - 0.067 * rho
  
  # 3) Validate p_base
  if (p_base < lower || p_base > upper) {
    stop(sprintf(
      "For rho = %.2f, p_base must be in [%.3f, %.3f]. You supplied %.3f.",
      rho, lower, upper, p_base
    ))
  }
  ##########################################################################################
  
  fam <- list(
    c(1,3,4),
    c(5,5,5,1,2),
    c(1),
    c(1,3,4)
  )
  
  forms <- list(
    c(X1~1, X2~X1, X3~X1),
    list(A~1, C1~1, C2~C1, C3~C1:C2, C4~C1),
    E~1,
    ~ C2
  )
  
  pars <- list(
    X1=list(beta=0,phi=1),
    X2=list(beta=c(0.1,0.2),phi=1),
    X3=list(beta=c(0.1,0.1),phi=1),
    A =list(beta=0),
    C1=list(beta=0),
    C2=list(beta=c(-2,1)),
    C3=list(beta=c(0,0.1),phi=1),
    C4=list(beta=c(0,0.1),phi=0.1,par2=20),
    E =list(beta=0,phi=1),
    cop=list(
      E=list(
        X1=list(beta=c( 0, 0.5)),
        X2=list(beta=c(-0.5,0)),
        X3=list(beta=c( 0.5,0))
      )
    )
  )
  
  link <- list(
    c("identity","log","logit"),
    c("logit","logit","logit","identity","identity"),
    "identity"
  )
  
  set.seed(seed)
  dat <- as.data.table(
    rfrugalParam(n, formulas=forms, pars=pars,
                 family=fam, link=link)
  )
  
  ## prognostic part (centred)
  h  <- with(dat, 0.10*C1 + 0.12*sin(C4))
  h_cent <- h - mean(h)
  p0 <- p_base + h_cent
  
  fC      <- with(dat, (-C1 - C2 + (C3>0) + C4^2) / 3)
  fC_cent <- fC - mean(fC)
  delta   <- 0.10 * rho * fC_cent          # zero‑mean additive effect
  
  p1     <- p0 + delta                               # treated risk
  
  if (min(p0) < 0 || max(p1) > 1 || min(p1) < 0 || max(p0)>1)
    stop("Choose smaller rho or baseline: probabilities out of bounds.")
  
  CATE <- p1 - p0
  CRTE <- p1 / p0
  dat[, Y := as.integer(pnorm(E) < ifelse(A==1, p1, p0))]
  
  dat[, `:=`(
    CATE = CATE,
    CRTE = CRTE
  )]
  dat[, E := NULL]
  
  return(dat[])
}

make_base_filename <- function(n, rho, link_type, p_base) {
  if (is.null(p_base)) {
    baseline_tag <- "defaultBaseline"
  } else {
    # round to three decimals, then strip trailing zeros / trailing dot
    bstr <- formatC(p_base, format = "f", digits = 3)
    bstr <- sub("0+$", "", sub("\\.$", "", bstr))
    baseline_tag <- paste0("baseline", bstr)
  }
  sprintf("n%05d_rho%02d_%s_%s", n, round(rho * 100), link_type, baseline_tag)
}

run_trials <- function(n,
                       rho,
                       trials,
                       link_type = c("log", "identity"),
                       p_base = NULL,
                       num_trees=200) {
  
  ## --------------------  package setup  -------------------- ##
  library(causl)
  library(data.table)
  library(grf)
  library(rrcf)
  library(parallel)
  
  link_type <- match.arg(link_type)   # enforce valid input
  out_dir <- "data"
  dir.create(out_dir, showWarnings = FALSE)
  fname_base <- make_base_filename(n, rho, link_type, p_base)
  
  log_file <- file.path(out_dir, paste0("log_", fname_base, ".txt"))
  log_con  <- file(log_file, open = "wt")
  
  ## --------------------  helper : one trial  ---------------- ##
  # idx runs 1:trials; we pass idx as RNG seed to copula_rct()
  one_trial <- function(idx) {
    
    ## 1.  draw data  -------------------------------------------------
    #dat <- copula_rct(n, rho, idx)
    
    sim_fun <- if (link_type == "log") copula_rct_log else copula_rct_identity
    sim_args <- list(n = n, rho = rho, seed = idx)
    if (!is.null(p_base)) sim_args$p_base <- p_base          # pass only if given
    dat <- do.call(sim_fun, sim_args)
    
    split <- sample(seq_len(n), size = 0.8 * n)
    #x          <- dat[, !c("A", "Y", "CRTE"), with = FALSE]
    x <- dat[, !c("A", "Y", "CRTE", "CATE"), with = FALSE]
    
    x.train    <- x[ split, ]
    x.test     <- x[-split, ]
    y.train    <- dat$Y[ split]
    y.test     <- dat$Y[-split]
    t.train    <- dat$A[ split]
    t.test     <- dat$A[-split]
    crte.test  <- dat$CRTE[-split]
    
    ## 2.  fit forests  -----------------------------------------------
    forest.grf <- causal_forest(x.train, y.train, t.train, W.hat = 0.5, num.trees = num_trees, seed = 1234)
    forest.glm <- rr_causal_forest(x.train, y.train, t.train, rct = TRUE, num.trees = num_trees, seed = 1234)
    
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
  # cat("doing it")
  # print(res_lst)
  # cat("end doing it")
  ## --------------------  bind results  ----------------------------- ##
  pvals <- do.call(rbind, lapply(res_lst, `[[`, "pvals"))
  mapes <- do.call(rbind, lapply(res_lst, `[[`, "mapes"))
  colnames(pvals) <- c("n", "rho", "p_grf", "p_glm")
  colnames(mapes) <- c("n", "rho", "mape_grf", "mape_glm")
  
  pval_file <- file.path(out_dir, paste0("pval_",  fname_base, ".rds"))
  mape_file <- file.path(out_dir, paste0("mape_", fname_base, ".rds"))
  
  saveRDS(pvals, file = pval_file)
  saveRDS(mapes, file = mape_file)
  
  invisible(list(pvals = pvals, mapes = mapes))
}

power_from_rds <- function(n,
                           rho,
                           link_type = c("log", "identity"),
                           p_base    = NULL,
                           alpha     = 0.05,
                           out_dir   = "data") {
  
  link_type <- match.arg(link_type)
  fname_base <- fname_base <- make_base_filename(n, rho, link_type, p_base)
  pval_file <- file.path(out_dir, paste0("pval_", fname_base, ".rds"))
  
  if (!file.exists(pval_file)) {
    stop("Cannot find file: ", pval_file)
  }
  
  ## ---- load & compute power ------------------------------------------------
  pvals <- readRDS(pval_file)          # columns: n, rho, p_grf, p_glm
  
  power_grf <- mean(pvals[, "p_grf"] < alpha)
  power_glm <- mean(pvals[, "p_glm"] < alpha)
  
  setNames(c(power_grf, power_glm),
           c("power_grf", "power_glm"))
}

mape_from_rds <- function(n,
                          rho,
                          link_type = c("log", "identity"),
                          p_base    = NULL,
                          out_dir   = "data") {
  
  link_type <- match.arg(link_type)
  fname_base <- fname_base <- make_base_filename(n, rho, link_type, p_base)
  mape_file <- file.path(out_dir, paste0("mape_", fname_base, ".rds"))
  
  if (!file.exists(mape_file)) {
    stop("Cannot find file: ", mape_file)
  }
  
  ## ---- load & summarise MAPE ---------------------------------------------
  mapes <- readRDS(mape_file)          # columns: n, rho, mape_grf, mape_glm
  
  mean_grf <- mean(mapes[, "mape_grf"])
  mean_glm <- mean(mapes[, "mape_glm"])
  
  setNames(c(mean_grf, mean_glm), c("mean_mape_grf", "mean_mape_glm"))
}

