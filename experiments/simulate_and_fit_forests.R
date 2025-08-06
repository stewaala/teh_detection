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

rd_predict <- function(object, newdata) {
  X <- object[["X.orig"]]
  Y <- object[["Y.orig"]]
  W <- object[["W.orig"]]

  length = dim(newdata)[1]
  output = numeric(length)
  for(i in 1:length){
    fw = get_forest_weights(object, newdata[i,])
    df = cbind.data.frame(W, fw[1:length(Y)], fw[1:length(Y)]*Y)
    tau1 = sum(subset(df, df[,1] == 1)[,3])/sum(subset(df, df[,1] == 1)[,2])
    tau0 = sum(subset(df, df[,1] == 0)[,3])/sum(subset(df, df[,1] == 0)[,2])
    output[i] = tau1-tau0
  }
  output
}

run_trials <- function(n,
                       rho,
                       trials,
                       link_type = c("log", "identity"),
                       p_base = NULL,
                       num_trees = 200) {
  
  ## ---  packages quietly --------------------------------------------
  suppressMessages({
    library(causl);  library(data.table)
    library(grf);    library(rrcf)
    library(parallel)
  })
  link_type <- match.arg(link_type)
  
  out_dir <- "data"; dir.create(out_dir, showWarnings = FALSE)
  fname_base <- make_base_filename(n, rho, link_type, p_base)
  
  ## --- (A) open log file --------------------------------------------
  log_file <- file.path(out_dir, paste0("log_", fname_base, ".txt"))
  log_con  <- file(log_file, open = "wt")
  on.exit(close(log_con), add = TRUE)
  
  ## --- helper --------------------------------------------------------
  one_trial <- function(idx) {
    tryCatch({
      
      ## 1. simulate ---------------------------------------------------
      sim_fun  <- if (link_type == "log") copula_rct_log else copula_rct_identity
      sim_args <- list(n = n, rho = rho, seed = idx)
      if (!is.null(p_base)) sim_args$p_base <- p_base
      dat <- suppressWarnings(do.call(sim_fun, sim_args))   # keep console quiet
      
      ## 2. split / fit -----------------------------------------------
      split <- sample.int(n, size = 0.8 * n)
      x <- dat[, !c("A", "Y", "CRTE", "CATE"), with = FALSE]
      x.train <- x[ split, ];  x.test <- x[-split, ]
      y.train <- dat$Y[ split]; y.test <- dat$Y[-split]
      t.train <- dat$A[ split]; t.test <- dat$A[-split]
      crte.test <- dat$CRTE[-split]
      
      forest.grf <- suppressMessages(
        causal_forest(x.train, y.train, t.train, W.hat = 0.5, num.trees = num_trees, seed = 1234)
      )
      forest.glm <- suppressMessages(
        rr_causal_forest(x.train, y.train, t.train, rct = TRUE, num.trees = num_trees, seed = 1234)
      )
      
      ## 3. predictions & metrics -------------------------------------
      rr_pred.grf <- rr_predict(forest.grf, x.test)
      rr_pred.glm <- rr_predict(forest.glm, x.test)
      
      rd_pred.grf <- predict(forest.grf, x.test)$predictions
      rd_pred.glm <- rd_predict(forest.glm, x.test)
 
      dat_test <- data.table::copy(dat[-split])                # only test rows
      dat_test[, `:=`(
        crte_hat_grf = rr_pred.grf,
        crte_hat_glm = rr_pred.glm,
        cate_hat_grf = rd_pred.grf,
        cate_hat_glm = rd_pred.glm,
        trial_id     = idx
      )]
      
      list(dat = dat_test)
    }, error = function(e) {
      cat(sprintf("trial %d: %s\n", idx, e$message), file = log_con)
      return(NULL)                 # skips this trial
    })
  }
  
  ## --- parallel execution -------------------------------------------
  res_lst <- mclapply(seq_len(trials), one_trial,
                      mc.cores = parallel::detectCores() - 1,
                      mc.preschedule = FALSE)
  
  ## --- bind, removing NULLs -----------------------------------------
  ok      <- vapply(res_lst, is.null, logical(1), USE.NAMES = FALSE)
  good    <- res_lst[!ok]
  
  if (length(good) == 0) stop("All trials failed; see ", log_file)
  
  test_bundle <- data.table::rbindlist(lapply(good, `[[`, "dat"))
  bundle_file <- file.path(out_dir, paste0("bundle_", fname_base, ".rds"))
  saveRDS(test_bundle, file = bundle_file)
  
  if (file.size(log_file) == 0) unlink(log_file)   # delete empty log
  else message("Finished with some errors. See ", log_file)
  
  invisible(NULL)
}

get_covariate_cols <- function(dt) {
  drop_cols <- c("Y", "A",
                 "CRTE", "CATE",
                 "crte_hat_grf", "crte_hat_glm",
                 "cate_hat_grf", "cate_hat_glm",
                 "trial_id")
  setdiff(names(dt), drop_cols)
}

power_from_bundle <- function(n,
                              rho,
                              link_type = c("log", "identity"),
                              p_base    = NULL,
                              alpha     = 0.05,
                              out_dir   = "data") {
  
  link_type <- match.arg(link_type)
  fname_base <- make_base_filename(n, rho, link_type, p_base)
  bundle_file <- file.path(out_dir, paste0("bundle_", fname_base, ".rds"))
  if (!file.exists(bundle_file))
    stop("Cannot find bundle file: ", bundle_file)
  
  dt <- readRDS(bundle_file)
  cov_cols <- get_covariate_cols(dt)
  
  trial_ids <- unique(dt$trial_id)
  p_grf <- numeric(length(trial_ids))
  p_glm <- numeric(length(trial_ids))
  
  for (i in seq_along(trial_ids)) {
    sub <- dt[trial_id == trial_ids[i]]
    
    y.test <- sub$Y
    t.test <- sub$A                     # <-- define it, as in run_trials()
    
    anova_df <- data.frame(
      y.test = y.test,
      t.test = t.test,
      sub[, ..cov_cols]
    )
    
    base_mod <- glm(y.test ~ ., family = poisson, data = anova_df)
    
    add_grf <- glm(
      y.test ~ .,
      family = poisson,
      data = cbind(anova_df,
                   t.test * log(sub$crte_hat_grf))
    )
    add_glm <- glm(
      y.test ~ .,
      family = poisson,
      data = cbind(anova_df,
                   t.test * log(sub$crte_hat_glm))
    )
    
    p_grf[i] <- 1 - pchisq(anova(base_mod, add_grf)$Deviance[2], 1)
    p_glm[i] <- 1 - pchisq(anova(base_mod, add_glm)$Deviance[2], 1)
  }
  
  setNames(c(mean(p_grf < alpha), mean(p_glm < alpha)),
           c("power_grf", "power_glm"))
}

mape_from_bundle <- function(n,
                             rho,
                             link_type = c("log", "identity"),
                             p_base    = NULL,
                             out_dir   = "data") {
  
  link_type  <- match.arg(link_type)
  fname_base <- make_base_filename(n, rho, link_type, p_base)
  bundle_file <- file.path(out_dir, paste0("bundle_", fname_base, ".rds"))
  
  if (!file.exists(bundle_file))
    stop("Cannot find bundle file: ", bundle_file)
  
  ## ---- load test‑set bundle ----------------------------------------
  dt <- readRDS(bundle_file)   # columns include CRTE, crte_hat_grf, crte_hat_glm, trial_id
  
  ## ---- per‑trial MAPE ---------------------------------------------
  MAPE <- dt[, .(
    mape_grf = mean(abs((CRTE - crte_hat_grf) / CRTE)),
    mape_glm = mean(abs((CRTE - crte_hat_glm) / CRTE))
  ),
  by = trial_id
  ]
  
  ## ---- average over trials (matches mape_from_rds) ----------------
  mean_grf <- mean(MAPE$mape_grf)
  mean_glm <- mean(MAPE$mape_glm)
  
  setNames(c(mean_grf, mean_glm),
           c("mean_mape_grf", "mean_mape_glm"))
}

power_rd_from_bundle <- function(n,
                                 rho,
                                 link_type = c("log", "identity"),
                                 p_base    = NULL,
                                 alpha     = 0.05,
                                 out_dir   = "data") {
  
  link_type <- match.arg(link_type)
  fname_base <- make_base_filename(n, rho, link_type, p_base)
  bundle_file <- file.path(out_dir, paste0("bundle_", fname_base, ".rds"))
  if (!file.exists(bundle_file))
    stop("Cannot find bundle file: ", bundle_file)
  
  ## ------------------ load bundle ------------------ ##
  dt <- readRDS(bundle_file)
  cov_cols <- get_covariate_cols(dt)
  
  trial_ids <- unique(dt$trial_id)
  p_grf <- numeric(length(trial_ids))
  p_glm <- numeric(length(trial_ids))
  
  for (i in seq_along(trial_ids)) {
    sub <- dt[trial_id == trial_ids[i]]
    
    ## ---------- construct centred interaction term ----------
    tau_bar_grf <- mean(sub$cate_hat_grf)
    tau_bar_glm <- mean(sub$cate_hat_glm)
    
    int_grf <- sub$A * (sub$cate_hat_grf - tau_bar_grf)
    int_glm <- sub$A * (sub$cate_hat_glm - tau_bar_glm)
    
    ## ---------- baseline model: Y ~ X + W -------------------
    base_df <- data.frame(
      y  = sub$Y,
      W  = sub$A,
      sub[, ..cov_cols]
    )
    base_mod <- lm(y ~ ., data = base_df)  # linear‑probability calibration
    
    ## ---------- augmented with GRF interaction --------------
    add_grf <- lm(y ~ .,
                  data = cbind(base_df, int_grf))
    
    ## ---------- augmented with GLM interaction --------------
    add_glm <- lm(y ~ .,
                  data = cbind(base_df, int_glm))
    
    ## ---------- 1‑df F‑tests --------------------------------
    p_grf[i] <- anova(base_mod, add_grf)$`Pr(>F)`[2]
    p_glm[i] <- anova(base_mod, add_glm)$`Pr(>F)`[2]
  }
  
  setNames(c(mean(p_grf < alpha), mean(p_glm < alpha)),
           c("power_grf", "power_glm"))
}

# ------------------------------------------------------------
# Run a grid of (n, ρ) settings for a fixed link_type & p_base
# ------------------------------------------------------------
run_trials_grid <- function(n_values,
                            rho_values,
                            link_type = c("log", "identity"),
                            p_base    = NULL,
                            trials    = 100,
                            num_trees = 200,
                            out_dir   = "data") {
  
  # example usage:
  #
  # run_trials_grid(
  #   n_values  = c(2500, 5000, 7500, 10000),
  #   rho_values = c(0.00, 0.25, 0.50, 0.75),
  #   link_type  = "log",
  #   p_base     = NULL   # default baseline
  # )
  
  link_type <- match.arg(link_type)          # validate input
  dir.create(out_dir, showWarnings = FALSE)  # ensure folder exists
  
  status <- list()                           # bookkeeping
  
  for (n   in n_values) {
    for (rho in rho_values) {
      
      # ---- does bundle already exist? ----------------------
      fname_base  <- make_base_filename(n, rho, link_type, p_base)
      bundle_file <- file.path(out_dir, paste0("bundle_", fname_base, ".rds"))
      
      if (file.exists(bundle_file)) {
        message(sprintf("SKIP  n=%d rho=%.2f (%s) — bundle exists",
                        n, rho, fname_base))
        status[[fname_base]] <- "skipped"
        next
      }
      
      # ---- run simulations + forests -----------------------
      message(sprintf("RUN   n=%d rho=%.2f (%s)", n, rho, fname_base))
      tryCatch({
        run_trials(n        = n,
                   rho      = rho,
                   trials   = trials,
                   link_type = link_type,
                   p_base   = p_base,
                   num_trees = num_trees)
        status[[fname_base]] <- "completed"
      }, error = function(e) {
        warning(sprintf("FAILED n=%d rho=%.2f : %s", n, rho, e$message))
        status[[fname_base]] <- paste("failed:", e$message)
      })
    }
  }
  
  invisible(status)   # return a named list of outcomes
}

