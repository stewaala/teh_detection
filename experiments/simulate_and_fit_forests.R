library(causl)
library(data.table)

.EEXP_H_CONST <- 1.20332119584
.LOW_BASELINE_THRESHOLD <- 0.02  # treat p_base <= 2% as "low baseline"
.DEFAULT_PBASE_LOG <- .EEXP_H_CONST * exp(-2)  # ≈ 0.163
.DEFAULT_PBASE_ID  <- 0.23
.TARGET_TEST_EVENTS <- 40
.TEST_FRAC_DEFAULT  <- 0.20
.TEST_FRAC_CAP      <- 0.50
.LOW_N_TRAIN_THRESHOLD <- 2000

choose_test_frac <- function(n,
                             link_type,
                             p_base,
                             default = .TEST_FRAC_DEFAULT,
                             target_events = .TARGET_TEST_EVENTS,
                             cap = .TEST_FRAC_CAP,
                             floor = 0.20) {
  # effective baseline risk if user didn't override
  p_eff <- if (is.null(p_base)) {
    if (link_type == "log") .DEFAULT_PBASE_LOG else .DEFAULT_PBASE_ID
  } else p_base
  
  # If default split already yields enough test events, keep it
  exp_events_default <- n * default * p_eff
  if (exp_events_default >= target_events) return(default)
  
  # Otherwise, enlarge split just enough (but not beyond cap)
  tf <- target_events / (n * max(p_eff, 1e-9))
  tf <- min(cap, max(floor, tf))
  
  if (tf > default + 1e-9) {
    message(sprintf(
      "auto test_frac=%.2f (target ~%d test events; p≈%.3f, n=%d)",
      tf, target_events, p_eff, n
    ))
  }
  tf
}

choose_forest_hparams <- function(link_type, p_base, num_trees, n_train) {
  ultra_low <- (!is.null(p_base)) && (link_type == "log") && (p_base <= .LOW_BASELINE_THRESHOLD)
  if (ultra_low) {
    return(list(
      grf  = list(num.trees = max(1000, num_trees),
                  min.node.size = 40,
                  sample.fraction = 0.95,
                  ci.group.size = 1L),
      rrcf = list(num.trees = max(1000, num_trees),
                  min.node.size = 40,
                  sample.fraction = 0.49)
    ))
  }
  
  # --- NEW: mild tweaks for small training sets (any link/baseline) ---
  if (n_train <= .LOW_N_TRAIN_THRESHOLD) {
    # scale leaf size with training size; keep sample.fraction default
    leaf <- max(10L, as.integer(round(0.03 * n_train)))  # ~3% of train
    return(list(
      grf  = list(num.trees = max(500, num_trees),
                  min.node.size = leaf),
      rrcf = list(num.trees = max(500, num_trees),
                  min.node.size = leaf)
    ))
  }
  
  # Default settings (covers n≥~2000 and 5% baseline case)
  list(
    grf  = list(num.trees = num_trees),
    rrcf = list(num.trees = num_trees)
  )
}

arg_filter <- function(fun, args) {
  fn_args <- names(formals(fun))
  args[names(args) %in% fn_args]
}

# # Pick a test split when baseline is ultra-low (≤ .LOW_BASELINE_THRESHOLD)
# choose_test_frac <- function(n,
#                              link_type,
#                              p_base,
#                              default = 0.20,
#                              target_events = 40,  # aim for ~40 test positives
#                              cap = 0.50,
#                              floor = 0.20) {
#   ultra_low <- (!is.null(p_base)) && (link_type == "log") && (p_base <= .LOW_BASELINE_THRESHOLD)
#   if (!ultra_low) return(default)
#   
#   # Expected test events ≈ n * test_frac * p_base  ⇒  test_frac ≈ target/(n*p_base)
#   tf <- target_events / (n * max(p_base, 1e-9))
#   tf <- min(cap, max(floor, tf))
#   
#   if (tf > default + 1e-9) {
#     message(sprintf("auto test_frac=%.2f to target ~%d test events at p_base=%.3f",
#                     tf, target_events, p_base))
#   }
#   tf
# }

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

# ------------------------------------------------------------------------------
# Format just the baseline tag
# ------------------------------------------------------------------------------
format_baseline_tag <- function(p_base) {
  if (is.null(p_base)) {
    "defaultBaseline"
  } else {
    # round to three decimals, then strip trailing zeros and trailing dot
    bstr <- formatC(p_base, format = "f", digits = 3)
    bstr <- sub("0+$", "", sub("\\.$", "", bstr))
    paste0("baseline", bstr)
  }
}

# ------------------------------------------------------------------------------
# Build the full filename base using the new helper
# ------------------------------------------------------------------------------
make_base_filename <- function(n, rho, link_type, p_base) {
  baseline_tag <- format_baseline_tag(p_base)
  sprintf("n%05d_rho%02d_%s_%s",
          n,
          round(rho * 100),
          link_type,
          baseline_tag)
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
                       num_trees = 200,
                       test_frac = NULL) {
  
  ## ---  packages quietly --------------------------------------------
  suppressMessages({
    library(causl);  library(data.table)
    library(grf);    library(rrcf)
    library(parallel)
  })
  link_type <- match.arg(link_type)
  
  out_dir <- "data"; dir.create(out_dir, showWarnings = FALSE)
  fname_base <- make_base_filename(n, rho, link_type, p_base)
  
  # after fname_base:
  test_frac_eff <- if (is.null(test_frac)) choose_test_frac(n, link_type, p_base) else test_frac
  n_train <- as.integer(round((1 - test_frac_eff) * n))  # optional: precompute once
  
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
      # # decide test split
      split <- sample.int(n, size = as.integer(round((1 - test_frac_eff) * n)))
      x <- dat[, !c("A", "Y", "CRTE", "CATE"), with = FALSE]
      x.train <- x[ split, ];  x.test <- x[-split, ]
      y.train <- dat$Y[ split]; y.test <- dat$Y[-split]
      t.train <- dat$A[ split]; t.test <- dat$A[-split]
      crte.test <- dat$CRTE[-split]

      # choose hparams
      hp <- choose_forest_hparams(link_type, p_base, num_trees, n_train)
      
      # GRF
      grf_args <- c(
        list(X = x.train, Y = y.train, W = t.train, W.hat = 0.5, seed = 1234),
        arg_filter(causal_forest, hp$grf)
      )
      forest.grf <- tryCatch(
        suppressMessages(do.call(causal_forest, grf_args)),
        error = function(e) {
          cat(sprintf("trial %d [GRF]: %s\n", idx, e$message), file = log_con)
          stop(e)
        }
      )
      
      # RRCF
      rrcf_args <- c(
        list(X = x.train, Y = y.train, W = t.train, rct = TRUE, seed = 1234),
        arg_filter(rr_causal_forest, hp$rrcf)
      )
      forest.glm <- tryCatch(
        suppressMessages(do.call(rr_causal_forest, rrcf_args)),
        error = function(e) {
          cat(sprintf("trial %d [RRCF]: %s\n", idx, e$message), file = log_con)
          stop(e)
        }
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

load_bundle <- function(n, rho, link_type=c("log", "identity"), p_base=NULL) {
  
  out_dir = "data"
  link_type <- match.arg(link_type)
  fname_base <- make_base_filename(n, rho, link_type, p_base)
  bundle_file <- file.path(out_dir, paste0("bundle_", fname_base, ".rds"))
  if (!file.exists(bundle_file))
    stop("Cannot find bundle file: ", bundle_file)
  
  readRDS(bundle_file)
}

stable_log <- function(rr, cap_q = 0.995, cap_abs = 12) {
  x <- log(rr)
  x[!is.finite(x)] <- NA_real_
  if (all(!is.finite(x))) return(x)
  Lq <- stats::quantile(abs(x[is.finite(x)]), cap_q, na.rm = TRUE)
  L  <- min(Lq, cap_abs)  # cap to, say, |log RR| <= 12 (~ e^12 ≈ 1.6e5)
  pmax(pmin(x,  L), -L)
}

power_from_bundle <- function(n, rho, link_type = c("log","identity"),
                              p_base = NULL, alpha = 0.05) {
  dt <- load_bundle(n, rho, link_type, p_base)
  cov_cols <- get_covariate_cols(dt)
  
  tri <- unique(dt$trial_id)
  p_grf <- rep(NA_real_, length(tri))
  p_glm <- rep(NA_real_, length(tri))
  
  ctrl <- glm.control(maxit = 100, epsilon = 1e-10)
  
  for (i in seq_along(tri)) {
    sub <- dt[trial_id == tri[i]]
    base_df <- data.frame(y = sub$Y, W = sub$A, sub[, ..cov_cols])
    
    ## ---- GRF path ----
    lg <- stable_log(sub$crte_hat_grf)
    mask <- is.finite(lg) & stats::complete.cases(base_df)
    if (sum(mask) >= 5 && length(unique(sub$A[mask])) == 2) {
      df <- cbind(base_df[mask, , drop = FALSE], log_crte_hat = lg[mask])
      
      # (optional) ensure a few positives in each arm to avoid separation
      pos1 <- sum(df$y == 1 & df$W == 1)
      pos0 <- sum(df$y == 1 & df$W == 0)
      if (pos1 >= 5 && pos0 >= 5) {
        base_g <- try(glm(y ~ . - log_crte_hat, family = poisson, data = df, control = ctrl), silent = TRUE)
        if (!inherits(base_g, "try-error")) {
          add_g  <- try(glm(y ~ . - log_crte_hat + I(W * log_crte_hat),
                            family = poisson, data = df, control = ctrl,
                            start = c(coef(base_g), 0)), silent = TRUE)
          if (!inherits(add_g, "try-error")) {
            dev_g <- anova(base_g, add_g)$Deviance
            p_grf[i] <- 1 - pchisq(dev_g[2], df = 1)
          }
        }
      }
    }
    
    ## ---- RRCF path (mirror) ----
    lgm <- stable_log(sub$crte_hat_glm)
    maskm <- is.finite(lgm) & stats::complete.cases(base_df)
    if (sum(maskm) >= 5 && length(unique(sub$A[maskm])) == 2) {
      dfm <- cbind(base_df[maskm, , drop = FALSE], log_crte_hat = lgm[maskm])
      pos1 <- sum(dfm$y == 1 & dfm$W == 1)
      pos0 <- sum(dfm$y == 1 & dfm$W == 0)
      if (pos1 >= 5 && pos0 >= 5) {
        base_m <- try(glm(y ~ . - log_crte_hat, family = poisson, data = dfm, control = ctrl), silent = TRUE)
        if (!inherits(base_m, "try-error")) {
          add_m  <- try(glm(y ~ . - log_crte_hat + I(W * log_crte_hat),
                            family = poisson, data = dfm, control = ctrl,
                            start = c(coef(base_m), 0)), silent = TRUE)
          if (!inherits(add_m, "try-error")) {
            dev_m <- anova(base_m, add_m)$Deviance
            p_glm[i] <- 1 - pchisq(dev_m[2], df = 1)
          }
        }
      }
    }
  }
  
  setNames(c(mean(p_grf < alpha, na.rm = TRUE),
             mean(p_glm < alpha, na.rm = TRUE)),
           c("power_grf", "power_glm"))
}

# mape_from_bundle <- function(n,
#                              rho,
#                              link_type = c("log", "identity"),
#                              p_base    = NULL,
#                              out_dir   = "data") {
#   
#   # link_type  <- match.arg(link_type)
#   # fname_base <- make_base_filename(n, rho, link_type, p_base)
#   # bundle_file <- file.path(out_dir, paste0("bundle_", fname_base, ".rds"))
#   # 
#   # if (!file.exists(bundle_file))
#   #   stop("Cannot find bundle file: ", bundle_file)
#   # 
#   # ## ---- load test‑set bundle ----------------------------------------
#   # dt <- readRDS(bundle_file)   # columns include CRTE, crte_hat_grf, crte_hat_glm, trial_id
#   dt <- load_bundle(n, rho, link_type, p_base)
#   
#   ## ---- per‑trial MAPE ---------------------------------------------
#   MAPE <- dt[, .(
#     mape_grf = mean(abs((CRTE - crte_hat_grf) / CRTE)),
#     mape_glm = mean(abs((CRTE - crte_hat_glm) / CRTE))
#   ),
#   by = trial_id
#   ]
#   
#   ## ---- average over trials (matches mape_from_rds) ----------------
#   mean_grf <- mean(MAPE$mape_grf)
#   mean_glm <- mean(MAPE$mape_glm)
#   
#   setNames(c(mean_grf, mean_glm),
#            c("mean_mape_grf", "mean_mape_glm"))
# }

mape_from_bundle <- function(n,
                             rho,
                             link_type = c("log", "identity"),
                             p_base    = NULL) {
  dt <- load_bundle(n, rho, link_type, p_base)
  
  per_trial <- dt[, {
    num   <- abs(CRTE - crte_hat_grf) / CRTE
    keepg <- is.finite(num)
    num_m <- abs(CRTE - crte_hat_glm) / CRTE
    keepm <- is.finite(num_m)
    
    .( mape_grf = if (sum(keepg) >= 10) mean(num[keepg]) else NA_real_,
       mape_glm = if (sum(keepm) >= 10) mean(num_m[keepm]) else NA_real_ )
  }, by = trial_id]
  
  setNames(c(mean(per_trial$mape_grf, na.rm = TRUE),
             mean(per_trial$mape_glm, na.rm = TRUE)),
           c("mean_mape_grf", "mean_mape_glm"))
}


power_rd_from_bundle <- function(n,
                                 rho,
                                 link_type = c("log", "identity"),
                                 p_base    = NULL,
                                 alpha     = 0.05) {
  dt <- load_bundle(n, rho, link_type, p_base)
  cov_cols <- get_covariate_cols(dt)
  
  trial_ids <- unique(dt$trial_id)
  p_grf <- rep(NA_real_, length(trial_ids))
  p_glm <- rep(NA_real_, length(trial_ids))
  
  rd_pval <- function(base_df, tau_hat) {
    mask <- is.finite(tau_hat) & stats::complete.cases(base_df)
    if (sum(mask) < 10) return(NA_real_)
    Wm <- base_df$W[mask]
    tab <- table(Wm)
    if (length(tab) < 2 || any(tab < 5)) return(NA_real_)  # need both arms, min 5 each
    
    tau_c <- tau_hat[mask] - mean(tau_hat[mask])
    int   <- Wm * tau_c
    if (sd(int) == 0) return(NA_real_)                     # no variation ⇒ no test
    
    df_b <- base_df[mask, , drop = FALSE]
    df   <- cbind(df_b, int = int)
    
    base <- try(lm(y ~ ., data = df_b), silent = TRUE)
    add  <- try(lm(y ~ . + int, data = df), silent = TRUE)
    if (inherits(base, "try-error") || inherits(add, "try-error")) return(NA_real_)
    
    # If adding 'int' didn’t increase rank, ANOVA p-value is meaningless
    if (ncol(model.matrix(add)) == ncol(model.matrix(base))) return(NA_real_)
    
    av <- anova(base, add)
    p  <- suppressWarnings(av$`Pr(>F)`[2])
    if (!is.finite(p)) NA_real_ else p
  }
  
  for (i in seq_along(trial_ids)) {
    sub <- dt[trial_id == trial_ids[i]]
    base_df <- data.frame(y = sub$Y, W = sub$A, sub[, ..cov_cols])
    
    p_grf[i] <- rd_pval(base_df, sub$cate_hat_grf)
    p_glm[i] <- rd_pval(base_df, sub$cate_hat_glm)
  }
  
  setNames(c(mean(p_grf < alpha, na.rm = TRUE),
             mean(p_glm < alpha, na.rm = TRUE)),
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

