library(causl)
library(data.table)
library(sandwich)
library(lmtest)
library(brglm2)
library(grf)
library(rrcf)
library(parallel)


.EEXP_H_CONST <- 1.20332119584

# ------------------------------------------------------------------------------
# The following function replicates the functionality of copula_rct() in 
# Shirvaikar's github, but it takes an additional p_base parameter that
# allows the user the specify the level of the expected untreated risk
# ------------------------------------------------------------------------------
copula_rct_log <- function(n, rho = 0, seed = 111, p_base = .EEXP_H_CONST * exp(-2)){
  
  alpha <- log(p_base) - log(.EEXP_H_CONST)
  #message("alpha set to ", alpha, '\n')
  
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

make_and_save_bundle <- function(n, rho, link_type = c("log", "identity"), num_trials=100, p_base = NULL, num_trees = 500, test_frac = 0.2, skip_if_file_exists=TRUE) {

  out_dir <- "data"
  fname_base <- make_base_filename(n, rho, link_type, p_base)
  bundle_file <- file.path(out_dir, paste0("bundle_", fname_base, ".rds"))

  if (file.exists(bundle_file) && skip_if_file_exists) {
    return()
  }
  bundle <- make_bundle(n, rho, link_type, num_trials, p_base, num_trees, test_frac)
  saveRDS(bundle, file = bundle_file)
  
}

make_bundle <- function(n, rho, link_type = c("log", "identity"), num_trials=100, p_base = NULL, num_trees = 500, test_frac = 0.2) {
  
  link_type <- match.arg(link_type)

  one_trial <- function(idx) {

    sim_fun  <- if (link_type == "log") copula_rct_log else copula_rct_identity
    sim_args <- list(n = n, rho = rho, seed = idx)
    if (!is.null(p_base)) sim_args$p_base <- p_base
    dat <- suppressWarnings(do.call(sim_fun, sim_args))
    
    x <- dat[, !c("A", "Y", "CRTE", "CATE"), with = FALSE]
    set.seed(idx)
    split <- sample(seq_len(n), size = n*(1 - test_frac))
    x.train <- x[ split, ];  x.test <- x[-split, ]
    y.train <- dat$Y[ split]; y.test <- dat$Y[-split]
    t.train <- dat$A[ split]; t.test <- dat$A[-split]
    crte.test <- dat$CRTE[-split]
    
    forest.grf <- causal_forest(X = x.train, Y = y.train, W = t.train, W.hat = 0.5, seed = 1234, num.trees=num_trees)
    forest.glm = rr_causal_forest(x.train, y.train, t.train, rct=TRUE, seed=1234, num.trees=num_trees)
    
    rr_pred.grf <- rr_predict(forest.grf, x.test)
    rr_pred.glm <- rr_predict(forest.glm, x.test)
    
    rd_pred.grf <- predict(forest.grf, x.test)$predictions
    rd_pred.glm <- rd_predict(forest.glm, x.test)
    
    dat_test <- data.table::copy(dat[-split])
    dat_test[, `:=`(
      crte_hat_grf = rr_pred.grf,
      crte_hat_glm = rr_pred.glm,
      cate_hat_grf = rd_pred.grf,
      cate_hat_glm = rd_pred.glm,
      trial_id     = idx
    )]
  }
  
  res_lst <- mclapply(seq_len(num_trials), one_trial, mc.cores = parallel::detectCores() - 1, mc.preschedule = FALSE)
  bundle <- do.call(rbind, res_lst)
  bundle
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

poisson_omnibus_pval <- function(pred, y.test, t.test, x.test, small_sample=FALSE) {
  
  anova.data = data.frame(cbind(y.test, t.test, x.test))
  
  if (small_sample) {
    if (length(unique(pred)) == 1) {
      pval <- 1
    } else {
      dt <- cbind(anova.data, t_log_pred = t.test*log(pred))
      ctrl <- brglm2::brglm_control(maxit = 2000, epsilon = 1e-10, slowit = 1)
      fit  <- glm(y.test ~ ., family = poisson(link="log"), data = dt, method = "brglmFit", control = ctrl)
      pval <- coeftest(fit, vcov = vcovHC(fit, type = "HC3"))["t_log_pred", "Pr(>|z|)"]
    }
  } else {
    model.base = glm(y.test ~ ., family = poisson, data = anova.data)
    model.hte = glm(y.test ~ ., family = poisson, data = cbind(anova.data, t.test*log(pred)))
    anova.grf = anova(model.base, model.hte)
    pval = 1 - pchisq(anova.grf$Deviance[2], df = 1)
  }
  return(pval)
}

# poisson_power_from_bundle <- function(bundle, forest_type = c("grf", "rrcf"), small_sample=FALSE) {
# 
#   forest_type <- match.arg(forest_type)
#   cov_cols <- get_covariate_cols(bundle)
#   trials <- unique(bundle$trial_id)
#   
#   get_pval <- function(trial) {
#     dt = bundle[bundle$trial_id == trial]
#     x.test <- dt[, cov_cols, with = FALSE]
#     pred <- if (forest_type == "grf") dt$crte_hat_grf else dt$crte_hat_glm
#     poisson_omnibus_pval(pred, dt$Y, dt$A, x.test, small_sample)
#   }
#   pvals <- vapply(trials, FUN=get_pval, FUN.VALUE=numeric(1))
#   mean(pvals < 0.05)
# }
# New: returns estimate + CI + counts
poisson_power_from_bundle_detail <- function(bundle,
                                             forest_type = c("grf", "rrcf"),
                                             small_sample = FALSE,
                                             alpha_test  = 0.05,
                                             conf_level  = 0.95,
                                             ci_method   = c("wilson", "jeffreys", "exact")) {
  forest_type <- match.arg(forest_type)
  ci_method   <- match.arg(ci_method)
  
  cov_cols <- get_covariate_cols(bundle)
  trials   <- unique(bundle$trial_id)
  
  get_pval <- function(trial) {
    dt <- bundle[bundle$trial_id == trial]
    x.test <- dt[, cov_cols, with = FALSE]
    pred <- if (forest_type == "grf") dt$crte_hat_grf else dt$crte_hat_glm
    poisson_omnibus_pval(pred, dt$Y, dt$A, x.test, small_sample)
  }
  
  pvals <- vapply(trials, FUN = get_pval, FUN.VALUE = numeric(1))
  n <- length(pvals)
  k <- sum(pvals < alpha_test)
  p_hat <- k / n
  
  alpha <- 1 - conf_level
  if (ci_method == "wilson") {
    z <- qnorm(1 - alpha/2)
    denom  <- 1 + (z^2)/n
    center <- (p_hat + (z^2)/(2*n)) / denom
    halfw  <- (z * sqrt(p_hat*(1 - p_hat)/n + (z^2)/(4*n^2))) / denom
    lower  <- max(0, center - halfw)
    upper  <- min(1, center + halfw)
  } else if (ci_method == "jeffreys") {
    lower <- qbeta(alpha/2,     k + 0.5, n - k + 0.5)
    upper <- qbeta(1 - alpha/2, k + 0.5, n - k + 0.5)
  } else { # exact (Clopper–Pearson)
    lower <- if (k == 0) 0 else qbeta(alpha/2,     k,     n - k + 1)
    upper <- if (k == n) 1 else qbeta(1 - alpha/2, k + 1, n - k)
  }
  
  list(
    estimate   = p_hat,
    ci_lower   = lower,
    ci_upper   = upper,
    n_trials   = n,
    k_reject   = k,
    alpha_test = alpha_test,
    conf_level = conf_level,
    ci_method  = ci_method
    # , pvals = pvals  # uncomment if you want to return raw p-values
  )
}

# Original name: remains backward-compatible (still returns a numeric)
poisson_power_from_bundle <- function(bundle,
                                      forest_type = c("grf", "rrcf"),
                                      small_sample = FALSE) {
  res <- poisson_power_from_bundle_detail(bundle,
                                          forest_type = forest_type,
                                          small_sample = small_sample)
  res$estimate
}


# mape_from_bundle <- function(bundle, forest_type = c("grf", "rrcf"), trial = NULL) {
#   
#   data <- bundle
#   if (!is.null(trial)) 
#     bundle <- bundle[bundle$trial_id == trial,]
#   
#   pred <- if (forest_type == "grf") bundle$crte_hat_grf else bundle$crte_hat_glm
#   mean(abs((bundle$CRTE - pred)/bundle$CRTE))
# }
# New: detailed MAPE with confidence intervals
# Robust, CI-returning helper (backward-compatible with earlier "detail" API)
mape_from_bundle_detail <- function(bundle,
                                    forest_type = c("grf", "rrcf"),
                                    trial = NULL,
                                    conf_level = 0.95,
                                    ci_method  = c("t", "bootstrap"),
                                    B = 2000,
                                    # NEW: robustness knobs
                                    zero_handling     = c("safe", "epsilon", "drop"),
                                    eps               = 1e-10,
                                    na_pred_handling  = c("drop", "error")) {
  
  forest_type     <- match.arg(forest_type)
  ci_method       <- match.arg(ci_method)
  zero_handling   <- match.arg(zero_handling)
  na_pred_handling<- match.arg(na_pred_handling)
  alpha <- 1 - conf_level
  
  # Choose prediction column (kept consistent with your previous code)
  get_pred <- function(dt) {
    if (forest_type == "grf") dt$crte_hat_grf else dt$crte_hat_glm
  }
  
  # Compute MAPE on a data.frame/data.table subset with robust handling
  point_mape <- function(dt, return_meta = FALSE) {
    truth <- dt$CRTE
    pred  <- get_pred(dt)
    
    # Non-finite handling
    finite_mask <- is.finite(truth) & is.finite(pred)
    n_total <- length(truth)
    n_nonfinite_pred  <- sum(!is.finite(pred))
    n_nonfinite_truth <- sum(!is.finite(truth))
    
    if (na_pred_handling == "error" && (n_nonfinite_pred > 0 || n_nonfinite_truth > 0)) {
      stop("Non-finite values detected in truth or predictions.")
    }
    
    dt_used_truth <- truth[finite_mask]
    dt_used_pred  <- pred[finite_mask]
    n_used_after_finite <- length(dt_used_truth)
    
    # Zero-denominator handling
    zeros_truth <- dt_used_truth == 0
    zero_zero_pairs <- zeros_truth & (dt_used_pred == 0)
    
    if (zero_handling == "drop") {
      keep <- !zeros_truth
      dt_used_truth <- dt_used_truth[keep]
      dt_used_pred  <- dt_used_pred[keep]
      zeros_truth   <- zeros_truth[keep]
      zero_zero_pairs <- zero_zero_pairs[keep]
    }
    
    denom <- abs(dt_used_truth)
    
    if (zero_handling %in% c("safe", "epsilon")) {
      # Prevent divide-by-zero
      denom <- pmax(denom, eps)
    }
    
    # Absolute percentage errors
    ape <- abs(dt_used_truth - dt_used_pred) / denom
    
    if (zero_handling == "safe") {
      # Explicitly declare 0/0 as zero error (truth==0 & pred==0)
      ape[zero_zero_pairs] <- 0
    }
    
    est <- if (length(ape) > 0) mean(ape) else NA_real_
    
    if (!return_meta) return(est)
    
    list(
      estimate            = est,
      n_units_total       = n_total,
      n_units_used        = length(ape),
      n_nonfinite_pred    = n_nonfinite_pred,
      n_nonfinite_truth   = n_nonfinite_truth,
      n_zero_truth_total  = sum(truth == 0, na.rm = TRUE),
      n_zero_zero_pairs   = sum(zero_zero_pairs, na.rm = TRUE),
      zero_handling       = zero_handling,
      eps                 = eps,
      na_pred_handling    = na_pred_handling
    )
  }
  
  # ---- Trial-specific path: unit-level bootstrap CI ----
  if (!is.null(trial)) {
    dt <- bundle[bundle$trial_id == trial, ]
    meta <- point_mape(dt, return_meta = TRUE)
    est  <- meta$estimate
    
    # Build a "used rows" view for bootstrap (respecting drops)
    # We reconstruct the used rows by filtering dt the same way as point_mape()
    truth <- dt$CRTE
    pred  <- get_pred(dt)
    finite_mask <- is.finite(truth) & is.finite(pred)
    dt_used <- dt[finite_mask, ]
    
    if (zero_handling == "drop") {
      dt_used <- dt_used[dt_used$CRTE != 0, ]
    }
    n_used <- nrow(dt_used)
    
    if (is.na(est) || is.null(n_used) || n_used < 2) {
      return(c(meta,
               list(ci_lower = NA_real_, ci_upper = NA_real_,
                    n_trials = 1L, conf_level = conf_level, ci_method = "unit-bootstrap")))
    }
    
    boot_means <- replicate(B, {
      idx <- sample.int(n_used, size = n_used, replace = TRUE)
      point_mape(dt_used[idx, ], return_meta = FALSE)
    })
    qs <- quantile(boot_means, probs = c(alpha/2, 1 - alpha/2), names = FALSE)
    
    return(c(meta,
             list(ci_lower = qs[1], ci_upper = qs[2],
                  n_trials = 1L, conf_level = conf_level, ci_method = "unit-bootstrap")))
  }
  
  # ---- Aggregate across trials path ----
  trials <- sort(unique(bundle$trial_id))
  
  # Point estimate: robust mean over all rows (non-finite dropped, zero handling applied)
  meta_all <- point_mape(bundle, return_meta = TRUE)
  est <- meta_all$estimate
  
  # If multiple trials, CI over trial-level MAPEs; else fall back to unit bootstrap
  if (length(trials) >= 2) {
    m_by_trial <- vapply(trials, function(tr) {
      dt_tr <- bundle[bundle$trial_id == tr, ]
      point_mape(dt_tr, return_meta = FALSE)
    }, numeric(1))
    
    # Drop trials where MAPE is NA (e.g., no usable rows)
    m_by_trial <- m_by_trial[is.finite(m_by_trial)]
    nT <- length(m_by_trial)
    
    if (nT >= 2 && ci_method == "t") {
      s <- sd(m_by_trial)
      tcrit <- qt(1 - alpha/2, df = nT - 1)
      halfw <- tcrit * s / sqrt(nT)
      lower <- est - halfw
      upper <- est + halfw
      method <- "t"
    } else if (nT >= 2 && ci_method == "bootstrap") {
      boot_means <- replicate(B, mean(sample(m_by_trial, size = nT, replace = TRUE)))
      qs <- quantile(boot_means, probs = c(alpha/2, 1 - alpha/2), names = FALSE)
      lower <- qs[1]; upper <- qs[2]; method <- "bootstrap"
    } else {
      # Not enough valid trials → unit bootstrap on all usable rows
      dt <- bundle
      truth <- dt$CRTE
      pred  <- get_pred(dt)
      finite_mask <- is.finite(truth) & is.finite(pred)
      dt_used <- dt[finite_mask, ]
      if (zero_handling == "drop") dt_used <- dt_used[dt_used$CRTE != 0, ]
      n_used <- nrow(dt_used)
      
      if (is.na(est) || is.null(n_used) || n_used < 2) {
        lower <- upper <- NA_real_
        method <- "unit-bootstrap"
      } else {
        boot_means <- replicate(B, {
          idx <- sample.int(n_used, size = n_used, replace = TRUE)
          point_mape(dt_used[idx, ], return_meta = FALSE)
        })
        qs <- quantile(boot_means, probs = c(alpha/2, 1 - alpha/2), names = FALSE)
        lower <- qs[1]; upper <- qs[2]; method <- "unit-bootstrap"
      }
    }
    
    return(c(meta_all,
             list(ci_lower = lower, ci_upper = upper,
                  n_trials = length(trials), n_trials_used = nT,
                  conf_level = conf_level, ci_method = method)))
  } else {
    # Single trial present → unit bootstrap CI
    dt <- bundle
    truth <- dt$CRTE
    pred  <- get_pred(dt)
    finite_mask <- is.finite(truth) & is.finite(pred)
    dt_used <- dt[finite_mask, ]
    if (zero_handling == "drop") dt_used <- dt_used[dt_used$CRTE != 0, ]
    n_used <- nrow(dt_used)
    
    if (is.na(est) || is.null(n_used) || n_used < 2) {
      return(c(meta_all,
               list(ci_lower = NA_real_, ci_upper = NA_real_,
                    n_trials = 1L, conf_level = conf_level, ci_method = "unit-bootstrap")))
    }
    
    boot_means <- replicate(B, {
      idx <- sample.int(n_used, size = n_used, replace = TRUE)
      point_mape(dt_used[idx, ], return_meta = FALSE)
    })
    qs <- quantile(boot_means, probs = c(alpha/2, 1 - alpha/2), names = FALSE)
    
    return(c(meta_all,
             list(ci_lower = qs[1], ci_upper = qs[2],
                  n_trials = 1L, conf_level = conf_level, ci_method = "unit-bootstrap")))
  }
}

# Backward-compatible wrapper: still returns only the point estimate (numeric)
mape_from_bundle <- function(bundle, forest_type = c("grf", "rrcf"), trial = NULL) {
  mape_from_bundle_detail(bundle, forest_type = forest_type, trial = trial)$estimate
}

# Detailed power with CIs for Athey-style omnibus test
athey_power_from_bundle_detail <- function(bundle,
                                           forest_type = c("grf", "rrcf"),
                                           alpha_test = 0.05,
                                           conf_level = 0.95,
                                           ci_method = c("wilson", "jeffreys", "exact")) {
  forest_type <- match.arg(forest_type)
  ci_method   <- match.arg(ci_method)
  alpha <- 1 - conf_level
  
  rd_pval <- function(base_df, tau_hat) {
    mask <- is.finite(tau_hat) & stats::complete.cases(base_df)
    if (sum(mask) < 10) return(NA_real_)
    Wm <- base_df$W[mask]
    tab <- table(Wm)
    if (length(tab) < 2 || any(tab < 5)) return(NA_real_)  # need both arms, min 5 each
    
    tau_c <- tau_hat[mask] - mean(tau_hat[mask])
    int   <- Wm * tau_c
    if (stats::sd(int) == 0) return(NA_real_)              # no variation ⇒ no test
    
    df_b <- base_df[mask, , drop = FALSE]
    df   <- cbind(df_b, int = int)
    
    base <- try(stats::lm(y ~ ., data = df_b), silent = TRUE)
    add  <- try(stats::lm(y ~ . + int, data = df), silent = TRUE)
    if (inherits(base, "try-error") || inherits(add, "try-error")) return(NA_real_)
    
    # If adding 'int' didn’t increase rank, ANOVA p-value is meaningless
    if (ncol(stats::model.matrix(add)) == ncol(stats::model.matrix(base))) return(NA_real_)
    
    av <- stats::anova(base, add)
    p  <- suppressWarnings(av$`Pr(>F)`[2])
    if (!is.finite(p)) NA_real_ else p
  }
  
  cov_cols  <- get_covariate_cols(bundle)
  trial_ids <- unique(bundle$trial_id)
  
  pvals <- rep(NA_real_, length(trial_ids))
  for (i in seq_along(trial_ids)) {
    sub <- bundle[trial_id == trial_ids[i]]
    base_df <- data.frame(y = sub$Y, W = sub$A, sub[, ..cov_cols])
    cate_hat <- if (forest_type == "grf") sub$cate_hat_grf else sub$cate_hat_glm
    pvals[i] <- rd_pval(base_df, cate_hat)
  }
  
  n_total <- length(pvals)
  finite_mask <- is.finite(pvals)
  n_used  <- sum(finite_mask)
  n_na    <- n_total - n_used
  
  if (n_used == 0) {
    return(list(
      estimate    = NA_real_,
      ci_lower    = NA_real_,
      ci_upper    = NA_real_,
      n_trials    = n_total,
      n_trials_used = n_used,
      k_reject    = NA_integer_,
      alpha_test  = alpha_test,
      conf_level  = conf_level,
      ci_method   = ci_method,
      note        = "No valid p-values (all NA/invalid)."
    ))
  }
  
  k <- sum(pvals[finite_mask] < alpha_test)
  p_hat <- k / n_used
  
  # Binomial proportion CI
  if (ci_method == "wilson") {
    z <- stats::qnorm(1 - alpha/2)
    denom  <- 1 + (z^2)/n_used
    center <- (p_hat + (z^2)/(2*n_used)) / denom
    halfw  <- (z * sqrt(p_hat*(1 - p_hat)/n_used + (z^2)/(4*n_used^2))) / denom
    lower  <- max(0, center - halfw)
    upper  <- min(1, center + halfw)
  } else if (ci_method == "jeffreys") {
    lower <- stats::qbeta(alpha/2,     k + 0.5, n_used - k + 0.5)
    upper <- stats::qbeta(1 - alpha/2, k + 0.5, n_used - k + 0.5)
  } else { # exact (Clopper–Pearson)
    lower <- if (k == 0) 0 else stats::qbeta(alpha/2,     k,     n_used - k + 1)
    upper <- if (k == n_used) 1 else stats::qbeta(1 - alpha/2, k + 1, n_used)
  }
  
  list(
    estimate      = p_hat,
    ci_lower      = lower,
    ci_upper      = upper,
    n_trials      = n_total,
    n_trials_used = n_used,
    k_reject      = k,
    alpha_test    = alpha_test,
    conf_level    = conf_level,
    ci_method     = ci_method,
    n_trials_na   = n_na
  )
}

# Original wrapper: stays backward-compatible (returns numeric point estimate)
athey_power_from_bundle <- function(bundle, forest_type = c("grf", "rrcf")) {
  athey_power_from_bundle_detail(bundle, forest_type = forest_type)$estimate
}

# Small helper for NULL-coalescing
`%||%` <- function(a, b) if (!is.null(a)) a else b

