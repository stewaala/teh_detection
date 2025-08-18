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

poisson_power_from_bundle <- function(bundle, forest_type = c("grf", "rrcf"), small_sample=FALSE) {

  forest_type <- match.arg(forest_type)
  cov_cols <- get_covariate_cols(bundle)
  trials <- unique(bundle$trial_id)

  get_pval <- function(trial) {
    dt = bundle[bundle$trial_id == trial]
    x.test <- dt[, cov_cols, with = FALSE]
    pred <- if (forest_type == "grf") dt$crte_hat_grf else dt$crte_hat_glm
    poisson_omnibus_pval(pred, dt$Y, dt$A, x.test, small_sample)
  }
  pvals <- vapply(trials, FUN=get_pval, FUN.VALUE=numeric(1))
  mean(pvals < 0.05)
}

mape_from_bundle <- function(bundle, forest_type = c("grf", "rrcf"), trial = NULL) {

  if (!is.null(trial))
    bundle <- bundle[bundle$trial_id == trial,]
  
  mape_from_bundle_ci(bundle, forest_type)$mape

}

mape_from_bundle_ci <- function(bundle, forest_type = c("grf", "rrcf")) {
  
  pred <- if (forest_type == "grf") bundle$crte_hat_grf else bundle$crte_hat_glm
  ape <- abs((bundle$CRTE - pred)/bundle$CRTE)
  q <- quantile(ape[is.finite(ape)], c(0.05, 0.95), na.rm=TRUE)
  list('mape'=mean(ape[is.finite(ape)], na.rm=TRUE), 'quantile_5'=q[1], 'quantile_95'=q[2])

}

mae_from_bundle_ci <- function(bundle, forest_type = c("grf", "rrcf")) {
  
  pred <- if (forest_type == "grf") bundle$cate_hat_grf else bundle$cate_hat_glm
  ae <- abs(bundle$CATE - pred)
  q <- quantile(ae, c(0.05, 0.95))
  list('mae'=mean(ae), 'quantile_5'=q[1], 'quantile_95'=q[2])
}

mae_from_bundle <- function(bundle, forest_type = c("grf", "rrcf")) {
  
  mae_from_bundle_ci(bundle, forest_type)$mae
  
}


