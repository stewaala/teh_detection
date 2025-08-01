library(causl)
library(data.table)

.EEXP_H_CONST <- 1.20332119584

alan_copula_rct <- function(n, rho = 0, seed = 111, p_base = .EEXP_H_CONST * exp(-2)){
  
  alpha <- log(p_base) - log(.EEXP_H_CONST)
  cat("alpha set to ", alpha, '\n')
  
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
  # so that alpha = log(p_base) - log(.EEXP_H_CONST) reproduces alpha = -2 by default.
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
  # ------------------------------------------------------------------------------
  # Validate output probabilities are within [0,1]
  # ------------------------------------------------------------------------------
  # After simulation, check that no probability exceeds 1, which can occur
  # under a log-link when η = α + h(C) + W·f(C) > 0, yielding μ = exp(η) > 1.
  # High p_base or strong heterogeneity may drive η positive. The 'causl'
  # library could silently clip these values back into [0,1], but we want to
  # trap any out-of-bounds predictions. We recompute μ_i = exp(η_i) for all
  # rows and stop() with an error if any μ_i > 1, ensuring invalid probabilities
  # are detected immediately.
  
  # 1. Rebuild the model matrix for Y’s linear predictor:
  Xy  <- model.matrix(forms[[3]], dat)
  # 2. Grab your Y‐beta vector (with the updated alpha):
  betay <- pars$Y$beta
  # 3. Compute eta and the implied mu = exp(eta):
  eta <- as.vector(Xy %*% betay)
  mu  <- exp(eta)
  
  # 4. Trap any invalid probabilities:
  if (any(mu > 1)) {
    first_bad <- which(mu > 1)[1]
    stop(sprintf(
      "Invalid probability: mu>1 on row %d (mu = %.3f).",
      first_bad, mu[first_bad]
    ))
  }
  #-----------------------------------------------------
  
  # construct a model matrix and calculate the true CATE
  dat2 <- copy(dat)
  mm <- model.matrix(forms[[3]],dat2[,A:=1])-model.matrix(forms[[3]],dat2[,A:=0]) 
  dat$CRTE <- exp(mm %*% pars$Y$beta)
  return(dat)
}