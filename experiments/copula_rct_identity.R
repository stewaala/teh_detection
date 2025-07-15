# This is a modified version of copula_rct(), the function used to simulate data
# in the Shirvaikar paper. copula_rct() generates treatment effects on the log 
# scale. copula_rct_identity() generates treatment effects on the linear scale.
# 
# copula_rct() uses a log link with the Bernoulli outcome. The linear model has been 
# designed so that its value is always negative and hence leads to a valid expected 
# outcome in the range [0, 1]
#
# To generate treatment effects on the linear scale, copula_rct_identity() removes 
# log link and the linear model has been designed so that its value is always in
# range [0, 1] and hence leads to a valid expected outcome
#
# However, unfortunately, causl does not permit the use of an identity link with
# a Bernoulli variable and so the generation of the outcome variable is done in 
# two steps:
# 1) A latent Gaussian variable E with distribution N(0, 1), is defined via a copula
#    jointly with X1, X2 and X3. This maintains the correlation structure of the 
#.   outcome with these covariates. All data except the outcome Y is simulated.
# 2) Y is simulated by comparing the inverse cumulative of the latent variable E
#.   with the linear model valued at X1, X2 and X3

library(causl)
library(data.table)

# Distribution family codes:
# 1  ⇒  Gaussian (normal)
# 2  ⇒  Student-t
# 3  ⇒  Gamma
# 4  ⇒  Beta
# 5  ⇒  Bernoulli

copula_rct_identity <- function(n, rho = 0, seed = 111, clip = TRUE) {
  
  ## 0.  helper: original pair‑copula loadings ------------------------
  pc_load <- list(
    X1 = list(beta = c( 0,  0.5)),
    X2 = list(beta = c(-0.5, 0)),
    X3 = list(beta = c( 0.5, 0))
  )
  
  ## 1.  four‑block set‑up -------------------------------------------
  fam <- list(
    c(1, 3, 4),                 # block‑1  : X1, X2, X3
    c(5, 5, 5, 1, 2),           # block‑2  : A, C1–C4
    c(1),                       # block‑3  : latent error  E  (Gaussian)
    c(1, 3, 4)                  # block‑4  : tilt
  )
  
  forms <- list(
    ## prognostic block
    c(X1 ~ 1,
      X2 ~ X1,
      X3 ~ X1),
    
    ## treatment + effect modifiers
    list(A  ~ 1,
         C1 ~ 1,
         C2 ~ C1,
         C3 ~ C1:C2,
         C4 ~ C1),
    
    ## latent error block  (response = E)
    E ~ 1,
    
    ## dependence tilt (unchanged)
    ~ C2
  )
  
  pars <- list(
    ## --------------- marginals ----------------
    X1 = list(beta = 0,           phi = 1),
    X2 = list(beta = c(0.1, 0.2), phi = 1),
    X3 = list(beta = c(0.1, 0.1), phi = 1),
    
    A  = list(beta = c(0)),
    C1 = list(beta = c(0)),
    C2 = list(beta = c(-2, 1)),
    C3 = list(beta = c(0, 0.1),   phi = 1),
    C4 = list(beta = c(0, 0.1),   phi = 0.1, par2 = 20),
    
    ## latent Gaussian error
    E  = list(beta = 0, phi = 1),
    
    ## --------------- copula -------------------
    ## only *responses* may appear here → E is the sole entry
    cop = list(E = pc_load)
  )
  
  link <- list(
    c("identity", "log", "logit"),               # block‑1
    c("logit", "logit", "logit",
      "identity", "identity"),                   # block‑2
    "identity"                                   # block‑3  (E)
  )
  
  ## 2.  draw everything jointly --------------------------------------
  set.seed(seed)
  dat <- as.data.table(
    rfrugalParam(n,
                 formulas = forms,
                 pars     = pars,
                 family   = fam,
                 link     = link)
  )
  
  ## 3.  convert latent error to uniform & threshold ------------------
  if (clip)                                     # keep theoretical bound
    dat[, C4 := pmax(pmin(C4, 0.5), -0.5)]
  
  dat[, p := 0.40 + 0.15 * C1 + 0.20 * sin(C4) +
        A * (-0.10 + 0.10 * rho *
               (C1 + C2 + (C3 > 0) + C4^2))]
  
  dat[, p := pmin(pmax(p, 0), 1)]               # safety clamp
  
  ## Φ maps N(0,1) → Uniform(0,1)
  dat[, U := pnorm(E)]
  dat[, Y := as.integer(U < p)]
  
  ## final clean‑up
  dat[, c("p", "E", "U") := NULL]
  return(dat[])
}
