library(causl)
library(data.table)

#’ -----------------------------------------------------------------------------
#’ copula_rct_sim: Simulate RCT data with heterogeneous treatment effects
#’
#’ @description
#’   Generates covariates, treatment, and a binary outcome Y under a vine‑copula DGP.
#’   You can choose either a multiplicative (“log link”) or additive (“identity link”)
#’   treatment‐effect model.  The function exposes two main “knobs”:
#’     • p_base  – Target marginal untreated risk P(Y=1 | A=0)
#’     • rho     – Scale of heterogeneity in δ(C) = f(C)*rho + offset
#’   Passing the same p_base and rho to both link types yields
#’   comparable data sets (same average baseline risk, same f(C) scaling).
#’
#’   Outputs include:
#’     • Y    – simulated 0/1 outcome
#’     • CATE – true conditional risk difference E[Y|A=1,C] − E[Y|A=0,C]
#’     • CRTE – true conditional risk ratio   E[Y|A=1,C] / E[Y|A=0,C]
#’
#’ @usage
#’   dat <- copula_rct_sim(
#’              n = 10000,
#’              rho = 0.5,
#’              p_base = 0.1,
#’              link_type = "log",
#’              seed = 123,
#’              clip_C4 = TRUE
#’          )
#’
#’ @param n         Integer. Number of observations to simulate.
#’ @param rho       Numeric. Heterogeneity scale: multiplies f(C) in δ(C).
#’ @param p_base    Numeric in [0,1]. Desired overall P(Y=1 | A=0).
#’ @param link_type Character. Either "log" (multiplicative effects) or
#’                  "identity" (additive effects).
#’ @param seed      Integer. RNG seed for reproducibility.
#’ @param clip_C4   Logical. If TRUE, C4 is clipped to [−0.5,0.5] to bound f(C).
#’
#’ @return
#’   A data.table with columns:
#’     X1, X2, X3, A, C1, C2, C3, C4, Y, CATE, CRTE
#’
#’ @details
#’   • Covariate & copula structure is identical across link types:
#’       – Block 1: X1∼Gaussian, X2∼Gamma, X3∼Beta
#’       – Block 2: A,C1,C2∼Bernoulli; C3∼Gaussian; C4∼Student‑t(20)
#’       – Block 3: latent Gaussian E (for Y–X copula links)
#’       – Block 4: tilt on C2
#’
#’   • Treatment A is randomized: logit P(A=1)=0 ⇒ P(A=1)=0.5
#’
#’   • Prognostic score h(C) is centered so its mean is zero, making p_base
#’     exactly the marginal untreated risk.
#’
#’   • δ(C) = 
#’       – log link:  −0.2 + rho * f(C)  (on log–RR scale)
#’       – identity: −0.10 + 0.10 * rho * f(C)  (on risk‑difference scale)
#’     where f(C) = C1 + C2 + I(C3>0) + C4².
#’
#’   • Outcome Y is drawn by thresholding the latent copula uniform:
#’       U = Φ(E) ∼ Uniform(0,1);
#’       Y = 1{ U < p(C,A) }.
#’
#’ @section Distribution codes:
#’   1 → Gaussian (normal)  
#’   2 → Student‑t       
#’   3 → Gamma         
#’   4 → Beta          
#’   5 → Bernoulli     
#’
#’ -----------------------------------------------------------------------------
copula_rct_sim <- function(n,
                           rho        = 0.50,   # heterogeneity scale
                           p_base     = 0.10,   # target untreated prevalence
                           link_type  = c("log", "identity"),
                           seed       = 111,
                           clip_C4    = TRUE) {
  
  link_type <- match.arg(link_type)
  
  ## 1. Generate covariates + latent Gaussian error -----------------
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
  
  ## 2. Optionally clip C4 so C4^2 ≤ 0.25 ---------------------------
  if (clip_C4)
    dat[, C4 := pmax(pmin(C4,0.5),-0.5)]
  
  ## 3. Prognostic score h(C) and centred version -------------------
  if (link_type == "log") {
    h  <- with(dat, 0.3*C1 + 0.4*sin(C4))
  } else {                     # identity
    h  <- with(dat, 0.15*C1 + 0.20*sin(C4))
  }
  h_cent <- h - mean(h)        # ensures E[h_cent]=0
  
  ## 4. Intercept α₀ to match marginal untreated prevalence ---------
  if (link_type == "log") {
    alpha0 <- log(p_base) - log(mean(exp(h_cent)))
  } else {                     # identity
    alpha0 <- p_base           # because mean(h_cent)=0
  }
  
  ## 5. Heterogeneous treatment effect δ(C) -------------------------
  fC    <- with(dat, C1 + C2 + (C3>0) + C4^2)
  delta <- if (link_type == "log") {
    -0.2 + rho * fC
  } else {
    -0.10 + 0.10 * rho * fC
  }
  
  ## 6. Untreated prob p0(C) & treated prob p1(C) -------------------
  p0 <- if (link_type == "log") exp(alpha0 + h_cent)
  else                    alpha0 + h_cent
  
  p1 <- if (link_type == "log") p0 * exp(delta)
  else                    p0 + delta
  
  ## Bounds for identity link (log variant already ≤1 by design)
  p0 <- pmin(pmax(p0,0),1)
  p1 <- pmin(pmax(p1,0),1)
  
  ## 7. Ground‑truth effects ----------------------------------------
  CATE <- p1 - p0
  CRTE <- ifelse(p0 > 0, p1 / p0, NA_real_)      # safe division
  
  ## 8. Draw Bernoulli outcome via copula latent uniform ------------
  dat[, Y := as.integer(pnorm(E) < ifelse(A==1, p1, p0))]
  
  ## 9. Attach ground‑truth columns & clean up ----------------------
  dat[, `:=`(
    CATE = CATE,
    CRTE = CRTE,
    p1 = p1,
    p0 = p0
  )]
  dat[, E := NULL]
  
  return(dat[])
}
