simulate_geometric <- function(n, K, seed = NULL) {
  # n:   sample size
  # K:   controls size of Z×A interaction (true CRTE = exp(2*log(K)*Z) = K^(2Z))
  # p0:  baseline success probability when A=0, Z=0
  # seed: optional RNG seed for reproducibility
  
  if (!is.null(seed)) set.seed(seed)
  
  # 1. Marginals for X1 and X2
  forms1 <- c(X1 ~ 1, X2 ~ 1)
  fam1   <- c(1, 1)                        # 1 = Gaussian
  link1  <- c("identity", "identity")
  
  # 2. Marginals for A (treatment) and Z (moderator)
  forms2 <- list(A ~ 1, Z ~ 1)
  fam2   <- c(5, 1)                        # 5 = Bernoulli, 1 = Gaussian
  link2  <- c("logit", "identity")
  
  # 3. Outcome model Y | A, Z: positive real via Gamma(log)
  forms3 <- Y ~ A + Z + A:Z
  fam3   <- c(3)                           # 3 = Gamma
  link3  <- "log"
  
  # 4. Copula args (as before)
  forms4 <- ~ Z
  fam4   <- c(1, 1)
  
  # Bundle
  forms <- list(forms1, forms2, forms3, forms4)
  fam   <- list(fam1,   fam2,   fam3,   fam4)
  link  <- list(link1,  link2,  link3)
  
  # 5. Parameter lists
  pars <- list(
    X1  = list(beta = 0,      phi = 1),
    X2  = list(beta = 0,      phi = 1),
    Z   = list(beta = 0,      phi = 1),
    A   = list(beta = 0),                # intercept only → p=0.5
    Y   = list(
      #beta = c(log(p0), -log(K), 0, 2 * log(K)),
      beta = c(log(0.2), -log(K), 2., 2 * log(K)),
      phi  = 1                            # dispersion for Gamma
    ),
    cop = list(
      Y = list(
        X1 = list(beta = c( 0,  0.5)),
        X2 = list(beta = c(-0.5, 0  ))
      )
    )
  )
  
  # 6. Simulate
  dat <- as.data.table(
    rfrugalParam(
      n        = n,
      formulas = forms,
      pars     = pars,
      family   = fam,
      link     = link
    )
  )
  
  # # 7. True potential outcomes & CRTE
  # dat2 <- copy(dat)
  # dat2[, A := 1]
  # M1 <- model.matrix(forms3, dat2)
  # 
  # dat2[, A := 0]
  # M0 <- model.matrix(forms3, dat2)
  # 
  # lp1 <- M1 %*% pars$Y$beta
  # lp0 <- M0 %*% pars$Y$beta
  # 
  # dat[, Y1hat := as.numeric(exp(lp1))]
  # dat[, Y0hat := as.numeric(exp(lp0))]
  # dat[, CRTE   := Y1hat / Y0hat]
  
  return(dat)
}
