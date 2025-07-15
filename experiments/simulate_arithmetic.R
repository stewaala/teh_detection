simulate_arithmetic <- function(n, K, seed = NULL) {
  # n: sample size
  # K: controls size of Z×A interaction (true CATE = 2*K*Z)
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
  
  # 3. Outcome model Y | A, Z
  forms3 <- Y ~ A + Z + A:Z
  fam3   <- c(1)                           # 1 = Gaussian
  link3  <- "identity"
  
  # 4. (Optional) Additional copula args
  forms4 <- ~ Z
  fam4   <- c(1, 1)
  
  # Bundle them up
  forms <- list(forms1, forms2, forms3, forms4)
  fam   <- list(fam1,   fam2,   fam3,   fam4)
  link  <- list(link1,  link2,  link3)
  
  # 5. Parameter list
  pars <- list(
    X1  = list(beta = 0,      phi = 1),
    X2  = list(beta = 0,      phi = 1),
    Z   = list(beta = 0,      phi = 1),
    A   = list(beta = 0),                # intercept only → p = 0.5
    #Y   = list(beta = c(25., 5., 7., 2 * K), phi = 1),
    Y   = list(beta = c(100, 5, 7, 3. * K), phi = 1),
    cop = list(
      Y = list(
        X1 = list(beta = c( 0,  0.5)),
        X2 = list(beta = c(-0.5, 0  ))
      )
    )
  )
  
  # 6. Simulate via causl
  dat <- as.data.table(
    rfrugalParam(
      n         = n,
      formulas  = forms,
      pars      = pars,
      family    = fam,
      link      = link
    )
  )
  
  return(dat)
}
