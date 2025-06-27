n = 10
K = 3.5
seed = 1
p0 = 0.2

forms1  <- c(X1 ~ 1, X2 ~ 1)
fam1    <- c(1, 1)
link1   <- c("identity", "identity")

forms2  <- list(A ~ 1, Z ~ 1)
# fam2    <- c(5, 5)
# link2   <- c("logit", "logit")
fam2    <- c(5, 1)
link2   <- c("logit", "identity")

forms3  <- Y ~ A + Z + A:Z
fam3    <- c(5)
link3   <- "log"

forms4 <- ~ Z
fam4 <- c(1, 1)

forms   <- list(forms1, forms2, forms3, forms4)
fam     <- list(fam1,   fam2,   fam3, fam4)
link    <- list(link1,  link2,  link3)

pars <- list(
  X1 = list(beta = 0,  phi = 1),
  X2 = list(beta = 0,  phi = 1),
  # Z  = list(beta = 0),
  Z  = list(beta = 0, phi = 1),
  A  = list(beta = 0),
  Y  = list(beta = c(log(p0), -log(K) , 0, 2 * log(K))),
  cop = list(Y=list(X1=list(beta=c(0,0.5)), X2=list(beta=c(-0.5, 0))))
)

set.seed(seed)
dat <- as.data.table(rfrugalParam(n, formulas = forms, pars = pars, family = fam, link = link))

dat2 <- copy(dat)
dat2[, A := 1]           
M1 <- model.matrix(forms[[3]], dat2)

dat2[, A := 0]           
M0 <- model.matrix(forms[[3]], dat2)

lp1        <- M1 %*% pars$Y$beta
lp0        <- M0 %*% pars$Y$beta

dat$Y1hat <- exp(lp1)
dat$Y0hat <- exp(lp0)
dat$CRTE <- dat$Y1hat / dat$Y0hat


