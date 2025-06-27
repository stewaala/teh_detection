n = 100
K = 1.5

forms1  <- c(X1 ~ 1, X2 ~ 1)
fam1    <- c(1, 1)
link1   <- c("identity", "identity")

forms2  <- list(A ~ 1, Z ~ 1)
# fam2    <- c(5, 5)
# link2   <- c("logit", "logit")
fam2    <- c(5, 1)
link2   <- c("logit", "identity")

forms3  <- Y ~ A + Z + A:Z
fam3    <- c(1)
link3   <- "identity"

forms4 <- ~ Z
fam4 <- c(1, 1)

forms   <- list(forms1, forms2, forms3, forms4)
fam     <- list(fam1,   fam2,   fam3, fam4)
link    <- list(link1,  link2,  link3)

pars <- list(
  X1 = list(beta = 0,  phi = 1),
  X2 = list(beta = 0,  phi = 1),
  # Z  = list(beta = 0),
  Z = list(beta = 0, phi = 1),
  A  = list(beta = 0),
  Y  = list(beta = c(81.0, 0.2 , 3.5, 2 * K), phi = 1),
  cop = list(Y=list(X1=list(beta=c(0,0.5)), X2=list(beta=c(-0.5, 0))))
)

set.seed(42)
dat <- as.data.table(rfrugalParam(n, formulas = forms, pars = pars, family = fam, link = link))

dat2 <- copy(dat)
dat2[, A := 1]           
M1 <- model.matrix(forms[[3]], dat2)

dat2[, A := 0]           
M0 <- model.matrix(forms[[3]], dat2)

lp1        <- M1 %*% pars$Y$beta
lp0        <- M0 %*% pars$Y$beta

dat$Y1hat <- lp1
dat$Y0hat <- lp0
dat$CATE <- dat$Y1hat - dat$Y0hat


