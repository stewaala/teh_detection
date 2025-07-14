## ------------------------------------------------------------------
##  SET‑UP -----------------------------------------------------------
## ------------------------------------------------------------------
library(data.table)

## ---- paste the copula_rct_identity() function here or source() it --
## (for brevity, assume it is already in your workspace)

## User‑tweakable settings ------------------------------------------
n_test <- 5e4       # sample size for the diagnostic run
rho    <- 0.5       # effect‑heterogeneity level to test
seed   <- 42        # reproducibility

## ------------------------------------------------------------------
##  1.  Generate a sample -------------------------------------------
## ------------------------------------------------------------------
set.seed(seed)
dat <- copula_rct_identity(n = n_test, rho = rho, seed = seed)

## Keep a copy of the mean 'p' for later checks
dat[, p_theory :=
      0.40 + 0.15*C1 + 0.20*sin(C4) +
      A * (-0.10 + 0.10*rho*(C1 + C2 + (C3 > 0) + C4^2)) ]

## ------------------------------------------------------------------
##  2.  Basic sanity checks -----------------------------------------
## ------------------------------------------------------------------
cat("----- BASIC SANITY CHECKS ------------------------------------\n")

## (a) Y is 0/1 and matches Bernoulli(p)
cat("Y unique values: ", paste(unique(dat$Y), collapse = ", "), "\n")

## (b) p is inside [0,1]
cat(sprintf("Range of theoretical p: [%.4f , %.4f]\n",
            min(dat$p_theory), max(dat$p_theory)))

## (c) Empirical mean of Y vs mean(p)
cat(sprintf("Mean(Y) = %.4f    Mean(p) = %.4f\n",
            mean(dat$Y), mean(dat$p_theory)))

## ------------------------------------------------------------------
##  3.  Check identity‑link *mean model* ----------------------------
## ------------------------------------------------------------------
cat("\n----- LINEAR‑PROBABILITY REGRESSION (should recover coefficients)\n")

## Fit the linear‑probability model the DGP was defined with
fit <- lm(Y ~ A + C1 + sin(C4) +
            A:C1 + A:C2 + A:I(C3>0) + A:I(C4^2),
          data = dat)

## True coefficient vector in the same order
true_beta <- c("(Intercept)" =  0.40,
               "A"           = -0.10,
               "C1"          =  0.15,
               "sin(C4)"     =  0.20,
               "A:C1"        =  0.10*rho,
               "A:C2"        =  0.10*rho,
               "A:I(C3 > 0)" =  0.10*rho,
               "A:I(C4^2)"   =  0.10*rho)

## Print table of estimated vs. true
tab <- cbind(Estimate = coef(fit)[names(true_beta)],
             Truth    = true_beta)
print(round(tab, 3))

## ------------------------------------------------------------------
##  4.  Additive treatment effect check -----------------------------
## ------------------------------------------------------------------
cat("\n----- AVERAGE TREATMENT EFFECT (ATE) -------------------------\n")

## Empirical ATE                     E[Y|A=1] - E[Y|A=0]
emp_ate <- dat[A==1, mean(Y)] - dat[A==0, mean(Y)]

## Theoretical ATE = E[δ(C)]         = E[ -0.10 + 0.10 ρ f(C) ]
dat[, delta_theory := -0.10 + 0.10*rho*
      (C1 + C2 + (C3 > 0) + C4^2) ]
theo_ate <- mean(dat$delta_theory)

cat(sprintf("Empirical ATE  = %.4f\n", emp_ate))
cat(sprintf("Theoretical ATE = %.4f\n", theo_ate))

## ------------------------------------------------------------------
##  5.  Copula‑induced rank correlations ----------------------------
## ------------------------------------------------------------------
cat("\n----- SPEARMAN CORRELATIONS  Y  vs  (X1,X2,X3) ---------------\n")

corr <- sapply(c("X1","X2","X3"),
               function(v) cor(dat$Y, dat[[v]], method = "spearman"))

print(round(corr, 3))

cat("\nExpected signs from copula loadings:  X1  +,   X2  –,   X3  +\n")

## ------------------------------------------------------------------
##  6.  Residual diagnostics (optional) -----------------------------
## ------------------------------------------------------------------
cat("\n----- OPTIONAL: RANDOMNESS OF RESIDUALS ----------------------\n")
## Pearson residuals should be uncorrelated with p_theory
resid <- dat$Y - dat$p_theory
cat(sprintf("Cor(residual, p_theory) = %.4f\n",
            cor(resid, dat$p_theory)))

## ------------------------------------------------------------------
cat("\n===== SUMMARY =====\n")
cat("• Coefficients should be within ~0.01 of truth with n ≥ 5e4.\n")
cat("• Empirical ATE should match theoretical ATE to ±0.005.\n")
cat("• Spearman correlations should have signs (+, –, +) and magnitude\n")
cat("  about 0.05–0.08 at ρ = 0.5 (monotone but purposely modest).\n")
