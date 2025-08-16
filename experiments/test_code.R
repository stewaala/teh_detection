######################################################################################################
# This file shows that my functions create_bundle, poisson_omnibus_pval and my_mape_from_bundle
# reproduce the output from Shirvaikar's SimRCT.R file to machine precision
# Note: I have added the line set.seed(trial) before sample in Shirvaikar's code to make it reproducible
# It validates pred.grf, pred.glm, mape.grf, mape.glm, p.grf and p.glm
######################################################################################################

n=500
rho=0.75
trial=5

######################################################################################################
# Shirvaikar's SimRCT.R
######################################################################################################
library(causl)
library(data.table)
library(grf)
library(rrcf)

source("copula_rct.R")

data = copula_rct(n, rho, trial)
set.seed(trial)
split <- sample(seq_len(n), size = n*0.8)
x = data[, !c("A", "Y", "CRTE"), with=FALSE]
x.train <- x[split, ]
x.test <- x[-split, ]
y.train = data$Y[split]
y.test = data$Y[-split]
t.train = data$A[split]
t.test = data$A[-split]
crte.test = data$CRTE[-split]

forest.grf <- causal_forest(x.train, y.train, t.train, W.hat=0.5, num.trees=500, seed=1234)
forest.glm = rr_causal_forest(x.train, y.train, t.train, rct=TRUE, num.trees=500, seed=1234)
pred.grf = rr_predict(forest.grf, x.test)
pred.glm = rr_predict(forest.glm, x.test)
mape.grf = mean(abs((crte.test-pred.grf)/crte.test))
mape.glm = mean(abs((crte.test-pred.glm)/crte.test))

anova.data = data.frame(cbind(y.test, t.test, x.test))
model.base = glm(y.test ~ ., family = poisson, data = anova.data)
model.grf = glm(y.test ~ ., family = poisson, data = cbind(anova.data, t.test*log(pred.grf)))
model.glm = glm(y.test ~ ., family = poisson, data = cbind(anova.data, t.test*log(pred.glm)))
anova.grf = anova(model.base, model.grf)
anova.glm = anova(model.base, model.glm)
p.grf = 1 - pchisq(anova.grf$Deviance[2], df = 1)
p.glm = 1 - pchisq(anova.glm$Deviance[2], df = 1)

######################################################################################################
# reproduce the predicted risk ratios, p-values and mape variables and validate them against SimRCT.R
######################################################################################################
source("simulate_and_fit_forests.R")

bundle <- make_bundle(n, rho, link_type="log", num_trials=trial, num_trees=500, test_frac=0.2)
dt <- bundle[bundle$trial_id == trial,]
all.equal(pred.grf, dt$crte_hat_grf)
all.equal(pred.glm, dt$crte_hat_glm)

cov_cols <- get_covariate_cols(dt)
x.test <- dt[, cov_cols, with = FALSE]
my_p.grf <- poisson_omnibus_pval(dt$crte_hat_grf, dt$Y, dt$A, x.test)
my_p.glm <- poisson_omnibus_pval(dt$crte_hat_glm, dt$Y, dt$A, x.test)
all.equal(p.grf, my_p.grf)
all.equal(p.glm, my_p.glm)

my_mape.grf <- mape_from_bundle(bundle, forest_type="grf", trial=trial)
my_mape.rrcf <- mape_from_bundle(bundle, forest_type="rrcf", trial=trial)

all.equal(mape.grf, my_mape.grf)
all.equal(mape.glm, my_mape.rrcf)

poisson_power.grf <- poisson_power_from_bundle(bundle, "grf")
poisson_power.rrcf <- poisson_power_from_bundle(bundle, "rrcf")




