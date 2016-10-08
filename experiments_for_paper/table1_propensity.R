library(causalForest)
library(mgcv)
library(randomForestCI)
library(FNN)
library(Hmisc)
library(xtable)

rm(list = ls())

n = 500
ntree = 1000
sigma = 1

n.test = 100

dvals = c(2, 5, 10, 15, 20, 30, 40)
simu.reps = 500

baseline = function(x) {
	2 * (x[1] - 0.5)
}

propensity = function(x) {
	0.25 + dbeta(x[1], 2, 4)/4
}

simu.fun = function(d) {

X = matrix(runif(n * d, 0, 1), n, d) # features
e = apply(X, 1, propensity)
W = rbinom(n, 1, e) #treatment condition

# no treatment effect
Y = apply(X, 1, baseline) + sigma * rnorm(n)

X.test = matrix(runif(n.test * d, 0, 1), n.test, d)


#
# random forest
#

forest = propensityForest(X, Y, W, num.trees = ntree, sample.size = n / 10, nodesize = 1)
predictions = predict(forest, X.test)
forest.ci = randomForestInfJack(forest, X.test, calibrate = TRUE)

# true effect is 0
se.hat = sqrt(forest.ci$var.hat)
rf.covered = mean(abs(predictions) <= 1.96 * se.hat)
rf.mse = mean(predictions^2)

k.small = 10
knn.0.mu = knn.reg(X[W==0,], X.test, Y[W==0], k = k.small)$pred
knn.1.mu = knn.reg(X[W==1,], X.test, Y[W==1], k = k.small)$pred

knn.0.mu2 = knn.reg(X[W==0,], X.test, Y[W==0]^2, k = k.small)$pred
knn.1.mu2 = knn.reg(X[W==1,], X.test, Y[W==1]^2, k = k.small)$pred

knn.0.var = (knn.0.mu2 - knn.0.mu^2) / (k.small - 1)
knn.1.var = (knn.1.mu2 - knn.1.mu^2) / (k.small - 1)

knn.tau = knn.1.mu - knn.0.mu
knn.se = sqrt(knn.0.var + knn.1.var)
knn.covered.v1 = mean(abs(knn.tau) <= 1.96 * knn.se)
knn.covered.v2 = mean(abs(knn.tau) <= 1.96 * sqrt(mean(knn.0.var + knn.1.var)))
knn.mse = mean(knn.tau^2)

k.big = 100
knnbig.0.mu = knn.reg(X[W==0,], X.test, Y[W==0], k = k.big)$pred
knnbig.1.mu = knn.reg(X[W==1,], X.test, Y[W==1], k = k.big)$pred

knnbig.0.mu2 = knn.reg(X[W==0,], X.test, Y[W==0]^2, k = k.big)$pred
knnbig.1.mu2 = knn.reg(X[W==1,], X.test, Y[W==1]^2, k = k.big)$pred

knnbig.0.var = (knnbig.0.mu2 - knnbig.0.mu^2) / (k.big - 1)
knnbig.1.var = (knnbig.1.mu2 - knnbig.1.mu^2) / (k.big - 1)

knnbig.tau = knnbig.1.mu - knnbig.0.mu
knnbig.se = sqrt(knnbig.0.var + knnbig.1.var)
knnbig.covered.v1 = mean(abs(knnbig.tau) <= 1.96 * knnbig.se)
knnbig.covered.v2 = mean(abs(knnbig.tau) <= 1.96 * sqrt(mean(knnbig.0.var + knnbig.1.var)))
knnbig.mse = mean(knnbig.tau^2)


c(rf.covered = rf.covered,
           rf.mse = rf.mse,
           knn.covered.v1 = knn.covered.v1,
           knn.covered.v2 = knn.covered.v2,
           knn.mse = knn.mse,
           knnbig.covered.v1 = knnbig.covered.v1,
           knnbig.covered.v2 = knnbig.covered.v2,
           knnbig.mse = knnbig.mse)
}

results.raw = lapply(dvals, function(d) {
	print(paste("NOW RUNNING:", d))
	res.d = sapply(1:simu.reps, function(iter) simu.fun(d))
	res.fixed = data.frame(t(res.d))
	print(paste("RESULT AT", d, "IS", colMeans(res.fixed)))
	res.fixed
})

results.condensed = lapply(results.raw, function(RR) {
	RR.mu = colMeans(RR)
	RR.var = sapply(RR, var) / (nrow(RR) - 1)
	rbind("mu"=RR.mu, "se"=sqrt(RR.var))
})

results.condensed

save.image("table1_propensity.RData")

results.parsed = lapply(results.condensed, function(RR) {
	apply(RR, 2, function(arg) {
		paste0(round(arg[1], 2), " (", round(100 * arg[2], 0), ")")
	})
})

results.table = data.frame(cbind(d=dvals, Reduce(rbind, results.parsed)))

results.table = results.table[,c(1, 3, 6, 9, 2, 4, 7)]
xtab = xtable(results.table)
print(xtab, include.rownames = FALSE)
