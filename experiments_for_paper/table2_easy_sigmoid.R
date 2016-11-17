library(causalForest)
library(mgcv)
library(randomForestCI)
library(FNN)
library(Hmisc)
library(xtable)

rm(list = ls())

n = 5000
s = 2500
ntree = 2000
sigma = 1

n.test = 1000

dvals = c(2, 3, 4, 5, 6, 8)
simu.reps = 40

effect = function(x) {
	(1 + 1 / (1 + exp(-20 * (x[1] - 1/3)))) * (1 + 1 / (1 + exp(-20 * (x[2] - 1/3))))
}

simu.fun = function(d) {

X = matrix(runif(n * d, 0, 1), n, d) # features
W = rbinom(n, 1, 0.5) #treatment condition
Y = (W - 0.5) * apply(X, 1, effect) + sigma * rnorm(n)

X.test = matrix(runif(n.test * d, 0, 1), n.test, d)
true.eff = apply(X.test, 1, effect)

#
# random forest
#

forest = causalForest(X, Y, W, num.trees = ntree, sample.size = s/2, nodesize = 1)
predictions = predict(forest, X.test)
forest.ci = randomForestInfJack(forest, X.test, calibrate = TRUE)

se.hat = sqrt(forest.ci$var.hat)
rf.cov = abs(predictions - true.eff) <= 1.96 * se.hat
rf.covered = mean(rf.cov)
rf.mse = mean((predictions - true.eff)^2)

k.small = 7
knn.0.mu = knn.reg(X[W==0,], X.test, Y[W==0], k = k.small)$pred
knn.1.mu = knn.reg(X[W==1,], X.test, Y[W==1], k = k.small)$pred

knn.0.mu2 = knn.reg(X[W==0,], X.test, Y[W==0]^2, k = k.small)$pred
knn.1.mu2 = knn.reg(X[W==1,], X.test, Y[W==1]^2, k = k.small)$pred

knn.0.var = (knn.0.mu2 - knn.0.mu^2) / (k.small - 1)
knn.1.var = (knn.1.mu2 - knn.1.mu^2) / (k.small - 1)

knn.tau = knn.1.mu - knn.0.mu
knn.se = sqrt(knn.0.var + knn.1.var)
knn.cov = abs(knn.tau - true.eff) <= 1.96 * knn.se
knn.covered = mean(knn.cov)
knn.mse = mean((knn.tau - true.eff)^2)

k.big = 50
knnbig.0.mu = knn.reg(X[W==0,], X.test, Y[W==0], k = k.big)$pred
knnbig.1.mu = knn.reg(X[W==1,], X.test, Y[W==1], k = k.big)$pred

knnbig.0.mu2 = knn.reg(X[W==0,], X.test, Y[W==0]^2, k = k.big)$pred
knnbig.1.mu2 = knn.reg(X[W==1,], X.test, Y[W==1]^2, k = k.big)$pred

knnbig.0.var = (knnbig.0.mu2 - knnbig.0.mu^2) / (k.big - 1)
knnbig.1.var = (knnbig.1.mu2 - knnbig.1.mu^2) / (k.big - 1)

knnbig.tau = knnbig.1.mu - knnbig.0.mu
knnbig.se = sqrt(knnbig.0.var + knnbig.1.var)
knnbig.cov = abs(knnbig.tau - true.eff) <= 1.96 * knnbig.se
knnbig.covered = mean(knnbig.cov)
knnbig.mse = mean((knnbig.tau - true.eff)^2)


c(rf.covered = rf.covered,
           rf.mse = rf.mse,
	       knn.covered = knn.covered, 
           knn.mse = knn.mse,
	       knnbig.covered = knnbig.covered, 
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

save.image("table2_easy_sigmoid.RData")

results.parsed = lapply(results.condensed, function(RR) {
	apply(RR, 2, function(arg) {
		paste0(round(arg[1], 2), " (", round(100 * arg[2], 0), ")")
	})
})

results.table = data.frame(cbind(d=dvals, Reduce(rbind, results.parsed)))
results.table

results.table = results.table[,c(1, 3, 5, 7, 2, 4, 6)]
xtab = xtable(results.table)
print(xtab, include.rownames = FALSE)
