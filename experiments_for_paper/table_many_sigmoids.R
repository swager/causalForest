rm(list = ls())

library(causalForest)
library(mgcv)
library(randomForestCI)
library(FNN)
library(Hmisc)
library(xtable)

n = 2000
ntree = 4000
sigma = 1

k.small = 10
k.big = 100

n.test = 1000
simu.reps = 20

k.vals = c(2, 4, 6)
d.vals = c(6, 12)

params = expand.grid(k=k.vals, d=d.vals)

treat.effect = function(x, k) {
	4/k * sum(sapply(1:k, function(kk) {
		1/(1 + exp(-12 * (x[kk] - 0.5))) - 0.5
	}))
}

simu.fun = function(iter) {

d = params[iter, "d"]
k = params[iter, "k"]

X = matrix(runif(n * d, 0, 1), n, d) # features
W = rbinom(n, 1, 0.5) #treatment condition

# no treatment effect
Y = W * apply(X, 1, function(x) treat.effect(x, k)) + sigma * rnorm(n)

X.test = matrix(runif(n.test * d, 0, 1), n.test, d)
true.eff = apply(X.test, 1, function(x) treat.effect(x, k))

#
# random forest
#

forest = causalForest(X, Y, W, num.trees = ntree, sample.size = n / 4, nodesize = 1)
predictions = predict(forest, X.test)
forest.ci = randomForestInfJack(forest, X.test, calibrate = TRUE)

rf.tau = forest.ci$y.hat
rf.se = sqrt(forest.ci$var.hat)
rf.cov = abs(rf.tau - true.eff) <= 1.96 * rf.se
rf.covered = mean(rf.cov)
rf.mse = mean((rf.tau - true.eff)^2)

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

results.raw = lapply(1:nrow(params), function(iter) {
	print(paste("NOW RUNNING:", iter))
	res.d = sapply(1:simu.reps, function(rr) simu.fun(iter))
	res.fixed = data.frame(t(res.d))
	print(paste("RESULT AT", iter, "IS", colMeans(res.fixed)))
	res.fixed
})

results.condensed = lapply(results.raw, function(RR) {
	RR.mu = colMeans(RR)
	RR.var = sapply(RR, var) / (nrow(RR) - 1)
	rbind("mu"=RR.mu, "se"=sqrt(RR.var))
})

results.condensed

save.image("table4_many_sigmoids.RData")

results.parsed = lapply(results.condensed, function(RR) {
	apply(RR, 2, function(arg) {
		paste0(round(arg[1], 2), " (", round(100 * arg[2], 0), ")")
	})
})

results.table = data.frame(cbind(k=params[,"k"], d=params["d"], Reduce(rbind, results.parsed)))
results.table

results.table = results.table[,c(1, 2, 4, 6, 8, 3, 5, 7)]
xtab = xtable(results.table)
print(xtab, include.rownames = FALSE)

