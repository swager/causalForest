rm(list = ls())

library(causalForest)

source("utils/knn.R")
source("utils/lasso.R")
source("utils/bart.R")
source("utils/causaltree.R")

n = 1000
d = 20
sigma = 1

n.test = 1000

baseline = function(x) {
	3 / (1 + exp(8 * (x[1] - 0.5))) - 1.5 +
	2 * (x[2] + x[3] + x[4] + x[5] - 3)  - x[1]^2 -
	x[1]^2
}

propensity = function(x) {
	0.25 + dbeta(x[1], 2, 4)/4
}

treat.effect = function(x) {
	0
}

X = matrix(runif(n * d, 0, 1), n, d) # features
e = apply(X, 1, propensity)
W = rbinom(n, 1, e) #treatment condition
Y = apply(X, 1, baseline) + W * apply(X, 1, treat.effect) + sigma * rnorm(n)

X.test = matrix(runif(n.test * d, 0, 1), n.test, d)
tau.test = apply(X.test, 1, treat.effect)

tau.knn = knn.cate(X, Y, W, X.test, k = 50)
mean(tau.knn)
mean((tau.knn - tau.test)^2)

tau.lasso = lasso.cate(X, Y, W, X.test)
mean(tau.lasso)
mean((tau.lasso - tau.test)^2)

tau.bart = bart.cate(X, Y, W, X.test)
mean(tau.bart)
mean((tau.bart - tau.test)^2)

cf.fit = causalForest(X, Y, W, 500, sample.size = round(n^0.8))
tau.cf = predict(cf.fit, newdata = data.frame(X=X.test))
mean(tau.cf)
mean((tau.cf - tau.test)^2)

pf.fit = propensityForest(X, Y, W, 500, sample.size = round(n^0.8))
tau.pf = predict(pf.fit, newdata = data.frame(X=X.test))
mean(tau.pf)
mean((tau.pf - tau.test)^2)

tau.ct = causaltree.cate(X, Y, W, X.test)
mean(tau.ct)
mean((tau.ct - tau.test)^2)