library(causalForest)

rm(list = ls())

n = 5000
s = 2500
ntree = 500
sigma = 1

n.test = 1000

d = 8
simu.reps = 100

leaf = 2^(0:9)

effect = function(x) {
	(1 + 1 / (1 + exp(-20 * (x[1] - 1/3)))) * (1 + 1 / (1 + exp(-20 * (x[2] - 1/3))))
}

simu.fun = function(minsz) {

X = matrix(runif(n * d, 0, 1), n, d) # features
W = rbinom(n, 1, 0.5) #treatment condition
Y = (W - 0.5) * apply(X, 1, effect) + sigma * rnorm(n)

X.test = matrix(runif(n.test * d, 0, 1), n.test, d)
true.eff = apply(X.test, 1, effect)

#
# random forest
#

forest = causalForest(X, Y, W, num.trees = ntree, sample.size = s/2, nodesize = minsz, honest = FALSE)
predictions = predict(forest, X.test)
rf.mse = mean((predictions - true.eff)^2)

forest2 = causalForest(X, Y, W, num.trees = ntree, sample.size = s/2, nodesize = minsz, honest = TRUE)
predictions2 = predict(forest2, X.test)
rf.mse2 = mean((predictions2 - true.eff)^2)

c(rf.mse, rf.mse2)
}

results.raw = lapply(leaf, function(minsz) {
	print(paste("NOW RUNNING:", minsz))
	res.d = sapply(1:simu.reps, function(iter) simu.fun(minsz))
	print(paste("RESULT AT", minsz, "IS", rowMeans(res.d)))
	res.d
})

save.image("table4_honest_vs_adapt.RData")

results.condensed = Reduce(rbind, lapply(results.raw, rowMeans))

results.condensed

round(results.condensed, 3)
