library(causalForest)

rm(list = ls())

n = 10000
s = 2000
ntree = 10000
sigma = 1

n.test = 1000

dvals = c(2, 3, 4, 5, 6, 8)
simu.reps = 40

effect = function(x) {
	4 / ((1 + exp(-12 * (x[1] - 0.5))) * (1 + exp(-12 * (x[2] - 0.5)))) 
}

simu.fun = function(d) {

X = matrix(runif(n * d, 0, 1), n, d) # features
W = rbinom(n, 1, 0.5) #treatment condition
Y = (W - 0.5) * apply(X, 1, effect) + sigma * rnorm(n)

X.test = matrix(runif(n.test * d, 0, 1), n.test, d)
true.eff = apply(X.test, 1, effect)

forest = causalForest(X, Y, W, num.trees = ntree, sample.size = s/2, nodesize = 1)
predictions = predict(forest, X.test)
rf.mse = mean((predictions - true.eff)^2)
rf.mse
}

results.raw = lapply(dvals, function(d) {
	print(paste("NOW RUNNING:", d))
	res.d = sapply(1:simu.reps, function(iter) simu.fun(d))
	print(paste("RESULT AT", d, "IS", mean(res.d)))
	res.d
})

results.condensed = unlist(lapply(results.raw, mean))

results.condensed

save.image("table3_hard_sigmoid_adaptMSE.RData")

round(results.condensed, 3)