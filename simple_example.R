library(causalForest)
library(randomForestCI)

n = 500
ntree = 1000
sigma = 1
d = 10

n.test = 1000

#
# This example has no treatment effect, but has a treatment propensity
# that is correlated with outcomes
#

baseline = function(x) {
	2 * (x[1] - 0.5)
}

propensity = function(x) {
	0.25 + dbeta(x[1], 2, 4)/4
}

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
se.hat = sqrt(forest.ci$var.hat)

up.lim = predictions + 1.96 * se.hat
down.lim = predictions - 1.96 * se.hat

# true tau is 0 here
coverage = mean(up.lim > 0 & down.lim < 0)
print(coverage)
