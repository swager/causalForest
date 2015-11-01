library(causalForest)
n = 500
p = 10
X = matrix(rnorm(n * p), n, p)
Y = rnorm(n)
W = rbinom(n, 1, 0.5)

ct = causalTree(Y ~ X, treatment=W, split.option="CT", cv.option="TOT")
#summary(ct)
min(predict(ct))
max(predict(ct))
class(ct) == "causalTree"

ct = causalTree(Y ~ X, treatment=W, split.option="TOT", cv.option="TOT")
#summary(ct)
min(predict(ct))
max(predict(ct))
class(ct) == "causalTree"

cf = causalForest(X, Y, W, num.trees = 10)
min(predict(cf, X))
max(predict(cf, X))
class(cf) == "causalForest"

pf = propensityForest(X, Y, W, num.trees = 10)
min(predict(pf, X))
max(predict(pf, X))
class(pf) == "causalForest"
