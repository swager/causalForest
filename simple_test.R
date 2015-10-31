library(causalForest)
n = 500
p = 10
X = matrix(rnorm(n * p), n, p)
Y = rnorm(n)
W = rbinom(n, 1, 0.5)

ct = causalTree(Y ~ X, treatment=W, split.option="CT", cv.option="TOT")
summary(ct)
min(predict(ct))
max(predict(ct))
