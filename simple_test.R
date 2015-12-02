library(causalForest)
n = 500
p = 10

X = matrix(rnorm(n * p), n, p)
colnames(X) = sample.int(p)
Y = rnorm(n)
DF = data.frame(X=X, Y=Y)

W = rbinom(n, 1, 0.5)

ct = causalTree(Y ~ ., data = DF, treatment=W, split.option="TOT", cv.option="TOT")
#summary(ct)
min(predict(ct))
max(predict(ct))
class(ct) == "causalTree"

rct = refit.causalTree(ct, newdata=DF, treatment=W, propensity=rep(0.5, n))
min(predict(rct))
max(predict(rct))
class(rct) == "causalTree"

ct = causalTree(Y ~ ., data = DF, treatment=W, split.option="CT", cv.option="TOT")
#summary(ct)
min(predict(ct))
max(predict(ct))
class(ct) == "causalTree"

rct = refit.causalTree(ct, newdata=DF, treatment=W)
min(predict(rct))
max(predict(rct))
class(rct) == "causalTree"

cf = causalForest(X, Y, W, num.trees = 10)
min(predict(cf, X))
max(predict(cf, X))
class(cf) == "causalForest"

pf = propensityForest(X, Y, W, num.trees = 10)
min(predict(pf, X))
max(predict(pf, X))
class(pf) == "causalForest"
