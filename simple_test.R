library(causalForest)
n = 500
p = 10

X = matrix(rnorm(n * p), n, p)
colnames(X) = sample.int(p)
Y = rnorm(n)
DF = data.frame(X=X, Y=Y)

W = rbinom(n, 1, 0.5)

#
# Test basic tree functionality
#

ct = causalTree(Y ~ ., data = DF, treatment=W, split.option="TOT", cv.option="TOT")
#summary(ct)
min(predict(ct))
max(predict(ct))
class(ct) == "causalTree"

rct = refit.causalTree(ct, newx=DF, newy=Y, treatment=W, propensity=rep(0.5, n))
min(predict(rct))
max(predict(rct))
class(rct) == "causalTree"

ct = causalTree(Y ~ ., data = DF, treatment=W, split.option="CT", cv.option="TOT")
#summary(ct)
min(predict(ct))
max(predict(ct))
class(ct) == "causalTree"

rct = refit.causalTree(ct, newx=DF, newy=Y, treatment=W)
min(predict(rct))
max(predict(rct))
class(rct) == "causalTree"


#
# Test refitting with small, potentially empty, leaves
#

n2 = 25
X2 = matrix(rnorm(n2 * p), n2, p)
colnames(X) = sample.int(p)
Y2 = rnorm(n2)
DF2 = data.frame(X=X2, Y=Y2)

W2 = rbinom(n2, 1, 0.5)

ct <- causalTree(Y ~ ., data = DF, treatment = W, method = "anova", cp = 0, minbucket = 1, cv.option = "matching", split.option = "CT", xval = 0)

rct = refit.causalTree(ct, newx=DF2, newy=Y2, treatment=W2)
sum(is.na(predict(rct))) == 0

rct.tot = refit.causalTree(ct, newx=DF2, newy=Y2, treatment=W2)
sum(is.na(predict(rct.tot))) == 0

#
# Test forests
#


cf = causalForest(DF[,1:p], DF$Y, W, num.trees = 10)
min(predict(cf, data.frame(X=X)))
max(predict(cf, data.frame(X=X)))
class(cf) == "causalForest"

pf = propensityForest(DF[,1:p], DF$Y, W, num.trees = 10)
min(predict(pf, data.frame(X=X)))
max(predict(pf, data.frame(X=X)))
class(pf) == "causalForest"