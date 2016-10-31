library(causalForest)
library(Hmisc)
library(mgcv)
library(ggplot2)
library(randomForestCI)

rm(list = ls())
source("utils/knn.R")

n = 10000
ntree = 10000
n.test = 10000

sigma = 1
d = 20

effect = function(x) {
  4/((1 + exp(-12 * (x[1] - 0.5))) * (1 + exp(-12 * (x[2] - 0.5))))
}

baseline = function(x) { 0 }

X = matrix(runif(n * d, 0, 1), n, d) # features
W = rbinom(n, 1, 0.5) #treatment condition
Y = apply(X, 1, baseline) +  (W - 0.5) * apply(X, 1, effect) + sigma * rnorm(n)

X.test = matrix(runif(n.test * d, 0, 1), n.test, d)
true.eff = apply(X.test, 1, effect)

#
# dry run
#


knn.cate(X, Y, W, X.test, k = 20)


#
# do d = 20 case
#

forest = causalForest(X, Y, W, num.trees = ntree, sample.size = n / 10)
tauhat.rf <- predict(forest, X.test)

kk = c(seq(3, 99, by = 3))
knn.mses = sapply(kk, function(k) {
        tauhat = knn.cate(X, Y, W, X.test, k = k)
        mean((true.eff - tauhat)^2)
})

print(knn.mses)

k.opt = kk[which.min(knn.mses)]
tauhat.knn = knn.cate(X, Y, W, X.test, k = k.opt)

print(min(knn.mses))
print(k.opt)

print(mean((true.eff - tauhat.rf)^2))
print(min(knn.mses) / mean((true.eff - tauhat.rf)^2))

minp = min(true.eff, tauhat.rf, tauhat.knn)
maxp = max(true.eff, tauhat.rf, tauhat.knn)
rngp = maxp - minp

ncol = 100

true.scl = pmax(ceiling(ncol * (true.eff - minp) / rngp), 1)
rf.scl = pmax(ceiling(ncol * (tauhat.rf - minp) / rngp), 1)
knn.scl = pmax(ceiling(ncol * (tauhat.knn - minp) / rngp), 1)

hc = heat.colors(ncol)

pdf("output/cate_true_20d.pdf")
plot(X.test[,1], X.test[,2], pch = 16, col = hc[true.scl], xlab = "", ylab = "")
dev.off()

pdf("output/cate_rf_20d.pdf")
plot(X.test[,1], X.test[,2], pch = 16, col = hc[rf.scl], xlab = "", ylab = "")
dev.off()

pdf("output/cate_knn_20d.pdf")
plot(X.test[,1], X.test[,2], pch = 16, col = hc[knn.scl], xlab = "", ylab = "")
dev.off()

#
# do d = 6 case
#

X = X[,1:6]
X.test = X.test[,1:6]

forest = causalForest(X, Y, W, num.trees = ntree, sample.size = n / 10)
tauhat.rf <- predict(forest, X.test)

kk = c(seq(3, 99, by = 3))
knn.mses = sapply(kk, function(k) {
        tauhat = knn.cate(X, Y, W, X.test, k = k)
        mean((true.eff - tauhat)^2)
})

print(knn.mses)

k.opt = kk[which.min(knn.mses)]
tauhat.knn = knn.cate(X, Y, W, X.test, k = k.opt)

print(min(knn.mses))
print(k.opt)

print(mean((true.eff - tauhat.rf)^2))
print(min(knn.mses) / mean((true.eff - tauhat.rf)^2))

# do not recompute!
# minp = min(true.eff, tauhat.rf, tauhat.knn)
# maxp = max(true.eff, tauhat.rf, tauhat.knn)
# rngp = maxp - minp

# ncol = 100
# hc = heat.colors(ncol)

true.scl = pmax(ceiling(ncol * (true.eff - minp) / rngp), 1)
rf.scl = pmax(ceiling(ncol * (tauhat.rf - minp) / rngp), 1)
knn.scl = pmax(ceiling(ncol * (tauhat.knn - minp) / rngp), 1)

pdf("output/cate_true_6d.pdf")
plot(X.test[,1], X.test[,2], pch = 16, col = hc[true.scl], xlab = "", ylab = "")
dev.off()

pdf("output/cate_rf_6d.pdf")
plot(X.test[,1], X.test[,2], pch = 16, col = hc[rf.scl], xlab = "", ylab = "")
dev.off()

pdf("output/cate_knn_6d.pdf")
plot(X.test[,1], X.test[,2], pch = 16, col = hc[knn.scl], xlab = "", ylab = "")
dev.off()

