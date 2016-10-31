library(randomForest)
library(causalForest)

rm(list = ls())

n = 1000
p = 20
beta = 2

xx = seq(-3, 3, by = 0.01)
X.test = matrix(0, length(xx), p)
X.test[,1] = xx

twoforest = replicate(10, {
	
X = matrix(rnorm(n * p), n, p)
W = rbinom(n, 1, 0.1)
Y = beta * X[,1] + rnorm(n)
	
rf.0 = randomForest(X[W==0,], Y[W==0], X.test, replace = FALSE, sampsize = sum(W==0)^0.8)
rf.1 = randomForest(X[W==1,], Y[W==1], X.test, replace = FALSE, sampsize = sum(W==1)^0.8)
rf.1$test$predicted - rf.0$test$predicted
})

singleforest = replicate(10, {
	
X = matrix(rnorm(n * p), n, p)
W = rbinom(n, 1, 0.1)
Y = beta * X[,1] + rnorm(n)
	
cf = causalForest(X, Y, W, num.trees = 500, sample.size = length(W)^0.8, honest=TRUE)
predict(cf, X.test)
})

pdf('~/git_local/causalForest-mirror/test/two_forests.pdf')
plot(NA, NA, xlim = range(xx), ylim = range(twoforest), xlab = "X1", ylab="tau.hat", main="RCT with no treatment effect \n There is a strong main effect, 1/10 observations are treated")
for(iter in 1:10) {
	lines(xx, singleforest[,iter], col = 1)
}
for(iter in 1:10) {
	lines(xx, twoforest[,iter], col = 2)
}
legend("topright", c("Two Forests", "Causal Forest"), lwd = 1, col=2:1)
dev.off()