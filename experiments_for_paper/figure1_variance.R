library(causalForest)
library(randomForestCI)

rm(list = ls())

sigma = 1
d = 20
nvals = c(100, 200, 400, 800, 1600)

tree_mult = 4
n.test = 1000
simu.reps = 50

baseline = function(x) {
	2 * (x[1] - 0.5)
}

propensity = function(x) {
	0.25 + dbeta(x[1], 2, 4)/4
}

X.test = matrix(runif(n.test * d, 0, 1), n.test, d)

single.run = function(n) {
	X = matrix(runif(n * d, 0, 1), n, d) # features
	e = apply(X, 1, propensity)
	W = rbinom(n, 1, e) #treatment condition
	
	# no treatment effect
	Y = apply(X, 1, baseline) + sigma * rnorm(n)
	
	#
	# random forest
	#
	
	forest = propensityForest(X, Y, W, num.trees = round(tree_mult * n), sample.size = n^(0.8), nodesize = 1)
	predictions = predict(forest, X.test)
	forest.ci = randomForestInfJack(forest, X.test, calibrate = TRUE)
	
	return(forest.ci)
}

simu.fun = function(n, reps) {
	
	results = lapply(1:reps, function(rr)single.run(n))
	preds = Reduce(cbind, lapply(results, function(xx) xx$y.hat))
	variances = apply(preds, 1, var)
	
	var.hat = Reduce(cbind, lapply(results, function(xx) xx$var.hat))
	err = var.hat - variances
	
	return(list(n=n, variances=variances, err=err))
}

results.raw = lapply(nvals, function(n) {
	print(paste("NOW RUNNING:", n))
	simu.fun(n, simu.reps)
})

save.image("figure2_variance.RData")

load("figure2_variance.RData")

variances = data.frame(Reduce(rbind, lapply(results.raw, function(xx)cbind(n=xx[[1]], V=xx[[2]]))))

pdf("output/variance_decay.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
boxplot(V ~ n, variances, ylim = c(0, max(variances$V)), xlab = "n", ylab = "forest sampling variance")
abline(h=0, lty = 3)
par=pardef
dev.off()

coef_var = data.frame(Reduce(rbind, lapply(results.raw, function(xx)cbind(n=xx[[1]], CV=sqrt(rowMeans(xx[[3]]^2)) / xx[[2]]))))

pdf("output/ij_coef_variation.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
boxplot(CV ~ n, coef_var, ylim = c(0, max(coef_var$CV)), xlab = "n", ylab = "inf. jack. coefficient of variation")
abline(h=0, lty = 3)
par=pardef
dev.off()


bias = data.frame(Reduce(rbind, lapply(results.raw, function(xx)cbind(n=xx[[1]], B=rowMeans(xx[[3]])))))
pdf("output/ij_bias.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
boxplot(B ~ n, bias, xlab = "n", ylab = "inf. jack. bias")
abline(h=0, lty = 3)
par=pardef
dev.off()

