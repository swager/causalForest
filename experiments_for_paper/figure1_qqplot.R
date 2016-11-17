library(causalForest)

rm(list = ls())

sigma = 1
d = 20
n = 800

tree_mult = 4
n.test = 1000
simu.reps = 20

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
	predictions
}


results = lapply(1:simu.reps, function(rr)single.run(n))

save.image("figure2b_gauss.RData")

preds = Reduce(cbind, results)
preds.std = t(apply(preds, 1, function(xx) (xx - mean(xx)) / sd(xx)))
ymax = max(abs(preds.std))

gauss.quantiles = qnorm((1:simu.reps) / (simu.reps + 1))

#plot(NA, NA, xlim = range(gauss.quantiles), ylim = c(-ymax, ymax), xlab = "theoretical quantiles", ylab = "sample quantiles")
#for(iter in 1:nrow(preds.std)) {
#	points(gauss.quantiles, sort(preds.std[iter,]))
#}
#abline(0, 1, lwd = 3, col = 2)

all.data = data.frame(Reduce(rbind, lapply(1:nrow(preds.std),
  function(iter) cbind(Q= gauss.quantiles, T=sort(preds.std[iter,]))
)))


pdf("output/gauss_qqplot.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(NA, NA, xlim = range(c(gauss.quantiles, -2, 2)), ylim = c(-ymax, ymax), xlab = "theoretical quantiles", ylab = "sample quantiles")
abline(0, 1, lwd = 2, lty = 3)
boxplot(T~Q, all.data, xaxt="n", at = gauss.quantiles, boxwex=0.2, add=TRUE)
xaxes = seq(-2, 2, by = 1)
axis(1, at=xaxes, label=xaxes)
par=pardef
dev.off()



