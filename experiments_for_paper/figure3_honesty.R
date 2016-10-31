rm(list = ls())

library(causalForest)

p = 10

ntree = 500
nvals = round(100 * 2^seq(0, 5, length.out=10))
reps = 40

results.all = lapply(nvals, function(n) {

    print(paste("NOW AT N", n))

	replicate(reps, {

		X = matrix(runif(n*p), n, p)
		W = rbinom(n, 1, 0.5)
		Y = 2 * rbinom(n, 1, 0.05) * W + 0.1 * rnorm(n)
		
		data = data.frame(X=X, Y=Y)
		
		pred.point = data.frame(X=matrix(0, 2, p))
		names(pred.point) = 1:p
		
		cf.h = causalForest(X, Y, W, num.trees= ntree, sample.size=n^0.8, honest=TRUE)
		cf.a = causalForest(X, Y, W, num.trees= ntree, sample.size=n^0.8, honest=FALSE)
		
		c(honest=predict(cf.h, pred.point)[1], adapt=predict(cf.a, pred.point)[1])
	})
})

CATE.true = 0.1
RMSE = sapply(1:length(nvals), function(idx) {
	sqrt(rowMeans((results.all[[idx]] - CATE.true)^2))
})
BIAS = sapply(1:length(nvals), function(idx) {
	rowMeans(results.all[[idx]] - CATE.true)
})

save.image(file="honesty_evaluation.RData")


labvec = c(100, 200, 400, 800, 1600, 3200)

pdf("output/honesty_evaluation.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(NA, NA, xaxt = "n", xlim = range(nvals), ylim = c(0, max(RMSE, BIAS) + 0.05), xlab = "n", ylab = "error", log = "x")
lines(nvals, RMSE[1,], col = 4, lwd = 2, lty = 2)
lines(nvals, RMSE[2,], col = 2, lwd = 2, lty = 2)
lines(nvals, BIAS[1,], col = 4, lwd = 2, lty = 1)
lines(nvals, BIAS[2,], col = 2, lwd = 2, lty = 1)
abline(h=0, lty=3)
legend("topleft", c("Honest bias", "Honest RMSE", "Adapt. bias", "Adapt. RMSE"), lty = c(1, 2, 1, 2), col = c(4, 4, 2, 2), lwd = 2, cex=1.5)
axis(1, at=labvec, labels=labvec)
par = pardef
dev.off()
