rm(list = ls())

library(causalForest)

p = 10

ntree = 500
n = 100 * 2^5
reps = 40
n.test = 400

tau.star = 0.1

#results.all = lapply(nvals, function(n) {

#    print(paste("NOW AT N", n))

	res = replicate(reps, {

		X = matrix(runif(n*p), n, p)
		W = rbinom(n, 1, 0.5)
		Y = 2 * rbinom(n, 1, 0.05) * W + 0.1 * rnorm(n)
		
		data = data.frame(X=X, Y=Y)
		
		X.test = matrix(runif(n.test * p), n.test, p)
		
		cf.h = causalForest(X, Y, W, num.trees= ntree, sample.size=n^0.8, honest=TRUE)
		honest=predict(cf.h, X.test)
		honest.mse = mean((honest - tau.star)^2)
		
		cf.a = causalForest(X, Y, W, num.trees= ntree, sample.size=n^0.8, honest=FALSE)
		adapt = predict(cf.a, X.test)
		adapt.mse=mean((adapt - tau.star)^2)
		
		c(honest.mse, adapt.mse)
	})
#})

rowMeans(res)

save.image(file="honesty_evaluation_globalMSE.RData")