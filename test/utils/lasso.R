library(glmnet)

lasso.cate = function(X, Y, W, X.test) {
	
	X.W = diag(2 * W - 1) %*% X
	X.all = cbind(W, X, X.W)
	glmnet.fit = cv.glmnet(X.all, Y, penalty.factor = c(0, rep(1, 2 * ncol(X))))
	pred.control = predict(glmnet.fit, newx=cbind(0, X, -X))
	pred.treat = predict(glmnet.fit, newx=cbind(1, X, X))
	
	as.numeric(pred.treat - pred.control)

}