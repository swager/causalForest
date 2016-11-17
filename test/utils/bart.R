library(BayesTree)

bart.cate = function(X, Y, W, X.test) {
	
	X.all = cbind(W, X)
	X.dummy = cbind(c(rep(0, nrow(X.test)), rep(1, nrow(X.test))), rbind(X.test, X.test))
	bart.fit = bart(X.all, Y, X.dummy)
	bart.preds = colMeans(bart.fit$yhat.test)
	bart.tau = bart.preds[nrow(X.test) + 1:nrow(X.test)] - bart.preds[1:nrow(X.test)]	
	bart.tau
	
}