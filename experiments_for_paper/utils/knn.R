library(FNN)

knn.cate = function(X, Y, W, X.test, k = 10) {
	yhat.0 = knn.reg(train=X[W==0,], test=X.test, y=Y[W==0], k = k)
	yhat.1 = knn.reg(train=X[W==1,], test=X.test, y=Y[W==1], k = k)
	return(yhat.1$pred - yhat.0$pred)
}