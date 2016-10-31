# requires the version on github at susanathey/causalTree
loadNamespace("causalTree")

causaltree.cate = function(X, Y, W, X.test) {

	train.idx = rbinom(ncol(X), 1, 0.5)
	
	data.all = data.frame(X = X, Y = Y)               
	data.train = data.all[train.idx==1,]   
	est_data = data.all[train.idx==0,]  
	treat.train = W[train.idx==1]
	est_treatment = W[train.idx==0] 
	       
	honestTree <- causalTree::honest.causalTree(Y ~ ., 
	                                data = data.train,
	                                treatment = treat.train,
	                                est_data = est_data,
	                                est_treatment = est_treatment,
	                                split.Rule = "CT",
	                                split.Honest = TRUE,
	                                HonestSampleSize = nrow(est_data),
	                                split.Bucket = TRUE,
	                                cv.option = "fit",
	                                cv.Honest = FALSE)
	                                
	opcp <- honestTree$cptable[,1][which.min(honestTree$cptable[,4])]
	opTree <- prune(honestTree, opcp)
	
	predict(opTree, data.frame(X=X.test))
}
