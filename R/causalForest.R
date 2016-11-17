causalForest <- function(X, Y, W, num.trees, sample.size = floor(length(Y) / 10), nodesize = 1, honest = TRUE) {
  
  if (any(is.na(X)) || any(is.na(Y)) || any(is.na(W))) {
    stop("There are missing values in the input.")
  }

  num.obs <-nrow(X)
  causalForest <- init.causalForest(X, Y, W, num.trees)
  sample.size <- min(sample.size, floor(num.obs / 2))
  
  print("Building trees ...")
  
  for (tree.index in 1:num.trees) {
    
    print(paste("Tree", as.character(tree.index)))
    
    if (honest) {
    	
	    full.idx <- sample.int(num.obs, 2 * sample.size, replace = FALSE)
	    train.idx <- full.idx[1:sample.size]
	    reestimation.idx <- full.idx[sample.size + (1:sample.size)]
	    
	    tree.DF = data.frame(X = X, Y = Y)
	    
	    tree.standard <- causalTree(Y ~ ., data = tree.DF[train.idx,], treatment = W[train.idx], method = "anova", cp = 0, minbucket = nodesize, cv.option = "matching", split.option = "CT", xval = 0)
	    
	    tree.honest <- refit.causalTree(tree.standard, newx=tree.DF[reestimation.idx,], newy = Y[reestimation.idx], treatment=W[reestimation.idx])
	    
	    causalForest$trees[[tree.index]] <- tree.honest
	    causalForest$inbag[full.idx, tree.index] <- 1
     
     } else {
     	
     	train.idx <- sample.int(num.obs, sample.size, replace = FALSE)
	    
	    tree.DF = data.frame(X = X, Y = Y)
	    
	    tree.standard <- causalTree(Y ~ ., data = tree.DF[train.idx,], treatment = W[train.idx], method = "anova", cp = 0, minbucket = nodesize, cv.option = "matching", split.option = "CT", xval = 0)
	    
	    causalForest$trees[[tree.index]] <- tree.standard
	    causalForest$inbag[train.idx, tree.index] <- 1
	    
	 }
  }
  
  return(causalForest)
}
