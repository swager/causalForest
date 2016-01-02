causalForest <- function(X, Y, W, num.trees, sample.size = floor(length(Y) / 10), mtry = ceiling(ncol(X)/3), nodesize = 1, cores = NULL, verbose = TRUE) {
  
  if (any(is.na(X)) || any(is.na(Y)) || any(is.na(W))) {
    stop("There are missing values in the input.")
  }

  num.obs <-nrow(X)
  causalForest.honest <- init.causalForest(X, Y, W, num.trees)
  sample.size <- min(sample.size, floor(num.obs / 2))
  
  if (verbose) print("Building trees ...")
  
  init.seed = sample.int(268435456, 1)
  
  if (is.null(cores)) {
  	
    forest.raw <- foreach::foreach(tree.index = 1:num.trees) %do%
	  causalForest.getTree(tree.index, X, Y, W, num.obs, sample.size, mtry, nodesize, init.seed, verbose)
	  
  } else {
  	
  	registerDoMC(cores=cores)
  	forest.raw <- foreach::foreach(tree.index = 1:num.trees) %dopar%
	  causalForest.getTree(tree.index, X, Y, W, num.obs, sample.size, mtry, nodesize, init.seed, verbose)

  }
  
  for (tree.index in 1:num.trees) {
  	causalForest.honest$trees[[tree.index]] <- forest.raw[[tree.index]][[2]]
  	causalForest.honest$inbag[forest.raw[[tree.index]][[1]], tree.index] <- 1
  }
  
  return(causalForest.honest)
}

causalForest.getTree <- function(tree.index, X, Y, W, num.obs, sample.size, mtry, nodesize, init.seed, verbose) {
	
	if(verbose) print(paste("Tree", as.character(tree.index)))
	set.seed(init.seed + tree.index)
    
    full.idx <- sample.int(num.obs, 2 * sample.size, replace = FALSE)
    train.idx <- full.idx[1:sample.size]
    reestimation.idx <- full.idx[sample.size + (1:sample.size)]
    
    tree.DF = data.frame(X = X, Y = Y)
    
    tree.standard <- causalTree(Y ~ ., data = tree.DF[train.idx,], treatment = W[train.idx], method = "anova", cp = 0, minbucket = nodesize, cv.option = "matching", split.option = "CT", xval = 0)
    
    tree.honest <- refit.causalTree(tree.standard, newx=tree.DF[reestimation.idx,], newy = Y[reestimation.idx], treatment=W[reestimation.idx])
    
    return(list(full.idx, tree.honest))
}