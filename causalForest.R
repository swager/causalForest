causalForest <- function(X, Y, W, weights = rep(1, nrow(X)), num.trees, sample.size = floor(length(Y) / 10), mtry = ceiling(ncol(X)/3), nodesize = 1, cores = NULL, verbose = TRUE, clustvar = NULL) {

  if (any(is.na(X)) || any(is.na(Y)) || any(is.na(W)) || any(is.na(weights))) {
    stop("There are missing values in the input.")
  }

  num.obs <-nrow(X)
  if (is.null(weights)){weights <- rep(1, num.obs)}
  causalForest.honest <- init.causalForest(X, Y, W, weights, num.trees)
  sample.size <- min(sample.size, floor(num.obs / 2))
    
  if (verbose) print("Building trees ...")
  
  init.seed = sample.int(268435456, 1)
  
  if (is.null(cores)) {
    forest.raw <- foreach::foreach(tree.index = 1:num.trees) %do%
	  causalForest.getTree(tree.index, X, Y, W, weights, num.obs, sample.size, mtry, nodesize, init.seed, verbose, clustvar)
  } else {
    if (Sys.info()[1] == "Windows"){
      registerDoParallel(cores=cores)
    } else {
      registerDoMC(cores=cores)
    }
  	forest.raw <- foreach::foreach(tree.index = 1:num.trees) %dopar%
	  causalForest.getTree(tree.index, X, Y, W, weights, num.obs, sample.size, mtry, nodesize, init.seed, verbose, clustvar)
  }
  
  for (tree.index in 1:num.trees) {
  	causalForest.honest$trees[[tree.index]] <- forest.raw[[tree.index]][[2]]
  	causalForest.honest$inbag[forest.raw[[tree.index]][[1]], tree.index] <- 1
  }
  
  return(causalForest.honest)
}

causalForest.getTree <- function(tree.index, X, Y, W, weights, num.obs, sample.size, mtry, nodesize, init.seed, verbose, clustvar = NULL) {
	if(verbose) print(paste("Tree", as.character(tree.index)))
	set.seed(init.seed + tree.index)

  if (is.null(clustvar)){
    full.idx <- sample.int(num.obs, 2 * sample.size, replace = FALSE)
    train.idx <- full.idx[1:sample.size]
    reestimation.idx <- full.idx[sample.size + (1:sample.size)]    
  } else {
    if (any(is.na(clustvar))){stop("Missing values in cluster indicator")}
    `%ni%` <- Negate(`%in%`)
    tab <- tapply(clustvar,clustvar, length)
    clusters <- data.frame(id = names(tab), size = tab)    
    clusters <- clusters[sample(1:nrow(clusters), size =nrow(clusters), replace = FALSE),]
    clusters$cs <- cumsum(clusters$size)
    full.cl <- clusters$id[clusters$cs <= (2*sample.size)]
    train.cl <- clusters$id[clusters$cs < (sample.size)]
    reestimation.cl <- full.cl[full.cl %ni% train.cl]
    full.idx <- c()
    for (i in 1:length(full.cl)){
      full.idx <- append(full.idx, which(clustvar %in% full.cl[i]))
      if (i == length(train.cl)){
        train.idx <- full.idx
      }
    }
    reestimation.idx <- full.idx[(length(train.idx)+1): length(full.idx)]
  }
   
  tree.DF = data.frame(X = X, Y = Y)
    
  tree.standard <- causalTree(Y ~ ., data = tree.DF[train.idx,], treatment = W[train.idx], weights = weights[train.idx], method = "anova", cp = 0, minbucket = nodesize, cv.option = "matching", split.option = "CT", xval = 0)


  tree.honest <- refit.causalTree(tree.standard, newx=tree.DF[reestimation.idx,], newy = Y[reestimation.idx], treatment=W[reestimation.idx], weights=weights[reestimation.idx])

  return(list(full.idx, tree.honest))
}
