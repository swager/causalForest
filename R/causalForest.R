causalForest <- function(X, Y, W, num.trees, sample.size = floor(length(Y) / 10), mtry = ceiling(ncol(X)/3), nodesize = 1) {
  
  if (any(is.na(X)) || any(is.na(Y)) || any(is.na(W))) {
    stop("There are missing values in the input.")
  }

  num.obs <-nrow(X)
  causalForest.honest <- init.causalForest(X, Y, W, num.trees)
  sample.size <- min(sample.size, floor(num.obs / 2))
  
  print("Building trees ...")
  
  for (tree.index in 1:num.trees) {
    
    print(paste("Tree", as.character(tree.index)))
    
    full.idx <- sample.int(num.obs, 2 * sample.size, replace = FALSE)
    train.idx <- full.idx[1:sample.size]
    reestimation.idx <- full.idx[sample.size + (1:sample.size)]
    
    tree.standard <- causalTree(Y ~ ., data = data.frame(X = X[train.idx,], Y = Y[train.idx]), treatment = W[train.idx], method = "anova", cp = 0, minbucket = nodesize, cv.option = "matching", split.option = "CT", xval = 0)
    
    tree.honest <- refit.causalTree(tree.standard, treatment=W[reestimation.idx], newdata=data.frame(X = X[reestimation.idx,], Y = Y[reestimation.idx]))
    
    causalForest.honest$trees[[tree.index]] <- tree.honest
    causalForest.honest$inbag[full.idx, tree.index] <- 1
  }
  
  return(causalForest.honest)
}
