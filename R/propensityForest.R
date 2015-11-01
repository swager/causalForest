propensityForest <- function(X, Y, W, num.trees, sample.size = floor(length(Y) / 10), mtry = ceiling(ncol(X)/3), nodesize = 1) {
  
  if (any(is.na(X)) || any(is.na(Y)) || any(is.na(W))) {
    stop("There are missing values in the input.")
  }
 
  num.obs <-nrow(X)
  causalForest.honest <- init.causalForest(X, Y, W, num.trees)
  sample.size <- min(sample.size, floor(num.obs / 2))
  
  print("Building trees ...")
  
  for (tree.index in 1:num.trees) {
    
    print(paste("Tree", as.character(tree.index)))
    
    full.idx <- sample.int(num.obs, sample.size, replace = FALSE)
    
    # Note: Here, we want to train a causal tree, but pick the splits by
    # maximizing classification accuracy for the treatment assignments.
    # (The goal is to have observations within a leaf all have the same treatment
    # propensity).
    #
    # Ideally, we would add an option "split = propensity" to the causalTree call;
    # for now, we just get the splits from rpart.
    tree.propensity <- rpart(W ~ ., data = data.frame(X = X[full.idx,], W = W[full.idx]), method = "class", minbucket = nodesize)
    class(tree.propensity) <- "causalTree"
    
    # This is a horrible hack to make predict.causalTree do the right thing. Sorry.
    attr(tree.propensity, "ylevels") = NULL
    
    tree.honest <- reestimate.tau(tree.propensity, Y[full.idx], data.frame(X = X[full.idx,]), W[full.idx])
    
    causalForest.honest$trees[[tree.index]] <- tree.honest
    causalForest.honest$inbag[full.idx, tree.index] <- 1
  }
  
  return(causalForest.honest)
}
