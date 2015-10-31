## TODO: 
## The notion we really want is that comparison forests implement
## a random forest like interface

init.causalForest <- function(x, y, w, num.trees) {
  num.obs <- nrow(x)
  trees <- vector("list", num.trees)
  inbag <- matrix(0, num.obs, num.trees)
  causalForest <- list(trees = trees, x = x, y = y, w = w, ntree = num.trees, inbag = inbag)
  class(causalForest) <- "causalForest"
  causalForest
}