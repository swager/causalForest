modify.estimations <- function(tree, XW, Y, model.name, propensity) {
  X <- XW$X
  W <- XW$W
  if (model.name == "CT") {
    estimated.leaves <- est.causalTree(tree, X)
    leaves <- which(tree$frame$var == "<leaf>")
    for (leaf in leaves) {
      obs.in.leaf <- which(estimated.leaves == leaf)
      relevant.W <- W[obs.in.leaf] 
      relevant.Y <- Y[obs.in.leaf]
      tree$frame$yval[leaf] <- sum(relevant.Y * relevant.W / sum(relevant.W)) - sum(relevant.Y * (1 - relevant.W) / sum(1 - relevant.W))
    }
  } else if (model.name == "TOT") {
    estimated.leaves <- est.causalTree(tree, X)
    transformed.Y <- Y * (W - propensity) / (propensity * (1 - propensity))
    leaves <- which(tree$frame$var == "<leaf>")
    for (leaf in leaves) {
      tree$frame$yval[leaf] <- mean(transformed.Y(which(estimated.leaves == leaf)))
    }
  } else if (model.name == "TT") {
    reestimate(tree, X, W, Y)
  } else if (model.name == "ST") {
    reestimate(tree, X, W, Y)
  } else {
    stop("model.name must be 'ST', 'TT', 'TOT', or 'CT'")
  }
  tree
}