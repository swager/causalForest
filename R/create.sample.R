create.sample <- function(Y, X, W, sample.size, replace) {
  num.obs <- nrow(X)
  sample.indices <- sample(1:num.obs, sample.size, replace = replace)
  sample.Y <- Y[sample.indices]
  sample.X <- X[sample.indices,]
  sample.W <- W[sample.indices]
  list(Y = sample.Y, X = sample.X, W = sample.W, indices = sample.indices)
}