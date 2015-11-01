#' Evalute the estimates for tau made by an already-trained causalForest.
#'
#' @param forest   the fitted causalForest object
#' @param newdata  the new test points at which the causalForest predictions
#'                 are to be evaluated
#'
#' @return estimates for tau, corresponding to each row of newdata

predict.causalForest <- function(forest, newdata, predict.all = FALSE) {
  if (!inherits(forest, "causalForest")) stop("Not a legitimate \"causalForest\" object")  
  test.data <- data.frame(X=newdata)
  individual <- sapply(forest$trees, function(tree.fit) {
  	predict(tree.fit, test.data)
  })
  
  aggregate <- rowMeans(individual)
  if (predict.all) {
    list(aggregate = aggregate, individual = individual)
  } else {
    aggregate
  }
}
