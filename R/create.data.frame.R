create.data.frame <- function(fit, X) {
  data <- data.frame(X)
  data.col.names <- rep("", ncol(X))
  for (col in 1:ncol(X)) {
    data.col.names[col] <- paste("X", as.character(col), sep = ".")
  }
  names(data) <- data.col.names
  data$y <- rep(0, nrow(X))
  if (is.null(attr(data, "terms"))) {
    Terms <- delete.response(fit$terms)
    data <- model.frame(Terms, data, na.action = na.pass, xlev = attr(fit, "xlevels"))
    if (!is.null(cl <- attr(Terms, "dataClasses")))
      .checkMFClasses(cl, data, TRUE)
  }
  data
}