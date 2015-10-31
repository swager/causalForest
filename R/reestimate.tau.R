# Reestimate the output at each leaf using the Y of a new set of observations.

reestimate.tau <- function(fit, Y, X, W) {
  where <- est.causalTree(fit, causalTree.matrix(create.data.frame(fit, X)))
  unique.leaves <- unique(where)
  treat <- which(W == 1)
  control <- which(W == 0)
  all.leaves <- as.numeric(row.names(fit$frame)[which(fit$frame$var == "<leaf>")])
  for (leaf in unique.leaves) {
    tau <- estimate.leaf.tau(where, treat, control, Y, all.leaves, leaf)
    if (!is.null(tau) && !is.nan(tau) && !is.na(tau)) {
      fit$frame[which(row.names(fit$frame) == as.character(leaf)), 'yval'] <- tau
    }
  }
  fit
}