# only suitable for full binary tree:
get.descendant.leaves <- function(parent, leaves) {
  max <- floor(log(max(leaves), 2))
  result <- get.descendant.leaves.helper(parent, leaves, 0, max)
  return (result)
}

get.descendant.leaves.helper <- function(parent, leaves, count, maxdepth) {
  if (count > maxdepth) {
    # rarely happens
    stop ("This node is not in the tree.")
  } else {
    if (!is.na(match(parent, leaves))) {
      return (parent)
    } else {
      left_son <- 2 * parent
      right_son <-  2 * parent + 1
      left_sons <- get.descendant.leaves.helper(left_son, leaves, count + 1, maxdepth)
      right_sons <- get.descendant.leaves.helper(right_son, leaves, count + 1, maxdepth)
      result <- c(left_sons, right_sons)
    }
    return(result)    
  } 
}

contains.all.treatment.levels <- function(index, treatment) {
  if (missing(treatment)) {
    TRUE
  } else {
    treat0.obs <- which(treatment == 0)
    treat1.obs <- which(treatment == 1)
    length(which(index %in% treat0.obs)) > 0 && length(which(index %in% treat1.obs)) > 0
  }
}

recursive.which.in.leaf <- function(leaf.assignments, leaf, leaves, treatment) {
  in.leaf <- which(leaf.assignments == leaf)
  parent <- leaf
  while(length(in.leaf) == 0 || !contains.all.treatment.levels(in.leaf, treatment)) {
    parent <- floor(parent / 2)
    descendant.leaves <- get.descendant.leaves(parent, leaves)
    in.leaf <- which(leaf.assignments %in% descendant.leaves)
  }
  in.leaf
}

estimate.leaf.tau <- function(leaf.assignments, treat, control, Y, leaves, leaf) {
  index <- which(leaf.assignments == leaf)
  index1 <- intersect(index, treat)
  index0 <- intersect(index, control)    
  tau <- mean(Y[index1]) - mean(Y[index0])
  parent <- leaf
  while(is.na(tau)){
    # go back to parent node who can compute the value:
    parent <- floor(parent / 2)
    descendant.leaves <- get.descendant.leaves(parent, leaves)
    obs.in.parent<- which(leaf.assignments %in% descendant.leaves)
    index1 <- intersect(obs.in.parent, treat)
    index0 <- intersect(obs.in.parent, control)
    tau <- mean(Y[index1]) - mean(Y[index0])
  }
  tau
}

#' Estimate the causal effects using honest tree model.
#' 
#' @param object A tree-structured fit object.
#' @param formula A regression formula.
#' @param data New observations.
#' @param treatment The weights status of new observations
#' @return The estimated causal effects of \code{data}.
#' Notice here when the leaf contains only treated or control cases, the function will trace back to the leaf's parent mnode recursively until the parent can be used to compute causal effect.
#' 
## estimate function for honest causal tree:
estimate.causalTree <- function(object, formula, data, treatment, na.action = na.pass)
{
  if (!inherits(object, "causalTree")) stop("Not a legitimate \"causalTree\" object")
  Call <- match.call()
  indx <- match(c("formula", "data", "treatment"),
                names(Call), nomatch = 0L)
  if (indx[1] == 0L) stop("a 'formula' argument is required")
  temp <- Call[c(1L, indx)]      # only keep the arguments we wanted
  temp$na.action <- na.action    # This one has a default
  temp[[1L]] <- quote(stats::model.frame) # change the function called
  m <- eval.parent(temp)
  Y <- model.response(m)
  n <- nrow(m)
  # check the treatment condition:
  if (missing(treatment)) stop("You should import the treatment status for data.")
  if (length(treatment) != n) 
    stop("The length of treatment status vector should be same as number
         of observations.")
  if (length(which(treatment == 0)) == 0 || length(which(treatment == 1)) == 0)
    stop("Can't make estimation since no treated cases or no control cases.")
  # get the leaf of the object
  leaves <- as.numeric(row.names(object$frame)[which(object$frame$var == "<leaf>")])
  
  # get the node id for each observation.
  where <- {
    if (is.null(attr(data, "terms"))) {
      Terms <- delete.response(object$terms)
      data <- model.frame(Terms, data, na.action = na.action,
                             xlev = attr(object, "xlevels"))
      if (!is.null(cl <- attr(Terms, "dataClasses")))
        .checkMFClasses(cl, data, TRUE)
    }
    est.causalTree(object, causalTree.matrix(data))
  }
  
  unique_leaf <- unique(where)
  causal_estimation <- rep(0, n)
  treat <- which(treatment == 1)
  control<- which(treatment == 0)
  for (leaf in unique_leaf) {
    index <- which(where == leaf)
    causal_estimation[as.numeric(index)] <- estimate.leaf.tau(where, treat, control, Y, leaves, leaf)
  }
  return(causal_estimation)  
}
