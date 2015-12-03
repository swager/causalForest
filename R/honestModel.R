#### honest model for ST, TT, TOT and CT:
# honest ST model 
honestST <- function(object, data, na.action = na.rpart) {

  warning( "The function honestST is deprecated, and will be removed. Please use refit.causalTree instead." )

  #if (!inherits(object, "rpart")) stop("Not a legitimate \"rpart\" object")
  
  if (!("treatment" %in% names(data))) stop("Function requires treatment column to be named 'treatment' in dataframe")
  
  nodes <- as.numeric(row.names(object$frame))
  num <- length(nodes)
  Terms <- object$terms
  #data <- model.frame(Terms, data, na.action = na.action,
  data <- model.frame(Terms, data, na.action = na.action, 
                      xlev = attr(object, "xlevels"))
  #print (data)
  
  if (!is.null(cl <- attr(Terms, "dataClasses")))
    .checkMFClasses(cl, data, TRUE)
  
  #treatment <- data$`(treatment)`
  n <- nrow(data)
  Y <- model.response(data)
  where <- est.causalTree(object, causalTree.matrix(data))
  #wt <- model.weights(data)
  
  
  
  ## initialize the yval and dev
  yval <- rep(NA, num)
  dev <- rep(NA, num)
  count <- rep(0, num)
  
  # for loop to insert values:
  for (i in 1:n) {
    node_id <- where[i]
    while (node_id > 0) {
      index <- which(nodes == node_id)
      yval[index] <- sum(yval[index], Y[i], na.rm = T)
      dev[index] <- sum(dev[index], Y[i] * Y[i], na.rm = T)
      count[index] <- count[index] + 1
      node_id <- floor(node_id / 2)
    }
  }
  
  # get the final results for each：
  #causal_effect <- y1 /count1 - y0/count0
  #deviance <- dev1 - y1^2/count1 + dev0 - y0^2/count0
  yval <- yval / count
  dev <- dev - count * yval^2 
  
  na_yval_indx = which(is.na(yval))
  if (length(na_yval_indx) > 0){
    for (i in 1:length(na_yval_indx)) {
      tmp = na_yval_indx[i]
      node_idx = nodes[tmp]
      while (is.na(yval[tmp])) {
        node_idx = floor(node_idx / 2)
        tmp = which(nodes == node_idx)
      }
      yval[na_yval_indx[i]] = yval[tmp]
    }
  }
  
  na_dev_indx = which(is.na(dev))
  if (length(na_dev_indx) > 0) {
    for (i in 1:length(na_dev_indx )) {
      tmp = na_dev_indx[i]
      node_idx = nodes[tmp]
      while (is.na(dev[tmp])) {
        node_idx = floor(node_idx / 2)
        tmp = which(nodes == node_idx)
      }
      dev[na_dev_indx[i]] = dev[tmp]
    }
  }
  
  new_object = object
  new_object$frame$dev <- dev
  new_object$frame$yval <- yval
  new_object$frame$n <- count
  
  # here we set weights as 1 for all data points
  new_object$frame$wt <- count  
  return (new_object)
}




# consider honest tree for TOT:
honestTOT <- function(object, data, treatment, p = 0.5, na.action = na.causalTree) {

  warning( "The function honestTOT is deprecated, and will be removed. Please use refit.causalTree instead." )

  # p is the propensity socre:
  #if (!inherits(object, "rpart")) stop("Not a legitimate \"rpart\" object")
  
  if (!("treatment" %in% names(data))) stop("Function requires treatment column to be named 'treatment' in dataframe")
  
  nodes <- as.numeric(row.names(object$frame))
  num <- length(nodes)
  Terms <- object$terms
  treatment <- treatment
  data <- model.frame(Terms, data, na.action = na.action, treatment = treatment, 
                      xlev = attr(object, "xlevels"))
  #print (data)
  
  if (!is.null(cl <- attr(Terms, "dataClasses")))
    .checkMFClasses(cl, data, TRUE)
  
  treatment <- data$`(treatment)`
  n <- nrow(data)
  Y <- model.response(data) / ( p - 1 + treatment)
  
  where <- est.causalTree(object, causalTree.matrix(data))
  
  ## initialize the yval and dev
  yval <- rep(NA, num)
  dev <- rep(NA, num)
  count <- rep(0, num)
  
  # for loop to insert values:
  for (i in 1:n) {
    node_id <- where[i]
    while (node_id > 0) {
      index <- which(nodes == node_id)
      yval[index] <- sum(yval[index], Y[i], na.rm = T)
      dev[index] <- sum(dev[index], Y[i] * Y[i], na.rm = T)
      count[index] <- count[index] + 1
      node_id <- floor(node_id / 2)
    }
  }
  
  # get the final results for each：
  #causal_effect <- y1 /count1 - y0/count0
  #deviance <- dev1 - y1^2/count1 + dev0 - y0^2/count0
  yval <- yval / count
  dev <- dev - count * yval^2 
  
  ### now replace the NA values:
  # now look for the na values and repalce them with parent values:
  na_yval_indx = which(is.na(yval))
  if (length(na_yval_indx) > 0){
    for (i in 1:length(na_yval_indx)) {
      tmp = na_yval_indx[i]
      node_idx = nodes[tmp]
      while (is.na(yval[tmp])) {
        node_idx = floor(node_idx / 2)
        tmp = which(nodes == node_idx)
      }
      yval[na_yval_indx[i]] = yval[tmp]
    }
  }
  
  na_dev_indx = which(is.na(dev))
  if (length(na_dev_indx) > 0) {
    for (i in 1:length(na_dev_indx )) {
      tmp = na_dev_indx[i]
      node_idx = nodes[tmp]
      while (is.na(dev[tmp])) {
        node_idx = floor(node_idx / 2)
        tmp = which(nodes == node_idx)
      }
      dev[na_dev_indx[i]] = dev[tmp]
    }
  }
  
  new_object = object
  new_object$frame$dev <- dev
  new_object$frame$yval <- yval
  new_object$frame$n <- count
  
  # here we set weights as 1 for all data points
  new_object$frame$wt <- count
  
  return (new_object)
}


## honest causal tree models:
honestCTree <- function(object, data, treatment, na.action = na.causalTree) {

  warning( "The function honestCTree is deprecated, and will be removed. Please use refit.causalTree instead." )

  #if (!inherits(object, "rpart")) stop("Not a legitimate \"rpart\" object")
  if (!("treatment" %in% names(data))) stop("Function requires treatment column to be named 'treatment' in dataframe")
  
  nodes <- as.numeric(row.names(object$frame))
  num <- length(nodes)
  Terms <- object$terms
  #data <- model.frame(Terms, data, na.action = na.action,
  treatment <- treatment
  data <- model.frame(Terms, data, na.action = na.action, treatment = treatment, 
                      xlev = attr(object, "xlevels"))
  #print (data)
  
  if (!is.null(cl <- attr(Terms, "dataClasses")))
    .checkMFClasses(cl, data, TRUE)
  
  treatment <- data$`(treatment)`
  n <- nrow(data)
  Y <- model.response(data)
  where <- est.causalTree(object, causalTree.matrix(data))
  #wt <- model.weights(data)
  
  ## begin to compute the yval and dev:
  # initialize:
  # where num is number of nodes
  y1 <- rep(NA, num)
  y0 <- rep(NA, num)
  dev1 <- rep(NA,num)
  dev0 <- rep(NA,num)
  count1 <- rep(0, num)
  count0 <- rep(0, num)
  
  # for loop to insert values:
  for (i in 1:n) {
    node_id <- where[i]
    while (node_id > 0) {
      index <- which(nodes == node_id)
      if (treatment[i] == 1) {
        y1[index] <- sum(y1[index], Y[i], na.rm = T)
        dev1[index] <- sum(dev1[index], Y[i] * Y[i], na.rm = T)
        count1[index] <- count1[index] + 1
      } else {
        y0[index] <- sum(y0[index], Y[i], na.rm = T)
        dev0[index] <- sum(dev0[index], Y[i] * Y[i], na.rm = T)
        count0[index] <- count0[index] + 1
      }
      node_id <- floor(node_id / 2)
    }
  }
  
  
  
  # get the final results for each：
  causal_effect <- y1 /count1 - y0/count0
  deviance <- dev1 - y1^2/count1 + dev0 - y0^2/count0
  
  # now look for the na values and repalce them with parent values:
  na_causal_indx = which(is.na(causal_effect))
  if (length(na_causal_indx) > 0) {
    for (i in 1:length(na_causal_index)) {
      tmp = na_causal_indx[i]
      node_idx = nodes[tmp]
      while (is.na(causal_effect[tmp])) {
        node_idx = floor(node_idx /2)
        tmp = which(nodes == node_idx)
      }
      causal_effect[na_causal_indx[i]] = causal_effect[tmp]
    }
  }
  
  # similar for deviance:
  na_dev_indx = which(is.na(deviance))
  if (length(na_dev_indx) < 0) {
    for (i in 1:length(na_dev_indx)) {
      tmp = na_dev_indx[i]
      node_idx = nodes[tmp]
      while (is.na(deviance[tmp])) {
        node_idx = floor(node_idx /2)
        tmp = which(nodes == node_idx)
      }
      deviance[na_dev_indx[i]] = deviance[tmp]    
    }
  }
  
  ## above finish athe recursive replacement:
  new_object = object
  new_object$frame$dev <- deviance
  new_object$frame$yval <- causal_effect 
  new_object$frame$n <- count1 + count0
  
  # here we set weights as 1 for all data points
  new_object$frame$wt <- count1 + count0
  return (new_object)
}
