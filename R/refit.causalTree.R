# Refit causal tree; used for getting honest trees
#
# object: the original causalTree object
# newx: the fresh training data features, for honest estimation
# newy: outcomes for new data
# treatment: treatment assignments for new data
# propensity: If NULL, we estimate a CT. If provided, estimate a TOT.

refit.causalTree <- function(object, newx, newy, treatment, na.action = na.causalTree, propensity = NULL) {

  if (!inherits(object, "causalTree")) stop("Not a legitimate \"causalTree\" object")
  
  if(length(newy) != length(treatment)) stop("The vectors newy and treatment must have the same length.")
  
  mode = "CT"
  if(!is.null(propensity)) {
  	mode = "TOT"
  	if(min(propensity) < 0 | max(propensity) > 1 | length(propensity) != nrow(newx)) {
  		stop("Invalid propensities")
  	}
  }
  
  nodes <- as.numeric(row.names(object$frame))
  num <- length(nodes)
  
  Terms <- object$terms
  data <- model.frame(Terms, newx, na.action = na.action, xlev = attr(object, "xlevels"))
  if (!is.null(cl <- attr(Terms, "dataClasses"))) {
    .checkMFClasses(cl, newx, TRUE)
  }
  
  where <- est.causalTree(object, causalTree.matrix(data))
  
  if (mode == "CT") {
    
    Y <- newy
    
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
    for (i in 1:nrow(newx)) {
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
    causal_effect <- y1/count1 - y0/count0
    deviance <- dev1 - y1^2/count1 + dev0 - y0^2/count0
    count <- count0 + count1
    
  } else if (mode == "TOT") {
  	
  	Y <- newy / ( propensity - 1 + treatment)
   
    yval <- rep(NA, num)
    dev <- rep(NA, num)
    count <- rep(0, num)
  
    # for loop to insert values:
    for (i in 1:nrow(newx)) {
      node_id <- where[i]
      while (node_id > 0) {
        index <- which(nodes == node_id)
        yval[index] <- sum(yval[index], Y[i], na.rm = T)
        dev[index] <- sum(dev[index], Y[i] * Y[i], na.rm = T)
        count[index] <- count[index] + 1
        node_id <- floor(node_id / 2)
      }
    }
   
    causal_effect <- yval / count
    deviance <- dev - count * yval^2
  
  } else {
  	stop("Invalid mode")
  }
  
  # now look for the na values and repalce them with parent values:
  na_causal_index = which(is.na(causal_effect))
  if (length(na_causal_index) > 0) {
    for (i in 1:length(na_causal_index)) {
      tmp = na_causal_index[i]
      node_idx = nodes[tmp]
      while (is.na(causal_effect[tmp])) {
        node_idx = floor(node_idx /2)
        tmp = which(nodes == node_idx)
      }
      causal_effect[na_causal_index[i]] = causal_effect[tmp]
    }
  }
  
  # similar for deviance:
  na_dev_index = which(is.na(deviance))
  if (length(na_dev_index) < 0) {
    for (i in 1:length(na_dev_index)) {
      tmp = na_dev_index[i]
      node_idx = nodes[tmp]
      while (is.na(deviance[tmp])) {
        node_idx = floor(node_idx /2)
        tmp = which(nodes == node_idx)
      }
      deviance[na_dev_index[i]] = deviance[tmp]    
    }
  }
  
  ## above finish athe recursive replacement:
  new_object = object
  new_object$frame$dev <- deviance
  new_object$frame$yval <- causal_effect 
  new_object$frame$n <- count
  
  # here we set weights as 1 for all data points
  new_object$frame$wt <- count
  return (new_object)
}
