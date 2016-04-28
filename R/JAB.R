

#Function computes the jackknife-after-bootstrap estimate of the covariance of tau(x)

load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

rep.col<-function(x,n){
   matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

#This function can take a CF object OR can take the combination of a `inbag` matrix and a prediction (with all=TRUE).  The latter has a lower memory footprint.

JAB <- function(f= NULL, inbag = NULL, taux = NULL, test.data, clusters = 1:nrow(test.data), parallel = TRUE, cores = ifelse(parallel == TRUE, detectCores()-1, NULL), verbose = FALSE, weights = rep(1, nrow(test.data)), useHardDisk = FALSE){
#f <- CF
#test.data <- x
#clusters <- 1:nrow(x)
####clusters <- 1:nrow(X)
#parallel = FALSE
#verbose = TRUE
#weights = rep(1, nrow(test.data))
  if (!is.null(f)){
    taux <- predict(f, newdata = test.data, predict.all = TRUE)
    inbag <- f$inbag
    rm(f)
    gc()
  } 
  if (parallel == TRUE){
    `%fun%` <- `%dopar%`
    registerDoMC(cores)
    if (useHardDisk == TRUE){
      system('mkdir JABtempmats')
#print(length(unique(clusters)))
      foreach(i = 1:length(unique(clusters))) %fun% {      
        if (verbose == TRUE) {print(paste("Calc for cluster",i))}
        cli <- which(clusters == unique(clusters)[[i]])
        inbag_cli <- which(colMeans(inbag[cli,,drop = F]) == 1)
        tauxi <- as.matrix(rowMeans(taux$individual[,inbag_cli]))
        wi <- weights[cli]
        covmat <- matrix(rep(0,length(tauxi)^2),length(tauxi))
        for (j in 1:length(tauxi)){
          covmat[,j] <- (taux$aggregate[j] - tauxi[j])*(taux$aggregate - tauxi)
        }
        covmat <- covmat*sum(wi)
        save(covmat, file = paste0(getwd(),'/JABtempmats/covmat',i))
        return('nothing to return')
      }
      fi <- list.files(,path = paste0(getwd(),'/JABtempmats'))
      varis <- load_obj(paste0(getwd(),'/JABtempmats/',fi[1]))
#print(dim(varis))
      for (i in 2:length(fi)){
        mat <- load_obj(paste0(getwd(),'/JABtempmats/',fi[i]))
#print(i)
        varis <- varis+mat
      }
      system('rm -r JABtempmats')
    } else {#If all of the matrices are to be stored in memory
      varis <- foreach(i = 1:length(unique(clusters)), .combine = '+') %fun% {
  #i=1
        if (verbose == TRUE) {print(paste("Calc for cluster",i))}
        cli <- which(clusters == unique(clusters)[[i]])
        inbag_cli <- which(colMeans(inbag[cli,,drop = F]) == 1)
        tauxi <- as.matrix(rowMeans(taux$individual[,inbag_cli]))
        wi <- weights[cli]
        covmat <- matrix(rep(0,length(tauxi)^2),length(tauxi))
        for (j in 1:length(tauxi)){
          covmat[,j] <- (taux$aggregate[j] - tauxi[j])*(taux$aggregate - tauxi)
        }
        covmat*sum(wi)
      }
    }
  } else {
    `%fun%` <- `%do%`
    varis <- matrix(rep(0, length(taux$aggregate)^2), length(taux$aggregate))
    for (i in 1:length(unique(clusters))){
      if (verbose == TRUE) {print(paste("Calc for cluster",i))}
      cli <- which(clusters == unique(clusters)[[i]])
      inbag_cli <- which(colMeans(inbag[cli,,drop = F]) == 1)
      tauxi <- as.matrix(rowMeans(taux$individual[,inbag_cli]))
      wi <- weights[cli]
      covmat <- matrix(rep(0,length(tauxi)^2),length(tauxi))
      for (j in 1:length(tauxi)){
        covmat[,j] <- (taux$aggregate[j] - tauxi[j])*(taux$aggregate - tauxi)
      }
      covmat <- covmat*sum(wi)
      varis <- varis+covmat
    }
  }
  finSampCorr <- (nrow(test.data)-1)/sum(weights)*length(unique(clusters))/(length(unique(clusters))-1)
  return(list(tauhat = taux$aggregate, vcv = varis*finSampCorr))
}

##This is an old version without the "useHardDisk" option
#JAB <- function(f= NULL, inbag = NULL, taux = NULL, test.data, clusters = 1:nrow(test.data), parallel = TRUE, cores = ifelse(parallel == TRUE, detectCores()-1, NULL), verbose = FALSE, weights = rep(1, nrow(test.data))){
##f <- CF
##test.data <- x
##clusters <- 1:nrow(x)
#####clusters <- 1:nrow(X)
##parallel = FALSE
##verbose = TRUE
##weights = rep(1, nrow(test.data))
#  if (!is.null(f)){
#    taux <- predict(f, newdata = test.data, predict.all = TRUE)
#    inbag <- f$inbag
#    rm(f)
#    gc()
#  } 
#  if (parallel == TRUE){
#    `%fun%` <- `%dopar%`
#    registerDoMC(cores)
#    varis <- foreach(i = 1:length(unique(clusters)), .combine = '+') %fun% {
##i=1
#      if (verbose == TRUE) {print(paste("Calc for cluster",i))}
#      cli <- which(clusters == unique(clusters)[[i]])
#      inbag_cli <- which(colMeans(inbag[cli,,drop = F]) == 1)
#      tauxi <- as.matrix(rowMeans(taux$individual[,inbag_cli]))
#      wi <- weights[cli]
#      covmat <- matrix(rep(0,length(tauxi)^2),length(tauxi))
#      for (j in 1:length(tauxi)){
#        covmat[,j] <- (taux$aggregate[j] - tauxi[j])*(taux$aggregate - tauxi)
#      }
#      covmat*sum(wi)
#    }
#  } else {
#    `%fun%` <- `%do%`
#    varis <- matrix(rep(0, length(taux$aggregate)^2), length(taux$aggregate))
#    for (i in 1:length(unique(clusters))){
#      if (verbose == TRUE) {print(paste("Calc for cluster",i))}
#      cli <- which(clusters == unique(clusters)[[i]])
#      inbag_cli <- which(colMeans(inbag[cli,,drop = F]) == 1)
#      tauxi <- as.matrix(rowMeans(taux$individual[,inbag_cli]))
#      wi <- weights[cli]
#      covmat <- matrix(rep(0,length(tauxi)^2),length(tauxi))
#      for (j in 1:length(tauxi)){
#        covmat[,j] <- (taux$aggregate[j] - tauxi[j])*(taux$aggregate - tauxi)
#      }
#      covmat <- covmat*sum(wi)
#      varis <- varis+covmat
#    }
#  }
#  finSampCorr <- (nrow(test.data)-1)/sum(weights)*length(unique(clusters))/(length(unique(clusters))-1)
#  return(list(tauhat = taux$aggregate, vcv = varis*finSampCorr))
#}

#    varis <- foreach(i = 1:length(unique(clusters)), .combine = '+') %dopar% {
#      if (verbose == TRUE) {print(paste("Calc for cluster",i))}
#      cli <- which(clusters == unique(clusters)[[i]])
#      inbag_cli <- which(colMeans(inbag[cli,,drop = F]) == 1)
#      tauxi <- as.matrix(rowMeans(taux$individual[,inbag_cli]))
#      covmat <- matrix(rep(0,length(tauxi)^2),length(tauxi))
#      for (j in 1:length(tauxi)){
#        covmat[,j] <- (taux$aggregate[j] - tauxi[j])*(taux$aggregate - tauxi)
#      }
#      covmat*length(cli)
#    }
#  } else {
#      covmat <- matrix(rep(0,length(tauxi)^2),length(tauxi))
#      for (j in 1:length(tauxi)){
#        covmat[,j] <- (taux$aggregate[j] - tauxi[j])*(taux$aggregate - tauxi)
#      }


#    }
#  }

#This is an old version that didn't take weights
#JAB <- function(f= NULL, inbag = NULL, taux = NULL, test.data, clusters = 1:nrow(test.data), parallel = TRUE, cores = ifelse(parallel == TRUE, detectCores()-1, NULL), verbose = FALSE){
#  if (!is.null(f)){
#    taux <- predict(f, newdata = test.data, predict.all = TRUE)
#    inbag <- f$inbag
#    rm(f)
#    gc()
#  } 
#  if (parallel == TRUE){
#    `%fun%` <- `%dopar%`
#    registerDoMC(cores)
#    varis <- foreach(i = 1:length(unique(clusters)), .combine = '+') %fun% {
#      if (verbose == TRUE) {print(paste("Calc for cluster",i))}
#      cli <- which(clusters == unique(clusters)[[i]])
#      inbag_cli <- which(colMeans(inbag[cli,,drop = F]) == 1)
#      tauxi <- as.matrix(rowMeans(taux$individual[,inbag_cli]))
#      covmat <- matrix(rep(0,length(tauxi)^2),length(tauxi))
#      for (j in 1:length(tauxi)){
#        covmat[,j] <- (taux$aggregate[j] - tauxi[j])*(taux$aggregate - tauxi)
#      }
#      covmat*length(cli)
#    }
#  } else {
#    `%fun%` <- `%do%`
#    varis <- matrix(rep(0, length(taux$aggregate)^2), length(taux$aggregate))
#    for (i in 1:length(unique(clusters))){
#      if (verbose == TRUE) {print(paste("Calc for cluster",i))}
#      cli <- which(clusters == unique(clusters)[[i]])
#      inbag_cli <- which(colMeans(inbag[cli,,drop = F]) == 1)
#      tauxi <- as.matrix(rowMeans(taux$individual[,inbag_cli]))
#      covmat <- matrix(rep(0,length(tauxi)^2),length(tauxi))
#      for (j in 1:length(tauxi)){
#        covmat[,j] <- (taux$aggregate[j] - tauxi[j])*(taux$aggregate - tauxi)
#      }
#      covmat <- covmat*length(cli)
#      varis <- varis+covmat
#    }
#  }
#  finSampCorr <- (nrow(test.data)-1)/nrow(test.data)*length(unique(clusters))/(length(unique(clusters))-1)
#  return(list(tauhat = taux$aggregate, vcv = varis*finSampCorr))
#}



####This is an old version that didn't return the covariance matrix
#JAB <- function(f, test.data, clusters, parallel = TRUE, cores = ifelse(parallel == TRUE, detectCores()-1, NULL), verbose = FALSE){
##f <- CF
##rm(comp_forest)
##test.data <- x
##clusters <- 1:nrow(x)
####clusters <- 1:nrow(X)
##parallel = FALSE
##verbose = TRUE
#  taux <- predict(f, newdata = test.data, predict.all = TRUE)
#  inbag <- f$inbag
#  rm(f)
#  gc()
#  if (parallel == TRUE){
#  registerDoMC(cores)
#    varis <- foreach(i = 1:length(unique(clusters)), .combine = cbind) %dopar% {
#      if (verbose == TRUE) {print(paste("Calc for cluster",i))}
#      cli <- which(clusters == unique(clusters)[[i]])
#      inbag_cli <- which(colMeans(inbag[cli,, drop = FALSE]) == 1)
#      tauxi <- rowMeans(taux$individual[,inbag_cli])
#      varcli <- (taux$aggregate - tauxi)^2
#      matrix(rep(varcli,length(cli)), length(taux$aggregate))
#    }
#  } else {
#    varis <- foreach(i = 1:length(unique(clusters)), .combine = cbind) %do% {
#      if (verbose == TRUE) {print(paste("Calc for cluster",i))}
#      cli <- which(clusters == unique(clusters)[[i]])
#      inbag_cli <- which(colMeans(inbag[cli,, drop = FALSE]) == 1)
#      tauxi <- rowMeans(taux$individual[,inbag_cli])
#      varcli <- (taux$aggregate - tauxi)^2
#      matrix(rep(varcli,length(cli)), length(taux$aggregate))
#    }
#  }
#  finSampCorr <- (length(test.data)-1)/length(test.data)*length(unique(clusters))/(length(unique(clusters))-1)
#  return(data.frame(y.hat = taux$aggregate, var.hat = rowSums(varis)*finSampCorr))
#}

