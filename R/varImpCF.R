
#Function to calculate variable importance scores
#Takes a causal forest object, and iterates through the variables in `vars`, permutng them and then estimating the difference in te (out-of-bag) predictions from the resulting forests

varImpCF <- function(f, test.data, parallel = TRUE, cores = ifelse(parallel == TRUE, detectCores()-1, NULL), scrambletimes = 5, vars = NULL, verbose = FALSE){
#f <- sfor
#test.data <- X
  registerDoMC(cores)
  pob <- predictOOB(f, test.data)
  if (is.null(vars)){
    vars <- colnames(f$x)
  }
  if (parallel == TRUE){
    `%fun%` <- `%dopar%`
    registerDoMC(cores)
  } else {
    `%fun%` <- `%do%`
  }
  imp <- foreach(i = 1:length(vars), .combine = rbind) %fun% {
    MSD <- foreach(j = 1:scrambletimes, .combine = c) %do% {  
      ds <- test.data
      ds[,vars[i]] <- sample(ds[,vars[i]])
      ps <- predictOOB(f, ds)                
      mean((pob-ps)^2)
    }
    if (verbose == TRUE) {print(i)}
    c(mean(MSD), sd(MSD))
  }
  imp <- as.data.frame(imp)
  imp <- cbind(vars, imp)
  colnames(imp) <- c('varname', 'importance', 'sd_imp')
  imp <- imp[with(imp, order(importance, decreasing = TRUE)),]
  return(imp)
}


