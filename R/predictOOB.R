
#Out-of-bag predictions for causal forests

predictOOB <- function(f, test.data){
#f <- CF
#test.data = x
  pred <- f$inbag*NA
  for (i in 1:ncol(f$inbag)){
#i=1
    oob <- test.data[f$inbag[,i] == 0,, drop = FALSE]
    colnames(oob) <- paste0('X.', colnames(oob))
    pred[f$inbag[,i] == 0,i] <- predict(f$tree[[i]], oob)
  }
  p <- rowMeans(pred, na.rm=TRUE)
  if (any(is.na(p))){stop('Some NA OOB predictions.  Ntree too small or s too big.')}
  return(p)
}


