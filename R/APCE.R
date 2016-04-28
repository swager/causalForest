
#Functions for calculating average partial causal effects, and plotting them.


load_obj <- function(f)
{
    env <- new.env()
    nm <- load(f, env)[1]
    env[[nm]]
}

#Function to pre-predict for dtauhat_dx, so that I don't need to feed the whole causal tree to that function
prepred <- function(m, dat, vars, ngrid = rep(5, length(vars)), vals = NULL, cond = NULL, verbose = T){
  m <- load_obj(m)
  if (verbose == T){print('loaded CF')}
  if (is.null(vals)){
    vals <- list()
    for (i in 1:length(vars)){
      if (class(dat[,vars[i]]) %in% c('integer','numeric')){
        if (length(unique(dat[,vars[i]]))<ngrid[i]){
          vals[[vars[i]]] <- sort(unique(dat[,vars[i]]))
        } else {
          vals[[vars[i]]] <- unique(quantile(dat[,vars[i]], probs = seq(.05, .95, length.out=ngrid[i])))
        }
      } else {
        vals[[vars[i]]] <- levels(dat[,vars[i]])
      }
    }
  }
  q <- expand.grid(vals)
  if (verbose == T){print(q)}
  if (!is.null(cond)){
    for (i in 1:length(cond)){
      dat[,names(cond)[i]] <- cond[[i]]
    }
  }
  predlist <- list()
  for (i in 1:nrow(q)){
    dt <- dat
    dt[,vars] <-  q[i,1:(length(vars))]
    if (verbose == T){print(i)}
    predlist[[i]] <- predict(m, newdata = dt, predict.all = TRUE)
  }
  rm(m)
  gc()
  return(predlist)
}



dtauhat_dx <- function(m = NULL, inbag = NULL, prepreds = NULL, dat, vars, ngrid = rep(5, length(vars)), plot_me = TRUE, clusters = NULL, parallel = FALSE, cores = NULL, verbose = FALSE, cond = NULL, weights = rep(1, nrow(dat)), useHardDisk = FALSE, vals = NULL){
  if (is.null(vals)){
    vals <- list()
    for (i in 1:length(vars)){
      if (class(dat[,vars[i]]) %in% c('integer','numeric')){
        if (length(unique(dat[,vars[i]]))<ngrid[i]){
          vals[[vars[i]]] <- sort(unique(dat[,vars[i]]))
        } else {
          vals[[vars[i]]] <- unique(quantile(dat[,vars[i]], probs = seq(.05, .95, length.out=ngrid[i])))
        }
      } else {
        vals[[vars[i]]] <- levels(dat[,vars[i]])
      }
    }
  }
  q <- expand.grid(vals)
  colnames(q) <- vars
  q$ciup <- q$mu <- q$cidown <- NA
  if (is.null(clusters)){
    clusters <- 1:nrow(dat)
  }
  if (!is.null(cond)){
    for (i in 1:length(cond)){
      dat[,names(cond)[i]] <- cond[[i]]
    }
  }
  for (i in 1:nrow(q)){
    if (verbose == TRUE){print(paste('calc for val',i))}
    dt <- dat
    dt[,vars] <-  q[i,1:(length(vars))]
    if (!is.null(prepreds)){
      p <- JAB(taux = prepreds[[i]], inbag = inbag, test.data = dt, clusters = clusters
        , parallel = parallel, cores = cores, verbose = verbose, weights = weights, useHardDisk = useHardDisk)
    } else {
      p <- JAB(f = m, test.data = dt, clusters = clusters
        , parallel = parallel, cores = cores, verbose = verbose, weights = weights, useHardDisk = useHardDisk)
    }
    se <- sqrt(sum(p[[2]])/(length(p[[1]])^2))
    mu <- mean(p[[1]])
    q[i,c('cidown','mu','ciup')] <- c(mu-se*1.96, mu, mu+se*1.96)
    if (verbose == T){print(q)}
  }
  if (plot_me == TRUE){
    plot_dtaudx(q, dat)
  }
  return(q)
}


plot_dtaudx <- function(q, dat, plotmat = c(1,3), main = '', boxsize = 2.5, wsize = 5, xlab = NULL){
  if (ncol(q)>4){#if bivariate
    numUnique <- apply(q[,1:2],2,function(x){length(unique(x))})
    mumat <- reshape(as.data.frame(q[,c(1,2,4)]), direction = 'wide', idvar = colnames(q)[1], timevar = colnames(q)[2])[,-1]
    upmat <- reshape(as.data.frame(q[,c(1,2,5)]), direction = 'wide', idvar = colnames(q)[1], timevar = colnames(q)[2])[,-1]
    dnmat <- reshape(as.data.frame(q[,c(1,2,3)]), direction = 'wide', idvar = colnames(q)[1], timevar = colnames(q)[2])[,-1]
    par(mfrow = plotmat)    
    contour(x = sort(unique(q[,1])), y = sort(unique(q[,2])), z = as.matrix(mumat)
      , xlab = colnames(q)[1], ylab = colnames(q)[2], main = expression(paste('E[',hat(tau),']'))
      , xlim = range(q[1]), ylim = range(q[2]))
    contour(x = sort(unique(q[,1])), y = sort(unique(q[,2])), z = as.matrix(upmat)
      , xlab = colnames(q)[1], ylab = colnames(q)[2], main = 'Upper CI', xlim = range(q[1]), ylim = range(q[2]))
    contour(x = sort(unique(q[,1])), y = sort(unique(q[,2])), z = as.matrix(dnmat)
      , xlab = colnames(q)[1], ylab = colnames(q)[2], main = 'Lower CI', xlim = range(q[1]), ylim = range(q[2]))
  } else {#univariate
    vals <- q[,1]
    vars <- colnames(q)[1]
    q <- q[,-1]
    if (class(dat[,vars]) %ni% c('integer','numeric') | length(vals) < 4){
      plot(NULL, xlim = c(.5,(length(vals)+.5)), ylim = range(q), xlab = ifelse(is.null(xlab), vars, xlab), ylab = expression(paste(tau,'(X)')), xaxt = 'n', main = main)
      axis(1, labels = as.character(vals), at = 1:length(vals))            
      segments(x0 = 1:length(vals), y0 = q[,1], y1 = q[,3])
      points(x = 1:length(vals), y = q[,2], cex = boxsize, pch = 15, col = 'red')
      points(x = 1:length(vals), y = q[,3], pch = '-', cex = wsize)
      points(x = 1:length(vals), y = q[,1], pch = '-', cex = wsize)
    } else {
      plot(NULL, xlim = range(vals), ylim = range(q), xlab = ifelse(is.null(xlab), vars, xlab), ylab = expression(paste(tau,'(X)')), main = main)
      polygon(x = c(vals,rev(vals)), y = c(q[,1], rev(q[,3])), col = rgb(1,0,0,.3), border = NA)
      points(vals, q[,2], pch = 19)
      lines(vals, q[,2], lty=2)
    }
  }
  return(q)
}



