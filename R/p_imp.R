
#Function calculates p-values for variable importance metrics, following Altmann et al 2010

load_obj <- function(f)
{
    env <- new.env()
    nm <- load(f, env)[1]
    env[[nm]]
}

p_imp <- function(x, y, w, weights, vars, parallel = TRUE, cores = ifelse(parallel == TRUE, detectCores()-1, NULL), n_iter = 100, ntree = 100, clustvar, VI, srate = round(nrow(x)*.45), tempsave = FALSE, usetempsave = FALSE, verbose = TRUE){
#w = x$W
#y=x$y
#x=x[grepl('V',colnames(x))]
#weights = rep(1, nrow(x))
#vars = as.character(VI_comp$varname[1:5])
#clustvar = 1:nrow(x)
#VI=VI
#n_iter = 100
#cores = detectCores()
#parallel=TRUE
#tempsave = TRUE
#verbose = TRUE
#usetempsave = FALSE
#srate = round(nrow(x)*.45)
#ntree = 500
  if (parallel == TRUE){
    `%fun%` <- `%dopar%`
    registerDoMC(cores)
  } else {
    `%fun%` <- `%do%`
  }
  if (tempsave == TRUE){
    if (usetempsave == TRUE) {
      ts_files <- list.files(,path = paste0(getwd(),'/ts'))
    } else {
      system('mkdir ts')
    }
  }
  if (usetempsave == TRUE){
      done <- gsub('ts','',ts_files)
      todo <- (1:n_iter)[1:n_iter %ni% done]
  } else {
    todo <- 1:n_iter
  }
  sdist <- foreach(i = todo) %fun% {
    set.seed(i)
    if (verbose == TRUE){print(paste('starting null tree',i))}
    pt <- proc.time()
    sfor <- causalForest(
        X = x
      , Y = sample(y)
      , W = w
      , num.trees = ntree
      , weights = weights
      , clustvar = clustvar
      , verbose = FALSE
      , sample.size = srate
      , cores = 1
    )
    if (verbose == TRUE){print(paste('fit null CF',i))}
    sVI <- varImpCF(sfor, x, scrambletimes = 1, parallel = FALSE, cores = 1, vars = vars)
    if (tempsave == TRUE){
      save(sVI, file = paste0(getwd(),'/ts/','ts',i))
      if (verbose == TRUE){
        print(paste('saved temp save file at', getwd()))
      }
    }
    ft <- proc.time()-pt
    if (verbose == TRUE){print(paste('finished null sim',i,'in',ft[3]))}
    return(sVI[,1:2])
  }
  #combine the saved ones into the new ones
  if (usetempsave == TRUE){
    if ('sdist' %ni% ls()){sdist <- vector('list',n_iter)}
    for (i in 1:length(ts_files)){
      num <- as.numeric(gsub('ts','',ts_files[i]))
      sdist[[num]] <- load_obj(paste0(getwd(),'/ts/ts',i))
    }
    print('WARNING!  You should probably delete the directory of temp save files now.  The code doesnt do it automatically because it seems risky.')
  }
  nullImpMat <- foreach(i = 1:n_iter, .combine = cbind) %fun%{
    x <- sdist[[i]]
    sorted <- x[with(x, order(varname)),]
    sorted[,2]
  }
  names <- as.character(sort(sdist[[1]]$varname))
  mdn <- ci90 <- ci95 <- pvals <- rep(NA,length(names))
  for (i in 1:length(names)){
    imp <- VI[VI$varname == names[i],'importance']
    q <- quantile(nullImpMat[i,], probs = c(.5,.9,.95))
    pvals[i] <- 1-sum(imp>nullImpMat[i,])/n_iter
    mdn[i] <- q[1]
    ci90[i] <- q[2]
    ci95[i] <- q[3]
  }
  return(as.data.frame(cbind(names, pvals, mdn, ci90, ci95)))
}


