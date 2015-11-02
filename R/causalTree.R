 #
#  The recursive partitioning function, for R
#

# TODO: add propensity score p in causalTree.R
# TODO: add cv.option in causalTree.R
causalTree <-
 #   function(formula, data, weights, subset, na.action = na.causalTree, method,
  function(formula, data, weights, treatment, subset, na.action = na.causalTree, method = "anova", 
           split.option, cv.option, minsize = 2L, model = FALSE, x = FALSE, y = TRUE, parms, p, control, alpha = 0.5, cost, ...) 
    # p := propensity score
    # split.option := CT or TOT, splitting rule
    # cv.option := cross validation option, TOT or matching
    # minsize := minimum leaf size for both control and treated cases, the default is 2.
{   
    mtemp <- model
    ## begin to call:
    Call <- match.call()
    if (is.data.frame(model)) {
        m <- model
        # set the temp model:
        #mtemp <- model
        model <- FALSE
    } else {
        indx <- match(c("formula", "data", "weights", "subset"),
                      names(Call), nomatch = 0L)
        if (indx[1] == 0L) stop("a 'formula' argument is required")
        temp <- Call[c(1L, indx)]      # only keep the arguments we wanted
        temp$na.action <- na.action    # This one has a default
        temp[[1L]] <- quote(stats::model.frame) # change the function called
        m <- eval.parent(temp)
    }
    
    Terms <- attr(m, "terms")
    if (any(attr(Terms, "order") > 1L))
      stop("Trees cannot handle interaction terms")
    
    Y <- model.response(m)
    wt <- model.weights(m)
    if (any(wt < 0)) stop("negative weights not allowed")
    if (!length(wt)) wt <- rep(1, nrow(m))
    offset <- model.offset(m)
    X <- causalTree.matrix(m)
    nobs <- nrow(X)
    nvar <- ncol(X)
    
    # requirement for split.option:
    # set this temporarily:
    if (missing(split.option)) {
      warning("The default splitting method is CT.")
      split.option <- "CT"
    } else {
      split.option.num <- pmatch(split.option, c("CT", "TOT"))
      if(is.na(split.option.num)) stop("Invalid spliting option.")
    }
    
    #requirement for treatment status
    if (missing(treatment)) 
      stop("You should input the treatment status vector for data:
           1 represent treated and 0 represent control.")
    if (sum(treatment %in% c(0,1)) != nobs)
      stop("The treatment status should be 1 or 0 only: 1 represent treated and 0 represent control.")
    if (sum(treatment) == 0 || sum(treatment) == nobs)
      stop("The data only contains treated cases or control cases, please check 'treatment' again.") 
    
    ## use original rpart for TOT method:
    if (split.option == "TOT") {
      ## set the tr vector for TOT method:
      if (missing(p)) {
        p = sum(treatment) / nobs
      } else if (p > 1 || p < 0) {
        stop("Propensity score should be between 0 and 1.")
      }
      trp <- 1 / (p - 1 + treatment)
      
      formula.elements <- strsplit(deparse(temp$formula), "~")[[1]]
      if (length(formula.elements) != 2)
        stop("Please enter valid formula.")
      formula.left <- formula.elements[1]
      formula.right <- formula.elements[2]
      formula.left <- paste(formula.left, "*", "trp")
      myformula <- as.formula(paste(formula.left, "~", formula.right))

      ## TODO: deal with subset
      if (missing(subset)) {
        ans <- rpart(formula = myformula, data = data, weights = wt, 
                   na.action = na.rpart, method = method, model = mtemp, x = x, y = y, 
                   parms = parms, control = control, cost = cost, ...)
      } else {
        ans <- rpart(formula = myformula, data = data, weights = wt, subset = subset,
                     na.action = na.rpart, method = method, model = mtemp, x = x, y = y, 
                     parms = parms, control = control, cost = cost, ...)
      }
      class(ans) <- "causalTree"
      return (ans)
    }

    ## get back to causalTree for CT method:
    # here we add the variance of the predictors of X:
    xvar <- apply(X, 2, var)

    if (missing(method)) {
      method <- if (is.factor(Y) || is.character(Y)) "class"
      else if (inherits(Y, "Surv")) "exp"
      else if (is.matrix(Y)) "poisson"
	    else "anova"
    }
    
    # requirement for cv.option in CT method:
    if (missing(cv.option)) {
      warning("The default corss validation method is TOT.")
      cv.option <-"TOT"
    } else {
      cv.option.num <- pmatch(cv.option, c("TOT", "matching"))
      if(is.na(cv.option.num)) stop("Invalid cv option.")
    }
    
    if (is.list(method)) {
        ## User-written split methods
      mlist <- method
	    method <- "user"
        ## Set up C callback.  Assign the result to a variable to avoid
        ## garbage collection
	    init <- if (missing(parms)) mlist$init(Y, offset, wt = wt)
	            else mlist$init(Y, offset, parms, wt)
      keep <- causalTreecallback(mlist, nobs, init)

    	method.int <- 4L             # the fourth entry in func_table.h
      #	numresp <- init$numresp
      #	numy <- init$numy
    	parms <- init$parms
    } else { 
      # not user function
	    method.int <- pmatch(method, c("anova", "poisson", "class", "exp", "anova2"))
	    if (is.na(method.int)) stop("Invalid method")
	    method <- c("anova", "poisson", "class", "exp", "anova2")[method.int]
    	if (method.int == 4L) method.int <- 2L

        ## If this function is being retrieved from the causalTree package, then
        ##   preferentially "get" the init function from there.  But don't
        ##   lock in the causalTree package otherwise, so that we can still do
        ##   standalone debugging.
  
     # consider differenct CV method.
      if (cv.option == "TOT") {
        if (missing(p)) {
          #stop("For TOT cross-validation test, propensity socre is needed.")
          #warning("For TOT cv test in CT, the defualt propensity score is proportion of treated cases.")
          p = sum(treatment) / nobs
        } else if (p > 1 || p < 0) {
          stop("Propensity score should be between 0 and 1.")
        }
      } else {
      # cv.option = "matching"                           
      p = -1
      }
        
	    init <- {
        if (missing(parms))
            get(paste("causalTree", method, sep = "."),
                 envir = environment())(Y, offset, , wt) 
        else
            get(paste("causalTree", method, sep = "."),
                 envir = environment())(Y, offset, parms, wt)
  	  }
        ## avoid saving environment on fitted objects
        ns <- asNamespace("causalForest")
        if (!is.null(init$print)) environment(init$print) <- ns
        if (!is.null(init$summary)) environment(init$summary) <- ns
        if (!is.null(init$text)) environment(init$text) <- ns
    }

    Y <- init$y

    xlevels <- .getXlevels(Terms, m)
    cats <- rep(0L, ncol(X))
    if (!is.null(xlevels))
      cats[match(names(xlevels), colnames(X))] <-
            unlist(lapply(xlevels, length))

    ## We want to pass any ... args to causalTree.control, but not pass things
    ##  like "dats = mydata" where someone just made a typo.  The use of ...
    ##  is simply to allow things like "cp = 0.05" with easier typing
    
    extraArgs <- list(...)
    if (length(extraArgs)) {
      controlargs <- names(formals(causalTree.control)) # legal arg names
	    indx <- match(names(extraArgs), controlargs, nomatch = 0L)
	    if (any(indx == 0L))
            stop(gettextf("Argument %s not matched",
                          names(extraArgs)[indx == 0L]),
                 domain = NA)
    }

    controls <- causalTree.control(...)
    if (!missing(control)) controls[names(control)] <- control

    xval <- controls$xval
    if (is.null(xval) || (length(xval) == 1L && xval == 0L) || method=="user") {
      xgroups <- 0L
    	xval <- 0L
    } else if (length(xval) == 1L) {
        ## make random groups
        ################ here we debug only:
        control_idx <- which(treatment == 0)
        treat_idx <- which(treatment == 1)
        xgroups <- rep(0, nobs)
        xgroups[control_idx] <- sample(rep(1L:xval, length = length(control_idx)), length(control_idx), replace = F)
        xgroups[treat_idx] <- sample(rep(1L:xval, length = length(treat_idx)), length(treat_idx), replace = F)  
        # for debug:
        #xgroups <- c(9, 4, 5, 2, 6, 1, 8, 7, 3, 10, 1, 2, 7, 3, 5, 10, 4, 8, 9, 6)
        #xgroups <- c(8, 3, 5, 2, 6, 6, 7, 4, 9, 1, 9, 2, 7, 4, 10, 1, 8, 3, 5, 
         #            10, 8, 1, 7, 6, 9, 8, 2, 4, 6, 5, 4, 2, 3, 7, 5,9, 3, 10, 10, 1) 
        #print("xgroups = ")
        #print(xgroups)
    } else if (length(xval) == nobs) {
        ## pass xgroups by xval 
	      xgroups <- xval
	      xval <- length(unique(xgroups))
    } else {
        ## Check to see if observations were removed due to missing
	    if (!is.null(attr(m, "na.action"))) {
            ## if na.causalTree was used, then na.action will be a vector
	    temp <- as.integer(attr(m, "na.action"))
	    xval <- xval[-temp]
	    if (length(xval) == nobs) {
        xgroups <- xval
		    xval <- length(unique(xgroups))
            } else stop("Wrong length for 'xval'")
        } else stop("Wrong length for 'xval'")
    }

    ##
    ## Incorprate costs
    ##
    if (missing(cost)) cost <- rep(1, nvar)
    else {
	  if (length(cost) != nvar)
            stop("Cost vector is the wrong length")
	  if (any(cost <= 0)) stop("Cost vector must be positive")
    }

    ##
    ## Have C code consider ordered categories as continuous
    ##  A right-hand side variable that is a matrix forms a special case
    ## for the code.
    ##
    tfun <- function(x)
	  if (is.matrix(x)) rep(is.ordered(x), ncol(x)) else is.ordered(x)
    labs <- sub("^`(.*)`$", "\\1", attr(Terms, "term.labels")) # beware backticks
    isord <- unlist(lapply(m[labs], tfun))

    storage.mode(X) <- "double"
    storage.mode(wt) <- "double"
    storage.mode(treatment) <- "double"
    temp <- as.double(unlist(init$parms))
    #temp <- as.integer(init$parms)
    minsize <- as.integer(minsize)
    if (!length(temp)) temp <- 0  # if parms is NULL pass a dummy
    ctfit <- .Call(C_causalTree,
                   ncat = as.integer(cats * !isord),
                   method = as.integer(method.int),
                   as.double(unlist(controls)),
                   temp, # parms
                   minsize, #minsize = min_node_size
                   as.double(p),
                   as.integer(xval),
                   as.integer(xgroups),
                   as.double(t(init$y)),
                   X,
                   wt,
                   treatment,
                   as.integer(init$numy),
                   as.double(cost),
                   as.double(xvar),
                   as.double(alpha))

    nsplit <- nrow(ctfit$isplit) # total number of splits, primary and surrogate
    ## total number of categorical splits
    ncat <- if (!is.null(ctfit$csplit)) nrow(ctfit$csplit) else 0L
#    nodes <- nrow(ctfit$inode)
    if (nsplit == 0L) xval <- 0L # No xvals were done if no splits were found

    numcp <- ncol(ctfit$cptable)
    temp <- if (nrow(ctfit$cptable) == 3L) c("CP", "nsplit", "rel error")
            else c("CP", "nsplit", "rel error", "xerror", "xstd")
    dimnames(ctfit$cptable) <- list(temp, 1L:numcp)

    tname <- c("<leaf>", colnames(X))
    splits <- matrix(c(ctfit$isplit[, 2:3], ctfit$dsplit), ncol = 5L,
                     dimnames = list(tname[ctfit$isplit[, 1L] + 1L],
                     c("count", "ncat", "improve", "index", "adj")))
    index <- ctfit$inode[, 2L]  # points to the first split for each node

    ## Now, make ordered factors look like factors again (a printout choice)
    nadd <- sum(isord[ctfit$isplit[, 1L]])
    if (nadd > 0L) { # number of splits at an ordered factor.
      newc <- matrix(0L, nadd, max(cats))
	    cvar <- ctfit$isplit[, 1L]
	    indx <- isord[cvar]             # vector of TRUE/FALSE
	    cdir <- splits[indx, 2L]        # which direction splits went
	    ccut <- floor(splits[indx, 4L]) # cut point
	    splits[indx, 2L] <- cats[cvar[indx]] # Now, # of categories instead
	    splits[indx, 4L] <- ncat + 1L:nadd # rows to contain the splits

        ## Next 4 lines can be done without a loop, but become indecipherable
	    for (i in 1L:nadd) {
	        newc[i, 1L:(cats[(cvar[indx])[i]])] <- -as.integer(cdir[i])
	        newc[i, 1L:ccut[i]] <- as.integer(cdir[i])
          }
    	catmat <- if (ncat == 0L) newc
        else {
            ## newc may have more cols than existing categorical splits
            ## the documentation says that levels which do no exist are '2'
            ## and we add 2 later, so record as 0 here.
            cs <- ctfit$csplit
            ncs <- ncol(cs); ncc <- ncol(newc)
            if (ncs < ncc) cs <- cbind(cs, matrix(0L, nrow(cs), ncc - ncs))
            rbind(cs, newc)
        }
  	  ncat <- ncat + nadd
      } else catmat <- ctfit$csplit

    ## NB: package adabag depends on 'var' being a factor.
    if (nsplit == 0L) {   
      # tree with no splits
	    frame <- data.frame(row.names = 1L,
			    var = "<leaf>",
			    n = ctfit$inode[, 5L],
			    wt = ctfit$dnode[, 3L],
			    dev = ctfit$dnode[, 1L],
			    yval = ctfit$dnode[, 4L],
			    complexity = ctfit$dnode[, 2L],
			    ncompete = 0L,
			    nsurrogate = 0L)
    } else {
    	temp <- ifelse(index == 0L, 1L, index)
	    svar <- ifelse(index == 0L, 0L, ctfit$isplit[temp, 1L]) # var number
	    frame <- data.frame(row.names = ctfit$inode[, 1L],
                            ## maybe better to specify tname as the level?
			    var = tname[svar + 1L],
			    n = ctfit$inode[, 5L],
			    wt = ctfit$dnode[, 3L],
			    dev = ctfit$dnode[, 1L],
          yval = ctfit$dnode[, 4L],
			    complexity = ctfit$dnode[, 2L],
			    ncompete = pmax(0L, ctfit$inode[, 3L] - 1L),
			    nsurrogate = ctfit$inode[, 4L])
    }
    if (method.int == 3L) {
        ## Create the class probability vector from the class counts, and
        ##   add it to the results
        ## Also scale the P(T) result
        ## The "pmax" 3 lines down is for the case of a factor y which has
        ##   no one at all in one of its classes.  Both the prior and the
        ##   count will be zero, which led to a 0/0.
        numclass <- init$numresp - 2L
        nodeprob <- ctfit$dnode[, numclass + 5L] / sum(wt) # see ginidev.c
        temp <- pmax(1L, init$counts)   # overall class freq in data
        temp <- ctfit$dnode[, 4L + (1L:numclass)] %*% diag(init$parms$prior/temp)
        yprob <- temp /rowSums(temp)    # necessary with altered priors
        yval2 <- matrix(ctfit$dnode[, 4L + (0L:numclass)], ncol = numclass + 1L)
	      frame$yval2 <- cbind(yval2, yprob, nodeprob)
    } else if (init$numresp > 1L)
        frame$yval2 <- ctfit$dnode[, -(1L:3L), drop = FALSE]
        #print("frame$yval2:")
        #print(frame$yval2)
        

    if (is.null(init$summary))
        stop("Initialization routine is missing the 'summary' function")
    functions <- if (is.null(init$print)) list(summary = init$summary)
                 else list(summary = init$summary, print = init$print)
    if (!is.null(init$text)) functions <- c(functions, list(text = init$text))
    if (method == "user") functions <- c(functions, mlist)

    where <- ctfit$which
    names(where) <- row.names(m)

    ans <- list(frame = frame,
             where = where,
             call = Call, terms = Terms,
             cptable = t(ctfit$cptable),
             method = method,
             parms = init$parms,
             control = controls,
             functions = functions,
             numresp = init$numresp)
    if (nsplit) ans$splits = splits
    if (ncat > 0L) ans$csplit <- catmat + 2L
    if (nsplit) ans$variable.importance <- importance(ans)
    if (model) {
	ans$model <- m
	if (missing(y)) y <- FALSE
    }
    if (y) ans$y <- Y
    if (x) {
	ans$x <- X
	ans$wt <- wt
    }
    ans$ordered <- isord
    if (!is.null(attr(m, "na.action"))) ans$na.action <- attr(m, "na.action")
    if (!is.null(xlevels)) attr(ans, "xlevels") <- xlevels
    if (method == "class") attr(ans, "ylevels") <- init$ylevels
    class(ans) <- "causalTree"
    ans
}
