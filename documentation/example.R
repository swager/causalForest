#must install devtools, rpart, and rpart.plot in order to install packages from github
#once packages are installed once, comment out the installation lines in your script
#you need only install once, but will need to load the libraries using the library()
#function everytime you start RStudio

install.packages("devtools")
install.packages("rpart", dependencies=TRUE, repos='http://cran.us.r-project.org')
install.packages("rpart.plot", dependencies=TRUE, repos='http://cran.us.r-project.org')
library(devtools)  #Load the installed "devtools" package
library(rpart)   #Load the installed "rparts" package
library(rpart.plot)

#again, these can be commented out once they have been run once
install_github("swager/randomForest")
install_github("swager/randomForestCI")
install_github("swager/causalForest", auth_token = "<insert your personal authentication token here>")
install.packages("rattle", dependencies=TRUE,repos='http://cran.us.r-project.org')
install.packages("knitr", dependencies=TRUE,repos='http://cran.us.r-project.org')
install.packages("ggplot2", dependencies=TRUE, repos='http://cran.us.r-project.org')

#load all packages needed 
library(randomForest)  
library(randomForestCI)  
library(rattle)
library(causalForest)
library(ggplot2)

#generating a sample dataset to on which to demonstrate functions.. ignore up to line 
#number of observations
NUM=500
#number of variables
p=10
#variance
var1=1
var2=1
#number of treated variables
treatsize=200
coveff=5
testsize=100
## this generates simulation data set for demonstrating the use of the functions
generate <- function(NUM, p, var1, var2, treatsize, coveff, testsize=-99) {
  if (testsize==-99) {
    testsize=NUM
  }
  ## for design 1:
  totalobs = NUM*2+testsize
  
  X <- matrix(rnorm(NUM * p), nrow = totalobs, ncol = p)
  
  w <- rep(c(0,1), totalobs/2)  
  wtreat <- as.integer(rep(1, each=totalobs))
  wcon <- as.integer(rep(0, each=totalobs))
  w_ <- 1-w
  
  
  y <- 2 * w - 1 + treatsize*(1+2 * w) *X[,1]+treatsize*(1+2 * w) *X[,2] + coveff*X[,3] + rnorm(NUM, 0, var2)
  y_ <- 2 * w_ - 1 + treatsize*(1+2 * w_) *X[,1]+treatsize*(1+2 * w_) *X[,2] + coveff*X[,3] + rnorm(NUM, 0, var2)
  numx <- p
  
  f <- ""
  nextx <- ""
  if (numx>1) {
    numxm1 <- numx-1
    for (ii in 1:numxm1) {
      nextx <- paste("x",ii, sep="")
      if (ii==1) {name <- nextx}
      if (ii>1) {name <- c(name, nextx)}
      f <- paste(f, nextx, "+", sep="")
    }
    f <- paste(f, "x", ii+1, sep="")
  } else if (numx==1) {
    f <- "x1"
  }
  
  for (ii in 1:p) {
    nextx <- paste("x",ii, sep="")
    if (ii==1) {name <- nextx}
    if (ii>1) {name <- c(name, nextx)}
  }
  
  name <- c( name,  "y", "w", "eff")
  
  # used for Q-infeasible for sample C:
  eff <- y_ - y
  A <- data.frame(X[1:NUM,], y[1:NUM], w[1:NUM], eff[1:NUM])
  B <- data.frame(X[(NUM+1):(2*NUM),], y[(NUM+1):(2*NUM)], w[(NUM+1):(2*NUM)], eff[(NUM+1):(2*NUM)])
  C <- data.frame(X[(2*NUM+1):totalobs,], y[(2*NUM+1):totalobs], w[(2*NUM+1):totalobs], eff[(2*NUM+1):(totalobs)])
  
  Atreat <- data.frame(X[1:NUM,], y[1:NUM], wtreat[1:NUM], eff[1:NUM])
  Acon <- data.frame(X[1:NUM,], y[1:NUM], wcon[1:NUM], eff[1:NUM])
  Btreat <- data.frame(X[(NUM+1):(2*NUM),], y[(NUM+1):(2*NUM)], wtreat[(NUM+1):(2*NUM)], eff[(NUM+1):(2*NUM)])
  Ctreat <- data.frame(X[(2*NUM+1):totalobs,], y[(2*NUM+1):totalobs], wtreat[(2*NUM+1):totalobs], eff[(2*NUM+1):(totalobs)])
  Bcon <- data.frame(X[(NUM+1):(2*NUM),], y[(NUM+1):(2*NUM)], wcon[(NUM+1):(2*NUM)], eff[(NUM+1):(2*NUM)])
  Ccon <- data.frame(X[(2*NUM+1):totalobs,], y[(2*NUM+1):totalobs], wcon[(2*NUM+1):totalobs], eff[(2*NUM+1):(totalobs)])
  
  names(A) <- name
  names(B) <- name
  names(C) <- name
  names(Atreat) <- name
  names(Btreat) <- name
  names(Ctreat) <- name
  names(Acon) <- name
  names(Bcon) <- name
  names(Ccon) <- name
  newlist <- list(A = A, B = B, C = C, Atreat=Atreat, Acon=Acon, Btreat=Btreat, 
                  Bcon=Bcon, Ctreat=Ctreat, Ccon=Ccon, f = f, w = w, numx = numx)
  return (newlist)
}


Data <- generate(NUM, p, var1, var2, treatsize, coveff)

#extract samples A, B and C
#tree building dataset
A<-Data$A
#leaf estimating dataset
B<-Data$B
#test dataset
C<-Data$C

#for this example, there are 10 covariates, named x1, x2, x3, x4 and so on, up to x10
#we build a formula using a string
sumx= paste(head(colnames(A),-3), collapse=" + ")
#this returns"x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10"
linear<-paste("y",paste("w",sumx, sep=" + "), sep=" ~ ")
#concatenates the full formula, returns "y ~ w + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10"
linear
#converts to full formula
linear<-as.formula(linear)

linear.caustree<-causalTree(linear, data=A, treatment=A$w, method = "anova", control=rpart.control(cp = 0, xval = 10), parms = 1,minsize = 10, 
                                    cv.option = "matching", p = 0.5)
linear.caustree.pruned<-prune(linear.caustree, 
                                     cp = linear.caustree$cptable[which.min(linear.caustree$cptable[,"rel error"]),"CP"])
linear.caustree$cptable

#the honest model functions require dataset to have a "treatment" column
A$treatment<-A$w
B$treatment<-B$w
C$treatment<-C$w

linear.tree<-rpart(linear, data=A, method = "anova", control=rpart.control(cp = 0, xval = 10), parms = 1)


#we now use the B sample, to do leaf estimation and build an honest tree
#can then view the subgroups, and/or predict on the test dataset C with 
#these trees
honestCTree<-honestCTree(object=linear.caustree, data=B, treatment=B$treatment)
honestCTree
honestST<-honestST(object=linear.tree, data=B)
honestST
honestTOT<-honestTOT(object=linear.caustree, data=B, treatment=B$treatment)
honestTOT


#computed an estimate of the treatment effect using the Q^TOT method
compute.ystar = function(Y, W){
  ystar<-NULL
  prob_treated<-mean(W)
  for (idx in (1:length(Y))) {
    #if treatment
    if (W[idx]==1) {
      ystar[idx]<-Y[idx]/prob_treated
      #if control 
    } else if (W[idx]==0){
      ystar[idx]<-(-Y[idx])/(1-prob_treated)
    }
  }
  return(ystar)
}

#computing standard errors
A.leaves <- predict(linear.caustree, A, type = 'vector')
B.leaves <- predict(linear.caustree, B, type = 'vector')
C.leaves <- predict(linear.caustree, C, type = 'vector')

A.leavesf <- factor(round(A.leaves,4))
B.leavesf <- factor(round(B.leaves,4))
C.leavesf <- factor(round(C.leaves,4))

#use Q infeasible method to get an estimate of treatment effects based on Y, needed for the error estimate
sampleA <- data.frame(leavesf = A.leavesf, y = A$y, w =A$w, eff=compute.ystar(A$y, A$w))
sampleB <- data.frame(leavesf = B.leavesf, y = B$y, w =B$w, eff=compute.ystar(B$y, B$w))
sampleC <- data.frame(leavesf = C.leavesf, y = C$y, w =C$w, eff=compute.ystar(C$y, C$w))


modelA <- lm(y~-1+leavesf+leavesf*w-w, data=sampleA)
modelB <- lm(y~-1+leavesf+leavesf*w-w, data=sampleB)
modelC <- lm(y~-1+leavesf+leavesf*w-w, data=sampleC)

#modelA estimates should be the same as the predictions using A.leaves
#these regressions give you the treatment effect estimates, and the standard error estimates 
summary(modelA)
summary(modelB)
summary(modelC)





#renaming the covariate names to be integers, as this is required for the causal forest function 
A_temp<-A
names(A_temp)[1:10]<-seq(1,10)
B_temp<-B
names(B_temp)[1:10]<-seq(1,10)
C_temp<-C
names(C_temp)[1:10]<-seq(1,10)

#Causal Forest handles 'honesty', on its own, so instead of the two step process, as is done above, we pass 
#causal forest 90% of the dataset, by combining the A and B datasets
cf_training<-rbind(A_temp, B_temp)
#set number of trees, and sample size
num.trees <- 500
sample.size <- 225
cv.option <- "matching"
#build causal forest, first parameter is the dataframe of covariates
forest<-causalForest(cf_training[,1:10], cf_training$y, cf_training$w, num.trees, sample.size, cv.option)
#make predictions on the test sample
predictions = predict(forest, C_temp)
#get confidence intervals on the predictions 
forest.ci = randomForestInfJack(forest, C_temp, calibrate = TRUE)









