#setwd("C:/Users/athey/Dropbox/SusanGuido/Datasets/Simulations")

##install.packages("devtools")
##install.packages("Rcpp")
##install.packages("sandwich")
#library(sandwich)

#library(rattle)
#library(ROCR)
#library(knitr)


##install_github("hadley/devtools")

#library(devtools)

##install_github("shuihu/causal_effect_estimation",  
##               auth_token = "44ec9832e53cfbb807d1faa9fcde10f29b297491"
##)

#source("honestTree.R")

##install.packages("C:/Users/athey/Documents/GitHub/causal_effect_estimation")

library(causalTree)
library(rpart)
library(rpart.plot)

#SA added


#### small simulation for CT methods:




## generate simulation data set:
generate <- function(NUM, p, design) {
  name <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "y", "w", "eff")
  ## for design 1:
  Xa <- matrix(rnorm(NUM * p), nrow = NUM, ncol = p)
  Xb <- matrix(rnorm(NUM * p), nrow = NUM, ncol = p)
  Xc <- matrix(rnorm(NUM * p), nrow = NUM, ncol = p)  
  
  wa <- rep(1:0, each = NUM/2)
  wb <- rep(1:0, each = NUM/2)
  wc <- rep(1:0, each = NUM/2)
  wtreat <- rep(1, each=NUM)
  wcon <- rep(0, each=NUM)
  wa_ <- rep(0:1, each = NUM/2)
  wb_ <- rep(0:1, each = NUM/2)
  wc_ <- rep(0:1, each = NUM/2)
  if (design == 1){
    ## y1 ~ N(1 - x1 + x2, 1) 
    ## y0 ~ N(x1 + x2, 1)
    ## f = y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10  
    ya <- wa + (1 - 2 * wa) *Xa[,1] + Xa[,2] + rnorm(NUM, 0, 1)
    yb <- wb + (1 - 2 * wb) *Xb[,1] + Xb[,2] + rnorm(NUM, 0, 1)
    yc <- wc + (1 - 2 * wc) *Xc[,1] + Xc[,2] + rnorm(NUM, 0, 1)
    ya_ <- wa_ + (1 - 2 * wa_) *Xa[,1] + Xa[,2] + rnorm(NUM, 0, 1)
    yb_ <- wb_ + (1 - 2 * wb_) *Xb[,1] + Xb[,2] + rnorm(NUM, 0, 1)
    yc_ <- wc_ + (1 - 2 * wc_) *Xc[,1] + Xc[,2] + rnorm(NUM, 0, 1)
    f <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10  
  } else if (design == 2) {
    ## y1 ~ N(1 - x1 + x2 + x3 + ... + x10, 1)
    ## y0 ~ N(x1 + x2 + x3 + ..., + x10, 1)
    ## f = y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10  
    ya <- wa + (1 - 2 * wa) *Xa[,1] + Xa[,2] + rowSums(Xa[,3:p]) + rnorm(NUM, 0, 1)
    yb <- wb + (1 - 2 * wb) *Xb[,1] + Xb[,2] + rowSums(Xb[,3:p]) + rnorm(NUM, 0, 1)
    yc <- wc + (1 - 2 * wc) *Xc[,1] + Xc[,2] + rowSums(Xc[,3:p]) + rnorm(NUM, 0, 1)
    ya_ <- wa_ + (1 - 2 * wa_) *Xa[,1] + Xa[,2] + rowSums(Xa[,3:p]) + rnorm(NUM, 0, 1)
    yb_ <- wb_ + (1 - 2 * wb_) *Xb[,1] + Xb[,2] + rowSums(Xb[,3:p]) + rnorm(NUM, 0, 1)
    yc_ <- wc_ + (1 - 2 * wc_) *Xc[,1] + Xc[,2] + rowSums(Xc[,3:p]) + rnorm(NUM, 0, 1)
    f <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10  
    
  } else if (design == 3) {
    ## y1 ~ N(1 - x1 + x2, 0.1)
    ## y0 ~ N(x1 + x2, 0.1)
    ## f = y ~ x1 + x2 
    ya <- wa + (1 - 2 * wa) *Xa[,1] + Xa[,2] + rnorm(NUM, 0, 0.1)
    ya_ <- wa_ + (1 - 2 * wa_) *Xa[,1] + Xa[,2] + rnorm(NUM, 0, 0.1)
    yb <- wb + (1 - 2 * wb) *Xb[,1] + Xb[,2] + rnorm(NUM, 0, 0.1)
    yb_ <- wb_ + (1 - 2 * wb_) *Xb[,1] + Xb[,2] + rnorm(NUM, 0, 0.1)
    yc <- wc + (1 - 2 * wc) *Xc[,1] + Xc[,2] + rnorm(NUM, 0, 0.1)
    yc_ <- wc_ + (1 - 2 * wc_) *Xc[,1] + Xc[,2] + rnorm(NUM, 0, 0.1)
    f <- y ~ x1 + x2  
  }else if(design == 4){
    ## y1 ~ N(1 - x1 + x2, 0.1)
    ## y0 ~ N(x1 + x2, 0.1)
    ## f = y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10  
    ya <- wa + (1 - 2 * wa) *Xa[,1] + Xa[,2] + rnorm(NUM, 0, 0.1)
    yb <- wb + (1 - 2 * wb) *Xb[,1] + Xb[,2] + rnorm(NUM, 0, 0.1)
    yc <- wc + (1 - 2 * wc) *Xc[,1] + Xc[,2] + rnorm(NUM, 0, 0.1)
    ya_ <- wa_ + (1 - 2 * wa_) *Xa[,1] + Xa[,2] + rnorm(NUM, 0, 0.1)
    yb_ <- wb_ + (1 - 2 * wb_) *Xb[,1] + Xb[,2] + rnorm(NUM, 0, 0.1)
    yc_ <- wc_ + (1 - 2 * wc_) *Xc[,1] + Xc[,2] + rnorm(NUM, 0, 0.1)
    f <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10
  }else if(design == 5){
    ## y1 ~ N(x1 + x2, 0.1)
    ## y0 ~ N(x1 + x2, 0.1)
    ## f = y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10  
    ya <- Xa[,1] + Xa[,2] + rnorm(NUM, 0, 0.1)
    yb <- Xb[,1] + Xb[,2] + rnorm(NUM, 0, 0.1)
    yc <- Xc[,1] + Xc[,2] + rnorm(NUM, 0, 0.1)
    ya_ <- Xa[,1] + Xa[,2] + rnorm(NUM, 0, 0.1)
    yb_ <- Xb[,1] + Xb[,2] + rnorm(NUM, 0, 0.1)
    yc_ <- Xc[,1] + Xc[,2] + rnorm(NUM, 0, 0.1)
    f <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10
  }
  # used for Q-infeasible for sample C:
  eff_a <- ya_ - ya
  eff_b <- yb_ - yb
  eff_c <- yc_ - yc
  eff_a[1:(NUM/2)] <- -eff_a[1:(NUM/2)]
  eff_b[1:(NUM/2)] <- -eff_b[1:(NUM/2)]
  eff_c[1:(NUM/2)] <- -eff_c[1:(NUM/2)]
  A <- data.frame(Xa, ya, wa, eff_a)
  B <- data.frame(Xb, yb, wb, eff_b)
  C <- data.frame(Xc, yc, wc, eff_c)
  Atreat <- data.frame(Xa, ya, wtreat, eff_a)
  Acon <- data.frame(Xa, ya, wcon, eff_a)  
  Btreat <- data.frame(Xb, yb, wtreat, eff_b)
  Bcon <- data.frame(Xb, yb, wcon, eff_b) 
  Ctreat <- data.frame(Xc, yc, wtreat, eff_c)
  Ccon <- data.frame(Xc, yc, wcon, eff_c)
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
                  Bcon=Bcon, Ctreat=Ctreat, Ccon=Ccon, f = f, w = wa)
  return (newlist)
}