############################################
source("honestModeltmp.R")
source("generate_data.R")
## set parameters
minsize <- 20 
minbucket <- 5
NUM = 500
p = 10
Iter <- 4
modnum_lo <- 1# TOT-TOT, CT-TOT, CT-matching not implement ST and TT 
modnum_hi <- 5
numdesign <- 3 # generate data design
numoutputs <- 16

ans <- array(0, c(modnum_hi, numdesign, numoutputs)) # model num, design num, outputs
# ans <- data.frame(ans)

leafcount <- array(0,c(modnum_hi,numdesign,Iter))
Q_dish <- array(0,c(modnum_hi,numdesign,Iter))
Q_h <- array(0,c(modnum_hi,numdesign,Iter))
In_dish <- array(0,c(modnum_hi,numdesign,Iter))
In_h <- array(0,c(modnum_hi,numdesign,Iter)) 
In_hdish <- array(0,c(modnum_hi,numdesign,Iter))
In_trtest_h <- array(0,c(modnum_hi,numdesign,Iter))
In_trdish_dish <- array(0,c(modnum_hi,numdesign,Iter))
In_trh_h <- array(0,c(modnum_hi,numdesign,Iter))
In_trtest_dish <- array(0,c(modnum_hi,numdesign,Iter))
In_trtest_h <- array(0,c(modnum_hi,numdesign,Iter))
In_trtest_dish <- array(0,c(modnum_hi,numdesign,Iter))

win_infeas_h <- array(0,c(modnum_hi,numdesign,Iter))
win_infeas_dish <- array(0,c(modnum_hi,numdesign,Iter))


## initialie matrix of leaf counts, Q-infeas, ... 

## model:
# train dishonest causal tree
# do cross-validation
# get the honoest values
# count leaves in tree A
# calcaulate the Q-infeas
for (k in 1:numdesign){
  design = k
  print(paste0("Current working design: ", k))
  for (i in 1:Iter){  
    Data <- generate(NUM, p, design)
    print(paste0("Current iteration:", i))
    for (modnum in modnum_lo:modnum_hi){
      if(modnum==3){
        splitmeth = "TOT"  
        cvmeth = "TOT"
      } else if(modnum==4){
        splitmeth = "CT"
        cvmeth = "TOT"
      } else if(modnum==5){
        splitmeth = "CT"
        cvmeth = "matching"
      }
      
      # dishonest tree model:
      #set.seed(Iter-i)
      #if (i%%100 == 0) print (i)
      
      if (modnum==4 | modnum==5 | modnum == 3) {
        treeA <- causalTree(as.formula(Data$f), data = Data$A, treatment = Data$A$w, split.option = splitmeth, 
                            cv.option = cvmeth, cp = 0, minsize = minsize, minbucket = minbucket)
      }
      
      
      if (modnum==1){  #fix TT model later
        # ST: y ~ x + w
        treeA <- rpart(as.formula(paste("y~",paste(Data$f,"+w")[3])), data=Data$A, 
                       method="anova",
                       control=rpart.control(minsplit = minsize*2, 
                                             minbucket = minbucket, cp = 0, 
                                             xval = 10))    
      }
      if (modnum == 2) {
        ################### <<<<<<<<<<<<<<<<<<<< TT TT TT TT TT TT  <<<<<<<<<<<<<
      }
      
      if (modnum != 2) {
        opcpid <- which.min(treeA$cp[,4])
        opcp <- treeA$cp[opcpid,1]
        treeA_prune <- prune(treeA, cp = opcp)   
        leafcount[modnum,k,i] = length(unique(treeA_prune$where))
        
        if (modnum == 1){
          A.leavesraw <- predict(treeA_prune, Data$A, type = 'vector')
          B.leavesraw <- predict(treeA_prune, Data$B, type = 'vector')
          C.leavesraw <- predict(treeA_prune, Data$C, type = 'vector')
            
          A.leavestr <- predict(treeA_prune, Data$Atreat, type = 'vector')
          B.leavestr <- predict(treeA_prune, Data$Btreat, type = 'vector')
          C.leavestr <- predict(treeA_prune, Data$Ctreat, type = 'vector')
          
          
          A.leavescon <- predict(treeA_prune, Data$Acon, type = 'vector')
          B.leavescon <- predict(treeA_prune, Data$Bcon, type = 'vector')
          C.leavescon <- predict(treeA_prune, Data$Ccon, type = 'vector')
          
          ## the way we compute the causal effect:
          A.leaves <- A.leavestr - A.leavescon
          B.leaves <- B.leavestr - B.leavescon
          C.leaves <- C.leavestr - C.leavescon
          
          
          A.leavesrawf <- factor(round(A.leavesraw,4))
          B.leavesrawf <- factor(round(B.leavesraw,4))
          C.leavesrawf <- factor(round(C.leavesraw,4))
          
          A.leavestrf <- factor(round(A.leavestr,4))
          B.leavestrf <- factor(round(B.leavestr,4))
          C.leavestrf <- factor(round(C.leavestr,4))
          
          
          A.leavesconf <- factor(round(A.leavescon,4))
          B.leavesconf <- factor(round(B.leavescon,4))
          C.leavesconf <- factor(round(C.leavescon,4))
          
          
          A.leavesf <- factor(round(A.leaves,4))
          B.leavesf <- factor(round(B.leaves,4))
          C.leavesf <- factor(round(C.leaves,4))
          
          dish_est <- C.leaves
          #dish_est <- predict(treeA_prune, Data$C)         
          
        } else if (modnum >= 3){ 
          
          A.leaves <- predict(treeA_prune, Data$A, type = 'vector')
          B.leaves <- predict(treeA_prune, Data$B, type = 'vector')
          C.leaves <- predict(treeA_prune, Data$C, type = 'vector')
          
          
          A.leavesf <- factor(round(A.leaves,4))
          B.leavesf <- factor(round(B.leaves,4))
          C.leavesf <- factor(round(C.leaves,4))
          
          dish_est <- C.leaves
        }
         
        sampleA <- data.frame(leavesf = A.leavesf, y = Data$A$y, w = Data$A$w, eff=Data$A$eff)
        sampleB <- data.frame(leavesf = B.leavesf, y = Data$B$y, w=Data$B$w, eff=Data$B$eff)
        sampleC <- data.frame(leavesf = C.leavesf, y = Data$C$y, w = Data$C$w, eff=Data$C$eff)
        
        
        sampleBraw <- data.frame(leavesf = B.leavesrawf, y = Data$B$y, w=Data$B$w, eff=Data$B$eff)
        
        sampleBtr <- data.frame(leavesf = B.leavestrf, y = Data$B$y, w=Data$B$w, eff=Data$B$eff)
        
        sampleBcon <- data.frame(leavesf = B.leavesconf, y = Data$B$y, w=Data$B$w, eff=Data$B$eff)
        
        
        sampleCtr <- data.frame(leavesf = C.leavestrf, y = Data$B$y, w=Data$B$w, eff=Data$B$eff)
        
        sampleCcon <- data.frame(leavesf = C.leavesconf, y = Data$B$y, w=Data$B$w, eff=Data$B$eff)
        
        
      }
      if (modnum==2){ 
        #fix this later for TT
        opcpid <- which.min(treeA$cp[,4])
        opcp <- treeA$cp[opcpid,1]
        treeA_prune <- prune(treeA, cp = opcp)
        dish_est <- predict(treeA_prune, Data$C)
        
        
        leafcount[modnum,k,i] = length(unique(treeA_prune$where))
        
      }
      
      dish_err <- sum((Data$C$eff - dish_est)^2)
      Q_dish[modnum,k,i] <- dish_err
      
      
      ######### honest tree model:
      
      if (modnum >3) {
        treeB <- honestCTree(treeA_prune, data = Data$B, treatment = Data$B$w)
        h_est <- predict(treeB, data = Data$C)
      }
      
      if (modnum == 3) {
        treeB <- honestTOT(treeA_prune, data = Data$B, treatment = Data$B$w)
        h_est <- predict(treeB, data = Data$C)
      }
      
      
      if (modnum == 1){ 
        treeB <- honestST(treeA_prune, data = Data$B)
        treatpred <- rep(0, 500)
        conpred <- rep(0, 500)
        for (j in 1:500){
          treatpred[j] <- predict(treeB, Data$Ctreat[i,])
          conpred[j] <- predict(treeB, Data$Ccon[i, ])
        }
        h_est <- treatpred - conpred     
        
        "
        if (length(table(sampleBraw$leavesf))>1) {
          honestSTpred <- lm(y~leavesf-1, data = sampleBraw)
        } else {
          honestSTpred <- lm(y~1, data = sampleBraw)   
        }
        treatpred <- predict(honestSTpred, data = sampleCtr)
        conpred <- predict(honestSTpred, data = samlepleCcon)
        h_est <- treatpred - conpred
        "
      }
      
      ################### <<<<<<<<<<<<<<<<<<<< TT TT TT TT TT TT 
      if (modnum==2){ #fill in TT later--just placeholder
        
        "
        if (length(table(sampleB$leavesf))>1) {
          honestSTpred <- lm(y~leavesf, data = sampleB)
        } else {
          honestSTpred <- lm(y~1, data = sampleB)  
        }
        wtemp <- Data$C$w
        Data$C$w <- rep(1,NUM)
        treatpred <- predict(treeA_prune, data = Data$C)
        Data$C$w <- rep(0,NUM)
        conpred <- predict(treeA_prune, data = Data$C)
        Data$C$w <- wtemp
        h_est <- treatpred - conpred "

      }
      
      
      h_err <- sum((Data$C$eff - h_est)^2)
      Q_h[modnum,k,i] <- h_err
      
      
      
      ############################
      
      # confidence interval, ONLY FOR CT, CT-TOT models
      
      if (modnum >= 4) {
        
        
        if (length(levels(A.leavesf)) == 1){
          
          modelA <- lm(y~w, data=sampleA)
          modelB <- lm(y~w, data=sampleB)
          modelC <- lm(y~w, data=sampleC)
          
          
          A.coeftr <- coef(lm(eff~w, data=sampleA))[2]
          B.coeftr <- coef(lm(eff~w, data=sampleB))[2]
          C.coeftr <- coef(lm(eff~w, data=sampleC))[2]
          
          coefnumh <- 2
          coefnuml <- 2
        } else{
          
          
          modelA <- lm(y~-1+leavesf+leavesf*w-w, data=sampleA)
          modelB <- lm(y~-1+leavesf+leavesf*w-w, data=sampleB)
          modelC <- lm(y~-1+leavesf+leavesf*w-w, data=sampleC)
          
          coefnumh <- length(coef(modelC))
          coefnuml <- length(coef(modelC))/2 + 1
          
          A.coeftr <- coef(lm(eff~leavesf-1, data=sampleA))
          B.coeftr <- coef(lm(eff~leavesf-1, data=sampleB))
          C.coeftr <- coef(lm(eff~leavesf-1, data=sampleC))
        }      
        
        ## only use the parameters of leavesf
        C.coef <- coef(modelC)[coefnuml:coefnumh]
        # SA adds B.coef
        B.coef <- coef(modelB)[coefnuml:coefnumh]
        A.coef <- coef(modelA)[coefnuml:coefnumh]
        
        
        A.stderr <- sqrt(diag(vcov(modelA)))[coefnuml:coefnumh]
        B.stderr <- sqrt(diag(vcov(modelB)))[coefnuml:coefnumh]
        C.stderr <- sqrt(diag(vcov(modelC)))[coefnuml:coefnumh]
        
        
        ## SA changed 500 to NUM      
        C.fraction <- table(sampleC$leavesf)/NUM
        ## SA changed for testing
        #  C.fraction <- rep(1/length(C.fraction),length(C.fraction))
        
        conffact <- 1.96 * 2^(1/2)
        
        # dishonest_confidence_interval:
        A.coef.up <- A.coef + conffact * A.stderr  
        A.coef.lo <- A.coef - conffact * A.stderr
        # honest_confidence_interval:
        B.coef.up <- B.coef + conffact * B.stderr
        B.coef.lo <- B.coef - conffact * B.stderr
        
        # compute the fraction of in's.
        In_dish[modnum,k,i] <- sum(C.fraction * (C.coef <= A.coef.up & C.coef >= A.coef.lo))
        In_h[modnum,k,i] <- sum(C.fraction * (C.coef <= B.coef.up & C.coef >= B.coef.lo))  
        In_hdish[modnum,k,i] <- sum(C.fraction * (B.coef <= A.coef.up & B.coef >= A.coef.lo))
        
        
        conffact <- 1.96 
        
        ## SA adds B.coef
        
        # dishonest_confidence_interval:
        A.coef.up <- A.coef + conffact * A.stderr  
        A.coef.lo <- A.coef - conffact * A.stderr
        # honest_confidence_interval:
        B.coef.up <- B.coef + conffact * B.stderr
        B.coef.lo <- B.coef - conffact * B.stderr
        
        # compute the fraction of in's.
        In_trdish_dish[modnum,k,i] <- sum(C.fraction * (A.coeftr <= A.coef.up & A.coeftr >= A.coef.lo))
        In_trh_h[modnum,k,i] <- sum(C.fraction * (B.coeftr <= B.coef.up & B.coeftr >= B.coef.lo))  
        In_trtest_dish[modnum,k,i] <- sum(C.fraction * (C.coeftr <= A.coef.up & C.coeftr >= A.coef.lo))
        In_trtest_h[modnum,k,i] <- sum(C.fraction * (C.coeftr <= B.coef.up & C.coeftr >= B.coef.lo))  
        
      } # end confidence interval code
      
    } #Modnum loop
    
    max_h <- max(Q_h[,design,i])
    max_dish <- max(Q_dish[,design,i])
    
    for (modnum in modnum_lo:modnum_hi) {
      win_infeas_h[modnum,design,i]=as.numeric(Q_h[modnum,design,i]==max_h)
      win_infeas_dish[modnum,design,i]=as.numeric(Q_dish[modnum,design,i]==max_dish)
    }
    
    #ties
    for (modnum in modnum_lo:modnum_hi) {
      win_infeas_h[modnum,design,i]=win_infeas_h[modnum,design,i]/sum(win_infeas_h[,design,i])
      win_infeas_dish[modnum,design,i]=win_infeas_dish[modnum,design,i]/sum(win_infeas_dish[,design,i])
    }
    
  } # Iter loop
  
  for (modnum in modnum_lo:modnum_hi) {
    
    ans[modnum,design,1] <- design
    ans[modnum,design,2] <- mean(leafcount[modnum,design,], na.rm = T)
    ##SA added median  
    ans[modnum,design,3] <- median(leafcount[modnum,design,], na.rm = T)
    ans[modnum,design,4] <- mean(Q_dish[modnum,design,], na.rm = T)/ NUM
    ans[modnum,design,5] <- median(Q_dish[modnum,design,], na.rm = T)/ NUM
    ans[modnum,design,6] <- mean(Q_h[modnum,design,], na.rm = T)/NUM
    ans[modnum,design,7] <- median(Q_h[modnum,design,], na.rm = T)/NUM
    ans[modnum,design,8] <- mean(In_dish[modnum,design,], na.rm = T)
    ans[modnum,design,9] <- mean(In_h[modnum,design,], na.rm = T)
    ##SA added honest in dishonest, truth in honest/dishonest
    ans[modnum,design,10] <- mean(In_hdish[modnum,design,], na.rm=T)
    ans[modnum,design,11] <- mean(In_trdish_dish[modnum,design,], na.rm = T)
    ans[modnum,design,12] <- mean(In_trh_h[modnum,design,], na.rm = T)
    ans[modnum,design,13] <- mean(In_trtest_dish[modnum,design,], na.rm = T)
    ans[modnum,design,14] <- mean(In_trtest_h[modnum,design,], na.rm = T)
    ans[modnum,design,15] <- mean(win_infeas_dish[modnum,design,], na.rm = T)
    ans[modnum,design,16] <- mean(win_infeas_h[modnum,design,], na.rm = T)
  } #modnum small loop
  
  
} #design loop

for (modnum in modnum_lo:modnum_hi) {
  
  
  
  anstemp <- ans[modnum,,]
  anstemp <- data.frame(anstemp)
  names(anstemp) <- c("design","mean_leaf_number", "med_leaf_number", "Q_dishonest_mean", "Q_dishonest_median","Q_honest_mean", 
                      "Q_honest_median","Fraction_dishonest", "Fraction_honest", "Frac_hon_in_dish", "Frac_trdish_dish", "Frac_trh_h", "Frac_trtest_dishonest", "Frac_trtest_honest",
                      "Frac_win_dish","Frac_win_h")
  anstemp <- t(anstemp)
  
  
  if(modnum==3){
    ansTOT <- anstemp
  } else if(modnum==4){
    ansCTTOT <- anstemp
  } else if(modnum==5){
    ansCT <- anstemp
  } else if(modnum==1) {
    ansST <- anstemp
  } else if(modnum==2) {
    ansTT <- anstemp
  }
  
} # modnum loop

ansTOT
ansCTTOT
ansCT
ansST

