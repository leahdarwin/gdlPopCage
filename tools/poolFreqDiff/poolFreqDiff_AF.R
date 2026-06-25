#####
# Helper functions for GLM test.
#
# Author: R. Axel W. Wiberg, Leah Darwin (added labels for my experiment)
#
#####

library(dplyr)

cat("# Loaded: poolFreqDiffTest.R\n")
# Source the G-test
cat("# Looking in: ",currdir," for G_test.R\n")
source(paste(currdir,"/G_test.R",sep=""))
#
#
# FUNCTION: Woolf-test
# The script comes from the help page for the mantelhaen.test()
# ?mantelhaen.test()
woolf.test <- function(x) {
  x <- x + 1 / 2
  k <- dim(x)[3]
  or <- apply(x, 3, function(x) (x[1,1]*x[2,2])/(x[1,2]*x[2,1]))
  w <-  apply(x, 3, function(x) 1 / sum(1 / x))
  woolf <- sum(w * (log(or) - weighted.mean(log(or), w)) ^ 2)
  df <- k-1
  p <- 1 - pchisq(woolf, df)
  dat <- c(woolf,df,p)
  names(dat) <- c("Woolf", "df", "p-value")
  dat
}

# FUNCTION: Get GLM results from array ##
#get GLM data-set from k-way table
get_glm_dat <- function(array,zeroes=1){
  if(zeroes == 1){
    if(any(array == 0)){
      array <- array+1
    }
  }
  A_Cnt <- vector(length = dim(array)[3]*dim(array)[1])
  Tot_Cnt <- vector(length = dim(array)[3]*dim(array)[1])
  tr_l <- vector(length = dim(array)[3]*dim(array)[1])
  rep <- vector(length = dim(array)[3]*dim(array)[1])
  j <-1
  for(k in seq(1,dim(array)[3],1)){
    for(i in seq(1,dim(array)[1],1)){
      #      print(c(i,j,k))
      A_Cnt[j]<-array[i,1,k]
      Tot_Cnt[j]<-sum(array[i,,k])
      tr_l[j] <- as.character(i)
      rep[j] <- as.character(k)
      j <- j + 1
    }
  }
  d<-data.frame("A_Cnt"=A_Cnt,"Tot_Cnt"=Tot_Cnt,"tr_l"=tr_l,"rep"=rep)
  mod <- anova(glm(
    cbind(d$A_Cnt,d$Tot_Cnt-d$A_Cnt)~d$rep+d$tr_l+d$tr_l:d$rep,
    family = "binomial"),test="LRT")
  return(mod)
}
# FUNCTION: convert array to data.frame
get_dat <- function(array,zeroes=1){
  if(zeroes == 1){
    if(any(array == 0)){
      array<-array+1
    }
  }
  A_Cnt <- vector(length = dim(array)[3]*dim(array)[1])
  Tot_Cnt <- vector(length = dim(array)[3]*dim(array)[1])
  #tr_l <- vector(length = dim(array)[3]*dim(array)[1])
  #rep <- vector(length = dim(array)[3]*dim(array)[1])
  j <-1
  for(k in seq(1,dim(array)[3],1)){
    for(i in seq(1,dim(array)[1],1)){
      #      print(c(i,j,k))
      A_Cnt[j]<-array[i,1,k]
      Tot_Cnt[j]<-sum(array[i,,k])
      #tr_l[j] <- as.character(i)
      #rep[j] <- as.character(k)
      j <- j + 1
    }
  }


  ##CODE FOR LEAH'S GDL POP CAGES 
  ##-----------------------------------------------

  ##hard code treatment ordering based on /data/sync_files.txt 
  ## {CONTROL,ROTENONE} ~ {1,2}
  tr_l = c(rep(c(1, 1, 1, 1, 2, 2, 2, 2), times = 12),rep(c(2,1), times = 12))
  tr_l = as.character(tr_l)

  ##hard code replicate ordering 
  ## {1A1, 1A2, ... , 2B2, 2B3} ~ {1,2,...,11,12}
  rep = c(rep(1:12, each = 8),rep(1:12, each = 2))
  rep = as.character(rep)

  ##hard code time point ordering
  ## {F20,F22,F40,F50} ~ {0,2,20,30}
  time = c(rep(c(2,20,30,40), times = 24), rep(0,24))

  ##-----------------------------------------------

  d<-data.frame("A_Cnt"=A_Cnt,"Tot_Cnt"=Tot_Cnt,"tr_l"=tr_l,"rep"=rep, "time"=time)

  return(d)
}
