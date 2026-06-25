####
# G - test for independence with 2x2xk contingency tables.
# Author: R. Axel W. Wiberg
#
# Notes:
# ALL users of this test should read Sokal and Rohlf (1981) to understand
# what is being done.
# 
# This script also provides a function for converting a 2x2xk contingency table
# into a data set for Binomial GLM.
#
# References:
#
# Sokal and Rohlf (1969, 1981, 2015) Biometry, 
# W.H.Frieman and Co. San Francisco
#
####

## Testing: G-test direct G computation
#####
cat("# Loaded: G_test.R\n")
#test data
###
#matrix <-array(c(1,1,40,42,1,1,35,41,111,1,1,87),dim = c(2,2,3),
#               dimnames = list(c("H","L"),c("A","T"),c("Pop1","Pop2","Pop3")))
###

#test code
###
## L = line
## P = population
## C = allele counts

##number of levels for each table
#l <- dim(matrix)[1]
#c <- dim(matrix)[2]
#p <- dim(matrix)[3]

#grand <- sum(matrix)*log(sum(matrix))

## Get PxL
#PxLa<-array(dim=c(p,l))
#for(i in seq(1,p,1)){
#  for(j in seq(1,l,1)){
#    PxLa[i,j] <- sum(matrix[j,,i])
#  }
#}
#PXLv <- vector(length=(l*p))
#for(i in seq(1,(l*p))){
#  #print(matrix[i]*log(matrix[i]))
#  PXLv[i] <- PxLa[i]*log(PxLa[i])
#}
#PxL <- sum(PXLv)
## Get PxC
#PxCa<-array(dim=c(p,c))
#for(i in seq(1,p,1)){
#  for(j in seq(1,c,1)){
#    PxCa[i,j] <- sum(matrix[,j,i])
#  }
#}
#PXCv <- vector(length=(p*c))
#for(i in seq(1,(p*c))){
#  #print(matrix[i]*log(matrix[i]))
#  PXCv[i] <- PxCa[i]*log(PxCa[i])
#}
#PxC <- sum(PXCv)
## Get LxC
#LXCa<-array(dim=c(l,c))
#for(i in seq(1,l,1)){
#  for(j in seq(1,c,1)){
#    LXCa[i,j]<-sum(matrix[i,j,])  
#  }
#}
#LXCv <- vector(length=(l*c))
#for(i in seq(1,(l*c))){
#  #print(matrix[i]*log(matrix[i]))
#  LXCv[i] <- LXCa[i]*log(LXCa[i])
#}
#LxC <- sum(LXCv)
## Get PxLxC
#PXLXCv <- vector(length=(l*c*p))
#for(i in seq(1,(l*c*p))){
#  #print(matrix[i]*log(matrix[i]))
#  PXLXCv[i] <- matrix[i]*log(matrix[i])
#}
#PxLxC <- sum(PXLXCv) 
#Pv <- vector(length=p)
#for(i in seq(1,p,1)){
#  Pv[i] <- sum(matrix[,,i])*log(sum(matrix[,,i]))
#}
#P <- sum(Pv)
#Lv <- vector(length=l)
#for(i in seq(1,l,1)){
#  Lv[i] <- sum(matrix[i,,])*log(sum(matrix[i,,]))
#}
#L <- sum(Lv)
#Cv <- vector(length=c)
#for(i in seq(1,c,1)){
#  Cv[i] <- sum(matrix[,i,])*log(sum(matrix[,i,]))
#}
#C <- sum(Cv)
#PxLi_G <- 2*(PxL-(P+L)+grand)
#PxCi_G <- 2*(PxC-(P+C)+grand)
#LxCi_G <- 2*(LxC-(L+C)+grand)
#PxLxCi_G <- 2*(PxLxC-(P+L+C)+(2*grand))
#PxLxCx_G <- 2*(PxLxC-(PxL+PxC+LxC)+(P+L+C)-grand)
#PxLi_df <- (p*l)-(p-1)-(l-1)-1 
#PxCi_df <- (p*c)-(p-1)-(c-1)-1 
#LxCi_df <- (l*c)-(l-1)-(c-1)-1
#PxLxCi_G <- (p*l*c)-(p-1)-(l-1)-(c-1)-1
#PxLxCx_G <- (p-1)*(l-1)*(c-1)
#####

## Function: G-test direct G computations
#####
G_test<-function(array,add=FALSE,correction="none"){
  if(add == TRUE & any(array == 0)){
    array<-array+1
  }
  if(correction=="cont"){
    f_hat<-G_test_fhat(array)$f_hat
    for(i in 1:length(array)){
      if(array[i] < f_hat[i]){
        array[i] <- array[i] + 0.5
      }
      if(array[i] > f_hat[i]){
        array[i] <- array[i] - 0.5    
      }
    }
  }
  # L = line
  # P = population
  # C = allele counts
  #number of levels for each table
  l <- dim(array)[1]
  c <- dim(array)[2]
  p <- dim(array)[3]
  grand <- sum(array)*log(sum(array))
  #cat(c("grand",grand,"\n")) # Test line
  # Get PxL
  PxLa<-array(dim=c(p,l))
  for(i in seq(1,p,1)){
    for(j in seq(1,l,1)){
      PxLa[i,j] <- sum(array[j,,i])
    }
  }
  PXLv <- vector(length=(l*p))
  for(i in seq(1,(l*p))){
    #print(array[i]*log(array[i]))
    PXLv[i] <- PxLa[i]*log(PxLa[i])
  }
  PxL <- sum(PXLv)
  #cat(c("PxL",PxL,"\n")) # Test line
  # Get PxC
  PxCa<-array(dim=c(p,c))
  for(i in seq(1,p,1)){
    for(j in seq(1,c,1)){
      PxCa[i,j] <- sum(array[,j,i])
    }
  }
  PXCv <- vector(length=(p*c))
  for(i in seq(1,(p*c))){
    #print(array[i]*log(array[i]))
    PXCv[i] <- PxCa[i]*log(PxCa[i])
  }
  PxC <- sum(PXCv)
  #cat(c("PxC",PxC,"\n")) # Test line
  # Get LxC
  LXCa<-array(dim=c(l,c))
  for(i in seq(1,l,1)){
    for(j in seq(1,c,1)){
      LXCa[i,j]<-sum(array[i,j,])  
    }
  }
  LXCv <- vector(length=(l*c))
  for(i in seq(1,(l*c))){
    #print(array[i]*log(array[i]))
    LXCv[i] <- LXCa[i]*log(LXCa[i])
  }
  LxC <- sum(LXCv)
  #cat(c("LxC",LxC,"\n")) # Test line
  # Get PxLxC
  PXLXCv <- vector(length=(l*c*p))
  for(i in seq(1,(l*c*p))){
    #print(array[i]*log(array[i]))
    PXLXCv[i] <- array[i]*log(array[i])
  }
  PxLxC <- sum(PXLXCv)
  #cat(c("PxLxC",PxLxC,"\n")) # Test line
  Pv <- vector(length=p)
  for(i in seq(1,p,1)){
    Pv[i] <- sum(array[,,i])*log(sum(array[,,i]))
  }
  P <- sum(Pv)
  #cat(c("P",P,"\n")) # Test line
  Lv <- vector(length=l)
  for(i in seq(1,l,1)){
    Lv[i] <- sum(array[i,,])*log(sum(array[i,,]))
  }
  L <- sum(Lv)
  #cat(c("L",L,"\n")) # Test line
  Cv <- vector(length=c)
  for(i in seq(1,c,1)){
    Cv[i] <- sum(array[,i,])*log(sum(array[,i,]))
  }
  C <- sum(Cv)
  #cat(c("C",C,"\n")) # Test line
  PxLi_G <- 2*(PxL-(P+L)+grand)
  PxCi_G <- 2*(PxC-(P+C)+grand)
  LxCi_G <- 2*(LxC-(L+C)+grand)
  PxLxCi_G <- 2*(PxLxC-(P+L+C)+(2*grand))
  PxLxCx_G <- 2*(PxLxC-(PxL+PxC+LxC)+(P+L+C)-grand)
  PxLi_df <- (p*l)-(p-1)-(l-1)-1 
  PxCi_df <- (p*c)-(p-1)-(c-1)-1 
  LxCi_df <- (l*c)-(l-1)-(c-1)-1
  PxLxCi_df <- (p*l*c)-(p-1)-(l-1)-(c-1)-1
  PxLxCx_df <- (p-1)*(l-1)*(c-1)
  PxLi_p <- pchisq(PxLi_G,df=PxLi_df,lower.tail=FALSE) 
  PxCi_p <- pchisq(PxCi_G,df=PxCi_df,lower.tail=FALSE)
  LxCi_p <- pchisq(LxCi_G,df=LxCi_df,lower.tail=FALSE)
  PxLxCi_p <- pchisq(PxLxCi_G,df=PxLxCi_df,lower.tail=FALSE)
  PxLxCx_p <- pchisq(PxLxCx_G,df=PxLxCx_df,lower.tail=FALSE)
  
  printarray <- matrix(c(PxLi_G,PxCi_G,LxCi_G,PxLxCi_G,PxLxCx_G,
                        PxLi_df,PxCi_df,LxCi_df,PxLxCi_df,PxLxCx_df,
                        PxLi_p,PxCi_p,LxCi_p,PxLxCi_p,PxLxCx_p),nrow=5,
                      dimnames=list(c("PxL independence: ","PxC independence: ",
                                 "LxC independence: ","PxLxC independence: ","PxLxC interaction: "),
                                 c("G","df","p-value")))
  return(printarray)
}
#####

## Testing: Get f_hat and f_dev2 matrices
#####
#matrix <- sr_matrix
#K <- dim(matrix)[3]
#K
#L <- dim(matrix)[2]
#L
#A <- dim(matrix)[1]
#A
#n2 <- sum(matrix)^2
#n2
#matrix_hat <- array(dim=dim(matrix))
#for(k in seq(1,K,1)) {
#  for(l in seq(1,L,1)){
#    for(a in seq(1,A,1)){
#      #cat(k,l,a)
#      f<-matrix[l,a,k]
#      f_hat<-(sum(matrix[l,,])*sum(matrix[,a,])*sum(matrix[,,k]))/n2
#      matrix_hat[l,a,k] <- f_hat
#      #f <- matrix3[f]
#    }
#  }
#}
#matrix
#Gv <- vector(length=prod(dim(matrix)))
#for(i in seq(1,prod(dim(matrix)),1)){
#  Gv[i] <- matrix3[i]*(log(matrix[i]/matrix_hat[i]))
#}
#G <- 2*sum(Gv)
#matrix
#matrix_hat
#G_test(matrix)
#####

## Function: Get f_hat and f_dev2 matrices
#####
G_test_fhat<-function(array,add=FALSE){
  if(add == TRUE & any(array == 0)){
    array<-array+1
  }
  K <- dim(array)[3]
  L <- dim(array)[2]
  A <- dim(array)[1]
  n2 <- sum(array)^2
  # Make empty array for the array of f_hat
  array_hat <- array(dim=dim(array),dimnames=dimnames(array))
  # Make empty array for the array of squared deviations from f_hat
  array_dev2 <- array(dim=dim(array),dimnames=dimnames(array))
  # Calculate f_hat and the deviations of f from f_hat
  for(k in seq(1,K,1)) {
    for(l in seq(1,L,1)){
      for(a in seq(1,A,1)){
        #cat(k,l,a)
        f<-array[l,a,k]
        f_hat<-(sum(array[l,,])*sum(array[,a,])*sum(array[,,k]))/n2
        array_hat[l,a,k] <- f_hat
        array_dev2[l,a,k] <- (f-f_hat)^2 
        #f <- matrix3[f]
      }
    }
  }
  Gv <- vector(length=prod(dim(array)))
  for(i in seq(1,prod(dim(array)),1)){
    Gv[i] <- array[i]*(log(array[i]/array_hat[i]))
  }
  df = (K*L*A)-(K-1)-(L-1)-(A-1)-1
  G <- c("G"=2*sum(Gv),"d.f."=df,"p"=pchisq(2*sum(Gv),df,lower.tail = FALSE))
  return(list("G"=G,"f_hat"=array_hat,"f_dev"=array_dev2))
}
#####
#

## Example data from Sokal and Rohlf (1969) multiway independence test.
#####
sr_matrix <- array(c(55,34,6,17,23,15,1,5,7,3,4,5,8,5,3,3),dim=c(2,2,4),
                   dimnames=list(c("Male","Female"),c("Healthy","Poisoned"),
                                 c("IM","AM","OW","OM")))

sr_res<-G_test(sr_matrix,correction="cont")
sr_fhat<-G_test_fhat(sr_matrix)$f_hat
#####

## Testing the functions
#####
#Test 1: Should give interaction:
#matrix <-array(c(1,1,33,90,1,1,36,41,156,1,13,103),dim = c(2,2,3),
#               dimnames = list(c("H","L"),c("A","T"),c("Pop1","Pop2","Pop3")))
#matrix
#g<-G_test(matrix)
#g_hat<-G_test_fhat(matrix)
#
#
#Test 2: Same as Test 1 but total alleles scaled to 40 for each line
#matrix <-array(c(1.18,0.44,38.82,39.56,1.08,0.95,38.92,39.05,36.92,0.38,3.08,39.62),dim = c(2,2,3),
#                dimnames = list(c("H","L"),c("A","T"),c("Pop1","Pop2","Pop3")))
#matrix
#G_test(matrix)
#G_test_fhat(matrix)
#
#Test 3: 
#matrix <-array(c(1.18,0.44,38.82,39.56,1.08,0.95,38.92,39.05),dim = c(2,2,2),
#                dimnames = list(c("H","L"),c("A","T"),c("Pop1","Pop2")))
#matrix
#G_test_fhat(matrix)
#####
#
## Under Construction:
#####
# William's (1976) Correction (pp 745-746, Box 17.8 in Sokal and Rohlf (1981))
# Continuity correction

# For sample sizes n < 25, work out exact probabilities

# Adjusted G values are not additive, so cannot be directly summed or subtracted.
#####