########################################################################
# Title: mtbp_functions.R
# Purpose: functions that generate different mean matrices for different networks
# Author: Leah Keating
# Date last edited: 1 July 2021
########################################################################

# newman-type network where every node is part of m 3-cliques
# and n 2-cliques.
mean_mat <- function(alpha, q1, m, n){
  if(m != 0 & n != 0){
    mat <- Matrix(0, nrow = 6, ncol = 6)
    mat[1,1] <- 2*(m-1)*(1-q1)
    mat[1,2] <- 2*q1*(1-q1)
    mat[1,3] <- (1-q1)**2
    mat[1,5] <- 2*n*(1-q1)
    mat[2,1] <- (m-1)*(1-q1*(1-alpha))
    mat[2,4] <- 1-q1*(1-alpha)
    mat[2,5] <- n*(1-q1*(1-alpha))
    mat[5,1] <- m*(1-q1)
    mat[5,5] <- (n-1)*(1-q1)
    mat[5,6] <- 1-q1
  }
  else if(m == 0 & n != 0){
    mat <- Matrix(0, nrow = 2, ncol = 2)
    mat[1,1] <- (n-1)*(1-q1)
    mat[1,2] <- 1-q1
  }
  else if(m!=0 & n==0){
    mat <- Matrix(0, nrow = 4, ncol = 4)
    mat[1,1] <- 2*(m-1)*(1-q1)
    mat[1,2] <- 2*q1*(1-q1)
    mat[1,3] <- (1-q1)**2
    mat[2,1] <- (m-1)*(1-q1*(1-alpha))
    mat[2,4] <- 1-q1*(1-alpha)
  }
  else{
    print("If m=0 and n=0 there can be no diffusion!")
    mat <- 0
  }
  return(mat)
}