########################################################################
# Title: mtbp_functions.R
# Purpose: functions that generate different mean matrices for different networks
# Author: Leah Keating
# Date last edited: 1 July 2021
########################################################################
library(Matrix)

##########################
# Generate the mean matrix
##########################

# alpha is the social reinforcement parameter
# q1 is the probability of not adopting after 1 exposure
# m is the number of triangles that each node is part of
# n is the number of single links that each node is part of
# cl_4 is the number of 4-cliques that each node is part of
# this code doesn't allow us to produce matrices which mix
# 4-cliques with other clique types

mean_mat <- function(alpha, q1, m, n, cl_4 = 0){
  if(cl_4 == 0){
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
  }else if(cl_4 != 0){
    if (m == 0 & n == 0){
      mat <- Matrix(0, nrow = 7, ncol = 7)
      mat[1,1] <- (cl_4-1)*(3*(1-q1)*(q1**2) +6*((1-q1)**2)*q1 + 3*((1-q1)**3))
      mat[1,2] <- 3*(1-q1)*(q1**2)
      mat[1,3] <- 3*((1-q1)**2)*q1
      mat[1,4] <- (1-q1)**3
      mat[2,1] <- (cl_4-1)*(2*q1*(1-alpha)*(1-q1*(1-alpha)) + 2*((1-q1*(1-alpha))**2))
      mat[2,5] <- 2*q1*(1-alpha)*(1-q1*(1-alpha))
      mat[2,6] <- (1-q1*(1-alpha))**2
      mat[3,1] <- (cl_4-1)*(1-((q1**2)*((1-alpha)**3)))
      mat[3,7] <- 1-(q1**2)*((1-alpha)**3)
      mat[5,1] <- (cl_4-1)*(1-q1*((1-alpha)**2))
      mat[5,7] <- 1-q1*((1-alpha)**2)
    }else{
      mat <- 0
      print("Invalid values, can't find matrix if cl_4 > 0 and m or n > 0")
    }
  }
  return(mat)
}

