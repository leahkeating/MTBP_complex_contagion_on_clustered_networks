########################################################################
# Title: mtbp_functions.R
# Purpose: functions that generate different mean matrices for different networks
# Author: Leah Keating
# Date last edited: 6 July 2021
########################################################################
library(Matrix)
library(tidyverse)
library(igraph)
library(doParallel)
library(cowplot)
theme_set(theme_cowplot())
library(latex2exp)

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

####################################
# Expected cascade size - analytic
####################################

# as described in Sec. IV A this function calculates the expected caascade size
# analytically using the mean matrix

expected_size <- function(alpha, q1, m, n, cl_4 = 0){
  if (cl_4 == 0){
    if (m!= 0 & n != 0 ){
      if(2*m + n != 6){print("This function only works for 6-regular networks.")}
      M <- mean_mat(alpha = alpha, q1 = q1, m = 1, n = 4)
      a <- c(0,1,2,1,0,1)
      z0 <- c(m,0,0,0,n,0) %>% Matrix()
      x_tilde <- Matrix::solve(diag(nrow(M))-M,Matrix(a))
      x <- 1 + t(z0)%*%  M%*%x_tilde
      return(as.numeric(x))
    }else if(n == 0){
      if(m!=3){print("This function only works for 6-regular networks.")}
      M <- mean_mat(q1 = q1, alpha = alpha, m = m, n = n)
      a <- c(0,1,2,1)
      z0 <- c(3,0,0,0) %>% Matrix()
      x_tilde <- Matrix::solve(diag(nrow(M))-M,Matrix(a))
      x <- 1 + t(z0)%*%  M%*%x_tilde
      return(as.numeric(x))
    }else if(m==0){
      if(n!=6){print("This function only works for 6-regular networks.")}
      p1 <- 1-q1
      x_offspring <- 1/(1-5*p1)
      x_seed <- 1+6*p1*x_offspring
      return(x_seed)
    }
  }else if(cl_4 == 2){
    if((m == 0 & n == 0)==FALSE){print("This function only works for 6-regular networks.")}
    M <- mean_mat(q1 = q1, alpha = alpha, m = 0, n = 0, cl_4 = cl_4)
    a <- c(0,1,2,3,1,2,1)
    z0 <- c(2,0,0,0,0,0,0) %>% Matrix()
    x_tilde <- Matrix::solve(diag(nrow(M))-M,Matrix(a))
    x <- 1 + t(z0)%*%  M%*%x_tilde
    return(as.numeric(x))
  }else{
    print("invalid input")
  }
}

#######################################################
# Branching process simulation code
#######################################################
# Multi-type BP simulations for triangle type network
#######################################################

# simulate result matrix
m_triangle_bp_sim <- function(q1 = 0.85, alpha = 0, m = 2){
  # function to test the CC BP for cliques
  # matrix: population of each z1, z2 ,z3 ,z4
  # q1: prob of a node note adopting 
  # alpha: social reinf
  br_res=matrix(NA, nrow = 10000, 4)
  colnames(br_res) <- c('x1','x2','x3','x4')
  matrix = br_res
  matrix[1,] <- c(m,0,0,0)
  for(i in 2:nrow(matrix)){ 
    # im using a for loop because I've preallocated the size of the matrix
    # should be a little but more time efficients
    z_curr <- matrix[i-1,] # the previous steps population
    z_new <- c(0,0,0,0) # vector for offspring
    
    z1_res <- rbinom(n = z_curr[1], size = 2, 1 - q1) # number off spring from z1
    z2_res <- rbinom(n = z_curr[2], size = 1, 1 - (q1 * (1 - alpha))) # number of off spring from z2
    
    # state of the z_vec in the next timestep
    z_new[1] <- z_new[1] + (m-1)*sum(z1_res == 1) + (2 *(m-1)*sum(z1_res == 2)) + (m-1)*sum(z2_res == 1)
    z_new[2] <- z_new[2] + sum(z1_res == 1)
    z_new[3] <- z_new[3] + sum(z1_res == 2)
    z_new[4] <- z_new[4] + sum(z2_res == 1)
    
    matrix[i,] <- z_new # save to current z_new to maxtrix
    
    # stop if we do not have any new offspring
    if(z_new[1] == 0 & z_new[2] == 0) break
  }
  
  # change to tibble and remove rows that we do not use
  matrix <- matrix %>% as_tibble() %>% filter(complete.cases(.))
  
  return(matrix)
}
# function for finding its size
m_triangle_size <- function(res) {
  totals <- colSums(res)
  return(sum(1,totals[2],2*totals[3],totals[4]))
}
# simulate many times with the result being a vector of the sizes
m_triangle_sim_sizes <- function(q1, alpha, m, total = 100){
  sizes <- numeric()
  for (i in 1:total) {
    sizes <- c(sizes,m_triangle_bp_sim(q1 = q1, alpha = alpha, m = m) %>% m_triangle_size())
  }
  return(sizes)
}
#m_triangle_sim_sizes(q1 = 0.9, alpha = 0.3, m = 3)

#######################################################
# Multi-type BP simulations for tree type network
#######################################################

l_edges_bp_sim <- function(q1 = 0.85, l = 6){
  br_res=matrix(NA, nrow = 10000, 2)
  colnames(br_res) <- c('x1','x2')
  matrix = br_res
  matrix[1,] <- c(l,0)
  for(i in 2:nrow(matrix)){ 
    # im using a for loop because I've preallocated the size of the matrix
    # should be a little but more time efficients
    z_curr <- matrix[i-1,] # the previous steps population
    z_new <- c(0,0) # vector for offspring
    
    z1_res <- rbinom(n = z_curr[1], size = 1, 1 - q1) # number off spring from z1
    
    # state of the z_vec in the next timestep
    z_new[1] <- z_new[1] + (l-1)*sum(z1_res == 1)
    z_new[2] <- z_new[2] + sum(z1_res == 1)
    
    matrix[i,] <- z_new # save to current z_new to maxtrix
    
    # stop if we do not have any new offspring
    if(z_new[1] == 0) break
  }
  
  # change to tibble and remove rows that we do not use
  matrix <- matrix %>% as_tibble() %>% filter(complete.cases(.))
  
  return(matrix)
}
l_edges_size <- function(res){
  totals <- colSums(res)
  return(sum(1,sum(totals[2])))
}
l_edges_sim_sizes <- function(q1, l, total = 100){
  sizes <- numeric()
  for (i in 1:total) {
    sizes <- c(sizes,l_edges_bp_sim(q1 = q1, l = l) %>% l_edges_size())
  }
  return(sizes)
}
#l_edges_sim_sizes(q1 = 0.9, l = 6)

###########################################################
# Multi-type BP simulations for m triangles l single edges
###########################################################

newman_bp_sim <- function(q1 = 0.85, alpha = 0, m = 1, l = 4){# m 3-cliques, l 2-cliques
  # function to test the CC BP for cliques
  # matrix: population of each z1, z2 ,z3 ,z4
  # q1: prob of a node note adopting 
  # alpha: social reinf
  br_res=matrix(NA, nrow = 10000, 6)
  colnames(br_res) <- c('x1','x2','x3','x4','x5','x6')
  matrix = br_res
  matrix[1,] <- c(m,0,0,0,l,0)
  for(i in 2:nrow(matrix)){ 
    # im using a for loop because I've preallocated the size of the matrix
    # should be a little but more time efficients
    z_curr <- matrix[i-1,] # the previous steps population
    z_new <- c(0,0,0,0,0,0) # vector for offspring
    
    z1_res <- rbinom(n = z_curr[1], size = 2, 1 - q1) # number off spring from z1
    z2_res <- rbinom(n = z_curr[2], size = 1, 1 - (q1 * (1 - alpha))) # number of off spring from z2
    z5_res <- rbinom(n = z_curr[5], size = 1, 1 - q1) # number of off spring from z5
    
    # state of the z_vec in the next timestep
    z_new[1] <- z_new[1] + (m-1)*sum(z1_res == 1) + (2 *(m-1)*sum(z1_res == 2)) + (m-1)*sum(z2_res == 1) + m*sum(z5_res == 1)
    z_new[2] <- z_new[2] + sum(z1_res == 1)
    z_new[3] <- z_new[3] + sum(z1_res == 2)
    z_new[4] <- z_new[4] + sum(z2_res == 1)
    z_new[5] <- z_new[5] + l*sum(z1_res == 1) + 2*l*sum(z1_res == 2) + l*sum(z2_res == 1) + (l-1)*sum(z5_res == 1)
    z_new[6] <- z_new[6] + sum(z5_res == 1)
    
    matrix[i,] <- z_new # save to current z_new to maxtrix
    
    # stop if we do not have any new offspring
    if(z_new[1] == 0 & z_new[2] == 0 & z_new[5] == 0) break
  }
  # change to tibble and remove rows that we do not use
  matrix <- matrix %>% as_tibble() %>% filter(complete.cases(.))
  return(matrix)
}
newman_size <- function(matrix){
  total_arr <- colSums(matrix)
  return(sum(1,total_arr[2],2*total_arr[3],total_arr[4],total_arr[6]))
}
newman_sim_sizes <- function(q1, alpha, m, l, total = 100){
  sizes <- numeric()
  for (i in 1:total) {
    sizes <- c(sizes,newman_bp_sim(q1 = q1, alpha = alpha, m = m, l = l) %>% newman_size())
  }
  return(sizes)
}
#newman_sim_sizes(q1 = 0.85, alpha = 0, m = 1, l = 4)

###########################################################
# Multi-type BP simulations for m 4-cliques
###########################################################

four_cl_bp_sim <- function(q1 = 0.85, alpha = 0, m = 2){
  # function to test the CC BP for cliques
  # matrix: population of each z1, z2 ,z3 ,z4
  # q1: prob of a node note adopting 
  # alpha: social reinf
  br_res_4cl=matrix(NA, nrow = 10000, 7)
  colnames(br_res_4cl) <- c('x1','x2','x3','x4','x5','x6','x7')
  matrix <- br_res_4cl
  matrix[1,] <- c(m,0,0,0,0,0,0)
  for(i in 2:nrow(matrix)){ 
    # im using a for loop because I've preallocated the size of the matrix
    # should be a little but more time efficients
    z_curr <- matrix[i-1,] # the previous steps population
    z_new <- c(0,0,0,0,0,0,0) # vector for offspring
    
    z1_res <- rbinom(n = z_curr[1], size = 3, 1 - q1) # number off spring from z1
    z2_res <- rbinom(n = z_curr[2], size = 2, 1 - (q1 * (1 - alpha))) # number of off spring from z2
    z3_res <- rbinom(n = z_curr[3], size = 1, 1-((q1**2)*((1-alpha)**3)))
    z5_res <- rbinom(n = z_curr[5], size = 1, 1-(q1*((1-alpha)**2)))
    
    # state of the z_vec in the next timestep
    z_new[1] <- z_new[1] + (m-1)*(sum(z1_res == 1) + (2 *sum(z1_res == 2)) + (3*sum(z1_res == 3)) + sum(z2_res == 1) + 2*sum(z2_res == 2) + sum(z3_res == 1) + sum(z5_res == 1))
    z_new[2] <- z_new[2] + sum(z1_res == 1)
    z_new[3] <- z_new[3] + sum(z1_res == 2)
    z_new[4] <- z_new[4] + sum(z1_res == 3)
    z_new[5] <- z_new[5] + sum(z2_res == 1)
    z_new[6] <- z_new[6] + sum(z2_res == 2)
    z_new[7] <- z_new[7] + sum(z3_res == 1) + sum(z5_res == 1)
    
    matrix[i,] <- z_new # save to current z_new to maxtrix
    
    # stop if we do not have any new offspring
    if(z_new[1] == 0 & z_new[2] == 0 & z_new[3] == 0 & z_new[5] == 0) break
  }
  
  # change to tibble and remove rows that we do not use
  matrix <- matrix %>% as_tibble() %>% filter(complete.cases(.))
  
  return(matrix)
}
four_cl_size <- function(matrix){
  total_arr <- colSums(matrix)
  return(sum(1,total_arr[2],2*total_arr[3],3*total_arr[4],total_arr[5], 2*total_arr[6], total_arr[7]))
} 
four_cl_sim_sizes <- function(q1, alpha, m, total = 100){
  sizes <- numeric()
  for (i in 1:total) {
    sizes <- c(sizes,four_cl_bp_sim(q1 = q1, alpha = alpha, m = m) %>% four_cl_size())
  }
  return(sizes)
}
#four_cl_sim_sizes(q1 = 0.85, alpha = 0, m = 2)

# simulation function that works for all 4 network types

MTBP_sim <- function(q1, alpha = 0, m, n, cl_4 = 0, total = 100){
  if(cl_4 == 0){
    if(m == 0 & n == 6){
      return(l_edges_sim_sizes(q1 = q1, l = 6, total = total))
    }else if(m == 3 & n == 0){
      return(m_triangle_sim_sizes(q1 = q1, alpha = alpha, m = 3, total = total))
    }else if(2*m + n == 6){
      return(newman_sim_sizes(q1 = q1, alpha = alpha, m = m, l = n, total = total))
    }else{
      print("The degree of the nodes must be 6.")
    }
  }else if(cl_4 == 2){
    if(m ==0 & n == 0){
      return(four_cl_sim_sizes(q1 = q1, alpha = alpha, m = 2, total = total))
    }else{
      print("Please do not mix 4-cliques with other clique sizes.")
    }
  }else{
    print("Invalid input.")
  }
}

############################################################
# Functions for the network-based simulations
############################################################

projection_from_bipartite_df <- function(vertex_clique.df, nodes){
  bipartite.net <- graph_from_data_frame(vertex_clique.df, directed = FALSE)
  V(bipartite.net)$type <- V(bipartite.net)$name %in% nodes
  return(bipartite.projection(bipartite.net)$proj2)
}

# function which generates network with l triangles and m links per node
newman_net <- function(n = 1000, l = 0, m = 6){
  if(l!=0 & m!=0){
    while (!((n%%l==0) & (n%%m==0) & (n%%3==0) & (n%%2==0))) {
      n <- n+1
    }
    nodes <- str_c("v_",1:n, sep="")
    three_cliques <- str_c("three_",1:(l*(n/3)), sep = "") %>% rep(., each = 3) %>% sample()
    two_cliques <- str_c("two_",1:(m*(n/2)), sep = "") %>% rep(., each = 2) %>% sample()
    vertex_clique.df <- tibble(vertex = rep(nodes,l+m), clique = c(three_cliques, two_cliques))
  }else if(l != 0){
    while (!((n%%l==0) & (n%%3==0))) {
      n <- n+1
    }
    nodes <- str_c("v_",1:n, sep="")
    three_cliques <- str_c("three_",1:(l*(n/3)), sep = "") %>% rep(., each = 3) %>% sample()
    vertex_clique.df <- tibble(vertex = rep(nodes,l), clique = three_cliques)
  }else{
    while (!((n%%m==0) & (n%%2==0))) {
      n <- n+1
    }
    nodes <- str_c("v_",1:n, sep="")
    two_cliques <- str_c("two_",1:(m*(n/2)), sep = "") %>% rep(., each = 2) %>% sample()
    vertex_clique.df <- tibble(vertex = rep(nodes,m), clique = two_cliques)
  }
  net <- projection_from_bipartite_df(vertex_clique.df, nodes)
  return(net)
}
# function that generates a network where each node is in m 4-cliques
four_cl_net_fn <- function(n = 1000, m = 2){
  while (!((n%%m==0)& (n%%4==0) & (n%%2==0))) {
    n <- n+1
  }
  nodes <- str_c("v_",1:n, sep="")
  four_cliques <- str_c("four_",1:(m*(n/4)), sep = "") %>% rep(., each = 4) %>% sample()
  vertex_clique.df <- tibble(vertex = rep(nodes,m), clique = four_cliques)
  net <- projection_from_bipartite_df(vertex_clique.df, nodes)
  return(net)
}

# complex contagion model
p_cc_threshold <- function(q1 = 0.998, alpha, n, m){
  vals <- numeric(length(n))
  vals[which(n<m)] <- 1-q1*(1-alpha[which(n<m)])**(n[which(n<m)]-1)
  vals[which(n>=m)] <- 1-q1*(1-alpha[which(n>=m)])**(m-1)
  return(vals)
}

# function for generating the cascades
generate_cc_cascades <- function(follower.net = follower.net, follower.adj = follower.adj, q1 = 0.998, alpha = 0.002, total = 1000){#, seed_ = NULL
  all_cascades.df <- tibble(parent = character(), child = character(), generation = numeric(), ID = numeric(), exposures = numeric())
  for (j in 1:total) {
    active <- numeric()
    inactive <- numeric()
    removed <- numeric()
    
    vertex_names <- vertex_attr(follower.net)$name
    seed <- sample(vertex_names,1)
    active <- seed
    inactive <- vertex_names[! vertex_names %in% seed]
    
    exposures <- numeric(gorder(follower.net))
    names(exposures) <- vertex_names
    cascade.df <- tibble(parent = character(), child = character(), generation = numeric())
    generation <- 1
    while (length(active)>0) {
      new_active <- character()
      # shuffle active
      if(length(active)>1){
        active <- sample(active)
      }
      for (i in active) {
        followers <- vertex_names[which(follower.adj[,i]==1)]
        potential_adopters <- followers[followers %in% inactive]
        exposures[potential_adopters] <- exposures[potential_adopters] + 1
        if(length(potential_adopters)>0){
          # fix this, problem is with having n a vector
          adopters <- potential_adopters[runif(length(potential_adopters)) < p_cc_threshold(q1 = q1, alpha = rep(alpha, length(potential_adopters)), n = exposures[potential_adopters], m = 5)]
          if(length(adopters)>0){
            new_active <- c(new_active, adopters)
            inactive <- inactive[! inactive %in% new_active]
            cascade.df <- cascade.df %>% add_row(parent = rep(i, length(adopters)), child = adopters, generation = rep(generation, length(adopters)))
          }
        }
      }
      generation <- generation + 1
      removed <- c(removed, active)
      active <- new_active
    }
    if(nrow(cascade.df)>0){
      all_cascades.df <- all_cascades.df %>% add_row(cascade.df %>% mutate(ID = rep(j, nrow(cascade.df)), exposures = exposures[cascade.df$child]))
    }
  }
  return(all_cascades.df)
}

# function for simulating the cascades in parallel
cascade_sim_par <- function(j, net = tree_like_net, adj = tree_like_adj, q1 = 0.9, alpha = 0.15, total = 1000){
  out.df <- generate_cc_cascades(follower.net = net, follower.adj = adj, q1 = q1, alpha = alpha, total = total)
  out.df <- out.df %>% mutate(ID = str_c(j,ID,sep = "_"))
  return(out.df)
}

# this is the simulation function that is also in MTBP_simulations.R

MTBP_sim <- function(q1, alpha = 0, m, n, cl_4 = 0, total = 100){
  if(cl_4 == 0){
    if(m == 0 & n == 6){
      return(l_edges_sim_sizes(q1 = q1, l = 6, total = total))
    }else if(m == 3 & n == 0){
      return(m_triangle_sim_sizes(q1 = q1, alpha = alpha, m = 3, total = total))
    }else if(2*m + n == 6){
      return(newman_sim_sizes(q1 = q1, alpha = alpha, m = m, l = n, total = total))
    }else{
      print("The degree of the nodes must be 6.")
    }
  }else if(cl_4 == 2){
    if(m ==0 & n == 0){
      return(four_cl_sim_sizes(q1 = q1, alpha = alpha, m = 2, total = total))
    }else{
      print("Please do not mix 4-cliques with other clique sizes.")
    }
  }else{
    print("Invalid input.")
  }
}