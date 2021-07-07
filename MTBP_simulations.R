########################################################################
# Title: MTBP_simulations.R
# Purpose: a simulation scheme for cascades using the MTBP
# Author: Leah Keating
# Date last edited: 7 July 2021
########################################################################

source("mtbp_functions.R")

# make one big function to work for all network types
# MTBP_sim simulates cascades using the MTBP method described in
# the paper

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

