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

# repropduce distribution in Fig. 4.

# get simulations for m = 0, m = 1 and m = 3
# n = 6 - 2n
# do more simulations for a nicer plot (total = 1000000 in paper)

bp_sims_m0 = MTBP_sim(q1 = 0.9, alpha = 0.1, m = 0, n = 6, total = 10000)
bp_sims_m1 = MTBP_sim(q1 = 0.9, alpha = 0.1, m = 1, n = 4, total = 10000)
bp_sims_m3 = MTBP_sim(q1 = 0.9, alpha = 0.1, m = 3, n = 0, total = 10000)

# put this in a data frame

bp_sizes.df <- tibble(size = bp_sims_m0, net = "m = 0, n = 6") %>%
  add_row(tibble(size = bp_sims_m1, net = "m = 1, n = 4")) %>%
  add_row(tibble(size = bp_sims_m3, net = "m = 3, n = 0"))

# get the size distribution for each network

bp_dist.df <- bp_sizes.df %>% filter(size != 1) %>% group_by(size, net) %>%
  # filter out size 1 because these do not show in the network-based simulations
  summarise(n = n()) %>% group_by(net) %>%
  arrange(size) %>% mutate(p = n/sum(n), cdf = cumsum(p), ccdf = 1-cdf)

bp_dist.df %>%
  ggplot(aes(x = size, y = ccdf, colour = net)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  theme(legend.position = c(0.12, 0.2), legend.box.background = element_rect(fill='#ffffff'),
        legend.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size = 11),
        legend.title=element_text(size=14),
        text=element_text(size=16, family="Times New Roman"),
        legend.background = element_blank())

############################################################
# Need to do network-based simulations for comparison
############################################################
# Note: These are much slower

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

tree_like_net <- newman_net(n = 10000, l=0,m=6)
tree_like_adj <- as_adj(tree_like_net)
intermediate_net <- newman_net(n = 10000, l=1,m=4)
intermediate_adj <- as_adj(intermediate_net)
triangle_net <- newman_net(n = 10000, l=3,m=0)
triangle_adj <- as_adj(triangle_net)
four_cl_net <- four_cl_net_fn(n = 10000, m = 2)
four_cl_adj <- as_adj(four_cl_net)

# need a function for simulating the cascades

cascade_sim_par <- function(j, net = tree_like_net, adj = tree_like_adj, q1 = 0.9, alpha = 0.15, total = 1000){
  out.df <- generate_cc_cascades(follower.net = net, follower.adj = adj, q1 = q1, alpha = alpha, total = total)
  out.df <- out.df %>% mutate(ID = str_c(j,ID,sep = "_"))
  return(out.df)
}

cascade_sim_par(j = 1, net = tree_like_net, adj = tree_like_adj, total = 100)
