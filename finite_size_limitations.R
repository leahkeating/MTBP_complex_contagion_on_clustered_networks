########################################################################
# Title: finite_size_limitations.R
# Purpose: to reproduce Fig. 6 in paper
# Author: Leah Keating
# Date last edited: 9 July 2021
########################################################################

source("mtbp_functions.R")

# generate networks where each node is in 3 3-cliques
# of size 100, 1000 and 10000
# we also need the adjacency matrix for the simulatiions
net_100 <- newman_net(n = 100, l = 3, m = 0)
adj_100 <- as_adj(net_100)
net_1000 <- newman_net(n = 1000, l = 3, m = 0)
adj_1000 <- as_adj(net_1000)
net_10000 <- newman_net(n = 10000, l = 3, m = 0)
adj_10000 <- as_adj(net_10000)

# set up the simulations in parallel
cluster <- makeCluster(6) # run on 6 cores
registerDoParallel(cluster)

fs_100_sims.df <- foreach(i=1:6, .combine = "rbind", .packages = c("igraph", "dplyr", "stringr")) %dopar%
  cascade_sim_par(i, net = net_100, adj = adj_100, q1 = 0.9, alpha = 0.1, total = 1667)

fs_1000_sims.df <- foreach(i=1:6, .combine = "rbind", .packages = c("igraph", "dplyr", "stringr")) %dopar%
  cascade_sim_par(i, net = net_1000, adj = adj_1000, q1 = 0.9, alpha = 0.1, total = 1667)

fs_10000_sims.df <- foreach(i=1:6, .combine = "rbind", .packages = c("igraph", "dplyr", "stringr")) %dopar%
  cascade_sim_par(i, net = net_10000, adj = adj_10000, q1 = 0.9, alpha = 0.1, total = 1667)

stopImplicitCluster()

# make a data frame with the distribution for each

fs_all_sims.df <- fs_100_sims.df %>% mutate(net_size = 100) %>%
  add_row(fs_1000_sims.df %>% mutate(net_size = 1000)) %>%
  add_row(fs_10000_sims.df %>% mutate(net_size = 10000))

fs_dist.df <- fs_all_sims.df %>% group_by(net_size, ID) %>%
  mutate(size = n_distinct(c(parent, child))) %>%
  ungroup() %>% select(ID, net_size, size) %>%
  unique() %>% group_by(size, net_size) %>%
  summarise(n = n()) %>% group_by(net_size) %>%
  arrange(size) %>%
  mutate(p = n/sum(n), cdf = cumsum(p), ccdf = 1-cdf) %>%
  ungroup() %>% mutate(net_size = as.factor(net_size))

# now we get the distribution for the MTBP simulations

bp_sizes <- MTBP_sim(q1 = 0.9, alpha = 0.1, m = 3, n = 0, total = 1000)

bp_dist.df <- tibble(size = bp_sizes, net_size = "bp") %>%
  filter(size != 1) %>% # filter out size 1 cascades because these are not present in the network-based simulations
  group_by(size) %>%
  summarise(n = n()) %>%
  arrange(size) %>%
  mutate(p = n/sum(n), cdf = cumsum(p), ccdf = 1-cdf)

# make a data frame with both simulation types
dist.df <- fs_dist.df %>% add_row(bp_dist.df %>% mutate(net_size = as.factor("branching process"))) 

# plot Fig. 6 - increase number of simulations above to make it more similar to the paper
dist1.df <- dist.df %>% filter(net_size != "branching process")
dist2.df <- dist.df %>% filter(net_size == "branching process")

require(scales)
ggplot(data = NULL, aes(x = size, y = ccdf)) +
  geom_point(data = dist2.df, shape = 8, colour = "black") +
  geom_line(data = dist1.df, aes(color = net_size)) +
  scale_color_brewer("network size",palette = "Set1")+
  scale_x_log10() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks() +
  theme(legend.text = element_text(size = 10), legend.position = c(0.15,0.35),
        legend.box.background = element_rect(fill='white'),
        text=element_text(size=16, family="Times New Roman"),
        legend.margin = margin(6, 6, 6, 6), legend.title = element_text(size = 11)) +
  xlab("cascade size")
