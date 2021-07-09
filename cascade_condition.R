########################################################################
# Title: cascade_condition.R
# Purpose: code to determine parameters which result in sub- and supercritical cascades
# Author: Leah Keating
# Date last edited: 7 July 2021
########################################################################

# load the functions - we need the mean_mat function in particular
source("mtbp_functions.R")

# initialise a data frame with the alpha, q1 values we are interested in
# m is the number of triangles, n is the number of single edges, cl_4 is the
# number of 4-cliques

criticality.df <- expand.grid(alpha = seq(from = 0, to = 1, by = 0.01),
            q1 = seq(from = 0.7, to = 0.9, by = 0.001)) %>%
  as_tibble()

# function to calculate the criticality
# returns True if subcritical, false if supercritical

criticality <- function(alpha, q1, m, n, cl_4 = 0){
  x <- mean_mat(alpha = alpha, q1 = q1, m = m, n = n, cl_4 = cl_4)
  lambda <- Re(eigen(x)$values) %>% max()
  if (lambda>1){
    return(FALSE)
  }else{
    return(TRUE)
  }
}

# next we want to calculate the criticality for each network type
# for all parameter combinations

# do this in parallel

cluster <- makeCluster(6) # run on 6 cores
registerDoParallel(cluster)

# 6 2-cliques

crit_6_2cl <- foreach(i=1:nrow(criticality.df), .combine = "c", .packages = c("Matrix", "dplyr")) %dopar%
  criticality(alpha = criticality.df$alpha[i], q1 = criticality.df$q1[i], m = 0, n = 6, cl_4 = 0)

# m = 1, n = 4

crit_4_2cl_1_3cl <- foreach(i=1:nrow(criticality.df), .combine = "c", .packages = c("Matrix", "dplyr")) %dopar%
  criticality(alpha = criticality.df$alpha[i], q1 = criticality.df$q1[i], m = 1, n = 4, cl_4 = 0)

# 3 3-cliques

crit_3_3cl <- foreach(i=1:nrow(criticality.df), .combine = "c", .packages = c("Matrix", "dplyr")) %dopar%
  criticality(alpha = criticality.df$alpha[i], q1 = criticality.df$q1[i], m = 3, n = 0, cl_4 = 0)

# 2 4-cliques

crit_2_4cl <- foreach(i=1:nrow(criticality.df), .combine = "c", .packages = c("Matrix", "dplyr")) %dopar%
  criticality(alpha = criticality.df$alpha[i], q1 = criticality.df$q1[i], m = 0, n = 0, cl_4 = 2)

stopImplicitCluster()

# create a data frame where for each network True indicates that those parameters lead to subcritical
# diffusion and False means that there is supercritical diffusion

criticality.df <- criticality.df %>% mutate(crit_6_2cl, crit_4_2cl_1_3cl, crit_3_3cl, crit_2_4cl)

# next we would like to plot the boundaries between the sub- and supercritical boundaries,
# we can do this by finding the minimum alpha value for each q1 for which we get subcritical
# cascades

# add in p1 = 1-q1

criticality.df <- criticality.df %>% mutate(p1 = 1-q1)

# for 6 2-cliques the critical transition occurs at p1 = 0.2
criticality.df %>% filter(crit_6_2cl == TRUE) %>%
  group_by(p1) %>% summarise(alpha = max(alpha))
# this is a straight line at p1 = 0.2 and we can't really capture it 
# this way so we'll do it another way

critical_boundaries.df <- criticality.df %>% filter(crit_4_2cl_1_3cl == TRUE) %>%
  group_by(p1) %>% summarise(alpha = max(alpha)) %>%
  mutate(net = "m = 1, n  = 4") %>%
  add_row(criticality.df %>% filter(crit_3_3cl == TRUE) %>%
            group_by(p1) %>% summarise(alpha = max(alpha)) %>%
            mutate(net = "m = 3, n = 0")) %>%
  add_row(criticality.df %>% filter(crit_2_4cl == TRUE) %>%
            group_by(p1) %>% summarise(alpha = max(alpha)) %>%
            mutate(net = "2 4-cliques")) %>%
  add_row(criticality.df %>% select(alpha) %>% mutate(p1 = 0.2, net = "m = 0, n = 6"))

# the following code removes the alpha = 1 entries for all except those
# which are on the critical boundary

critical_boundaries.df <- critical_boundaries.df %>%
  filter(alpha != 1) %>%
  add_row(critical_boundaries.df %>% filter(alpha == 1) %>%
            group_by(alpha, net) %>% summarise(p1 = max(p1)))

# plot this to see what it looks like
# this gives us a line for each network

critical_boundaries.df %>%
  ggplot(aes(x = p1, y = alpha, linetype = net)) +
  geom_line() +
  xlab(TeX("$p_{1}$")) +
  ylab(TeX("$\\alpha$")) +
  theme(legend.position = c(0.12, 0.2), legend.box.background = element_rect(fill='#ffffff'),
                              legend.margin = margin(6, 6, 6, 6),
                              legend.text = element_text(size = 11),
                              legend.title=element_text(size=14),
                              text=element_text(size=16, family="Times New Roman"),
                              legend.background = element_blank())



