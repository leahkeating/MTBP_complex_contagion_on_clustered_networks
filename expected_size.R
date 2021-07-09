################################################################################
# Title: expected_size.R
# Purpose: to compare the expected size of cascades from simulations and theory
# Author: Leah Keating
# Date last edited: 9 July 2021
################################################################################

source("mtbp_functions.R")
library("boot")

# We initialise a data frame with alpha and q1 values to calculate the
# expected cascade size analyticallly

analytic_size <- numeric()
analytic.df <- expand.grid(alpha = seq(from = 0, to = 0.2, by = 0.001),
                           q1 = c(from = 0.81, 0.83, 0.85)) %>%  as_tibble()
for (i in 1:nrow(analytic.df)) {
  analytic_size <- c(analytic_size, expected_size(m = 3, n = 0, alpha = as.numeric(analytic.df$alpha[i]), q1 = as.numeric(analytic.df$q1[i])))
}

analytic.df <- analytic.df %>% mutate(size = analytic_size, p1 = as.factor(1-q1))

# simulate 1 million cascades for 15 different values and use bootstrapping
# to get the 2.5th and 97.5th percentiles

expected_size.df <- expand.grid(alpha = seq(from = 0, to = 0.2, by = 0.05),
                                q1 = c(0.81, 0.83, 0.85)) %>% as_tibble()

cluster <- makeCluster(6)
registerDoParallel(cluster)

# size_sims holds the values for the size of 1000 cascades of each type
# we will use these 1000 values to find CIs for the mean using bootstrapping
# note 1 million simulation were used for the graphs in the paper
size_sims <- foreach(i=1:nrow(expected_size.df), .combine = "c", .packages = c("igraph", "dplyr", "stringr")) %dopar%
  MTBP_sim(q1 = as.numeric(expected_size.df$q1[i]), alpha = as.numeric(expected_size.df$alpha[i]), m = 3, n = 0, total = 1000)

stopImplicitCluster()

# bootstrapping to find CIs

bs.mean <- function(data, indices){
  return(mean(data[indices]))
}

lq <- numeric()
uq <- numeric()
for (i in 0:14) {
  # change the 1000 below if increasing the number of simulations
  sizes_i <- size_sims[(1+1000*i):(1000*(i+1))]
  boot_mean_i <- boot(sizes_i, statistic = bs.mean, R = 1000, parallel = "multicore")
  CI <- boot.ci(boot_mean_i, conf = 0.95, type = "perc")$percent[4:5]
  lq <- c(lq, CI[1])
  uq <- c(uq, CI[2])
  print(paste0(i+1, " of 15 complete"))
}

# lq is the 2.5th percentile and uq is the 97.5th percentile
expected_size.df <- expected_size.df %>% mutate(lq, uq)

expected_size.df %>% mutate(p1 = 1-q1, p1 = as.factor(p1)) %>%
  ggplot(aes(x = alpha, color = p1)) +
  geom_line(data = analytic.df, aes(y = size)) +
  geom_errorbar(aes(ymin = lq, ymax = uq), colour = "black", width = 0.002) +
  scale_linetype_manual("",values = c("solid","dotted")) +
  scale_shape_manual("",values = c(1,8)) +
  scale_color_brewer(TeX("$p_{1}$"),palette = "Set2") +
  #scale_color_viridis_d(option = "inferno") +
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #             labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks() +
  scale_y_log10() +
  theme(legend.text = element_text(size = 10), legend.position = c(0.08,0.8),
        legend.box.background = element_rect(fill='white'), legend.direction = "horizontal",
        legend.margin = margin(6, 6, 6, 6),
        text=element_text(size=16, family="Times New Roman"),
        legend.title = element_text(size = 11)) +
  ylab("expected cascade size") +
  xlab(TeX("$\\alpha$"))
