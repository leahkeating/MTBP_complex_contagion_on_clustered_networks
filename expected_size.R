################################################################################
# Title: expected_size.R
# Purpose: to compare the expected size of cascades from simulations and theory
# Author: Leah Keating
# Date last edited: 8 July 2021
################################################################################

source("mtbp_functions.R")

# We initialise a data frame with alpha and q1 values to calculate the
# expected cascade size analyticallly

analytic_size <- numeric()
analytic.df <- expand.grid(alpha = seq(from = 0, to = 0.2, by = 0.001),
                           q1 = c(from = 0.81, 0.83, 0.85)) %>%  as_tibble()
for (i in 1:nrow(analytic.df)) {
  analytic_size <- c(analytic_size, expected_size(m = 3, n = 0, alpha = as.numeric(analytic.df$alpha[i]), q1 = as.numeric(analytic.df$q1[i])))
}

analytic.df <- analytic.df %>% mutate(size = analytic_size)

