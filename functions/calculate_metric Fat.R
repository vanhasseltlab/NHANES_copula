# helper functions for the calculation of marginal metrics in Fat imputation analysis
##~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 0 - Load libraries   ----
##~~~~~~~~~~~~~~~~~~~~~~~~

library(tibble)
library(tidyr)
library(ks)
library(sf)
library(dplyr)
library(pracma) 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 1 - calculate the statistics   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simple version for only one variable (Fat)
get_statistics_f <- function(data_set, columns = NULL) {
  
  stats <- function(x) {
    x <- x[!is.na(x)]
    c(mean = mean(x), median = median(x), sd = sd(x), min = min(x), max = max(x),
      Q5 = quantile(x, 0.05), Q95 = quantile(x, 0.95))
  }
  univariate <- as.data.frame(apply(data_set, 2, stats))
  long_format <- univariate %>% 
    rownames_to_column("statistic") %>% 
    pivot_longer(-statistic, names_to = "covariate")
  
  return(long_format)
}

# simple version for only one variable (Fat)
get_statistics_multiple_f <- function(data_set, m, columns = NULL, type = NULL) {
  data_set <- data_set[, colSums(is.na(data_set)) != nrow(data_set)]
  
  if (is.null(columns)) {
    columns <- setdiff(colnames(data_set), c("simulation_nr", "simulation","data_set"))
  }
  n_statistics <- length(columns)*7
  full_results <- as.data.frame(matrix(nrow = m*n_statistics, ncol = 3))
  names(full_results) <- c("statistic", "covariate", "value")
  full_results$simulation_nr <- rep(1:m, each = n_statistics)
  
  
  
  for(i in unique(data_set$simulation_nr)) {
    sim_results <- get_statistics_f(data_set[data_set$simulation_nr == i, ], columns = columns) %>% filter(covariate != "simulation_nr")
    full_results[full_results$simulation_nr == i, c("statistic", "covariate", "value")] <- sim_results
  }
  if (!is.null(type)) {
    full_results$type <- type
  }
  
  return(full_results)
}
