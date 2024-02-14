# load packages
library(tidyr)
library(ks)
library(sf)
library(dplyr)
library(pracma) 

get_overlap_multiple_a <- function(data_obs, data_set, m, columns = NULL, type = NULL) {
  data_set <- data_set[, colSums(is.na(data_set)) != nrow(data_set)]
  
  if (is.null(columns)) {
    columns <- setdiff(colnames(data_set), c("simulation_nr", "population"))
    # columns <- setdiff(colnames(data_set), "simulation_nr")
  }
  n_statistics <- choose(length(columns), 2)
  full_results <- as.data.frame(matrix(nrow = m*n_statistics, ncol = 5))
  names(full_results) <- c("statistic", "covariate", "value", "cov1", "cov2")
  full_results$simulation_nr <- rep(1:m, each = n_statistics)
  
  
  for(i in unique(data_set$simulation_nr)) {
    overlap_results <- get_overlap_a(obs_data = data_obs, sim_data = data_set[data_set$simulation_nr == i, ],columns = NULL, per =1)
    
    full_results[full_results$simulation_nr == i, c("statistic", "covariate", "value",  "cov1", "cov2")] <-  overlap_results
  }
  if (!is.null(type)) {
    full_results$type <- type
  }
  
  return(full_results)
}

get_overlap_a <- function(obs_data, sim_data, columns = NULL, per = 1) {
  
  if (is.null(columns)) {
    columns <- setdiff(colnames(sim_data), "simulation_nr")
  }
  
  percent <- c("5%","1%","10%")
  # m <- 1 # indicator for percentile
  cor_mat <- cor(obs_data[, columns], use = "pairwise.complete.obs",method = "pearson")
  test <- replace(cor_mat, lower.tri(cor_mat, TRUE), NA) 
  # TURE indicates the diagonal was included, replace the lower triangle to be NA
  cor_data <- reshape2::melt(test, na.rm = TRUE)
  # reshape the data to be a data frame, with all the values in one column
  No <- c(1:nrow(cor_data))
  pair_data <- cor_data[,c(1,2)] %>% 
    mutate(Area_sim = NA) %>% # create the column to store the area of 90th contour for sim data
    cbind(No)
  
  
  for (i in 1 : nrow(pair_data)){
    odata <- obs_data %>% dplyr::select(pair_data[i,1],pair_data[i,2]) %>% na.omit() # extract data
    sdata <- sim_data %>% dplyr::select(pair_data[i,1],pair_data[i,2]) %>% na.omit()
    
    kd_obs <- ks::kde(odata, compute.cont=TRUE) # estimate kernel density
    kd_sim <- ks::kde(sdata, compute.cont=TRUE)
    
    contour_obs <- with(kd_obs, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                             z=estimate, levels=cont[percent[per]]))
    contour_sim <- with(kd_sim, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                             z=estimate, levels=cont[percent[per]])) 
    # to make sure it's closed polygon
    for(b in 1: length(contour_obs)){
      contour_obs[[b]] <- as.data.frame(contour_obs[[b]])
      contour_obs[[b]][nrow(contour_obs[[b]])+1,] <- c(contour_obs[[b]][1,1], contour_obs[[b]][1,2], contour_obs[[b]][1,3])
    }
    contour_obs_long <- bind_rows(contour_obs,.id="circle")
    
    for(b in 1: length(contour_sim)){
      contour_sim[[b]] <- as.data.frame(contour_sim[[b]])
      contour_sim[[b]][nrow(contour_sim[[b]])+1,] <- c(contour_sim[[b]][1,1], contour_sim[[b]][1,2], contour_sim[[b]][1,3])
    }
    contour_sim_long <- bind_rows(contour_sim,.id="circle")
    
    contour_obs_list <- create_polygon(contour_obs_long)
    contour_sim_list <- create_polygon(contour_sim_long)
    
    union_area <- st_union(contour_obs_list,contour_sim_list)
    inter_area <- st_intersection(contour_obs_list,contour_sim_list)
    
    Areaunion = sf::st_area(union_area)
    AreaIntersec = st_area(inter_area)
    
    pair_data[i,3] <- 100 * AreaIntersec/Areaunion  # assign the overlap percentage(over obs area) to the pair_data
    
  }
  overlap <- pair_data %>% 
    mutate(statistic = "overlap", covariate = paste(Var1, Var2, sep = "_")) %>% 
    rename("value"= "Area_sim") %>% 
    dplyr::select(-No) %>% 
    mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>% 
    dplyr::select(statistic, covariate, value, Var1, Var2)
  return(overlap)
}


# original version wrote by Laura Zwep
create_polygon <- function(ci_data) {
  polygon_list <- list()
  for (cr in unique(ci_data$circle)) {
    dat_matrix <- list(ci_data %>% filter(circle == cr) %>% select(x, y) %>% as.matrix)
    polygon_list[[as.character(cr)]] <- sf::st_polygon(dat_matrix)
  }
  
  if (length(polygon_list) >= 2) {
    p_union <- polygon_list[[1]]
    for (cr in 2:length(polygon_list)) {
      p_union <- st_union(p_union, polygon_list[[cr]])
    }
    
    polygon_pairs <- t(combn(1:length(polygon_list), 2))
    p_intersect <- list()
    for (cr in 1:nrow(polygon_pairs)) {
      p_intersect[[cr]] <- st_intersection(polygon_list[[polygon_pairs[cr, 1]]], polygon_list[[polygon_pairs[cr, 2]]])
    }
    
    if (length(p_intersect) > 1){
      p_intersect_full <- p_intersect[[1]]
      for (cr in 2:length(p_intersect)) {
        p_intersect_full <- st_union(p_intersect_full, p_intersect[[cr]])
      }
    } else {
      p_intersect_full <- p_intersect[[1]]
    }
    
    p_full <- st_difference(p_union, p_intersect_full)
    
  } else {
    p_full <- polygon_list[[1]]
  }
  
  return(p_full)
}

