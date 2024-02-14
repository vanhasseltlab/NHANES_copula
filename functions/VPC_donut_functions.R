#VPC donut simulation function, geom_creation function and plot method

#depends on:
# dplyr
# ggplot
# ks
# sf


#' Restructure calculated contours in data.frame
#'
#' Output from contourLines is converted from a list to a data.frame for 
#' further use in other functions.
#'
#' @param contour_list the output of the contourLInes functoin
#' @param contour_levels numeric vector with density levels at which the 
#' percentiles are calculated
#' @param sim_nr value for identifier which simulation the contours are 
#' calculated, can be specified
#' @param pair vector with two strings with the names of the covariates
#'
#' @return data.frame with coordinates for a certain density containing a 
#' percentile of the observations
#'

extract_contour_df <- function(contour_list, contour_levels, sim_nr = NULL, pair = NULL) {
  cont_df <- NULL
  for (i in 1:length(contour_list)) {
    cont_df <- rbind.data.frame(cont_df,
                                contour_list[[i]] %>% as.data.frame() %>% 
                                  mutate(circle = i, percentile = names(contour_levels[contour_levels == level[1]])))
  }
  if (!is.null(sim_nr)) {
    cont_df <- cont_df %>% mutate(sim_nr = sim_nr)
  }
  if (!is.null(pair)) {
    cont_df <- cont_df %>% mutate(var1 = pair[1], var2 = pair[2])
  }
  
  return(cont_df)
}

#dependency function: create_polygon
#' Extract the polygons from contour calculations
#' 
#' Create a sf polygon object from (multiple) contour(s)
#'
#'
#' @param ci_data output from extract_contour function, data.frame with 
#' coordinates for a certain density containing a percentile of the observations
#' 
#' @return sf polygon object with the polygon(s) of a contour
#' 

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

# in case some of the polygons are not closed
# force them to be closed

create_polygon_c <- function(ci_data) {
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

#' Get contours from simulated data
#'
#' Simulate or add simulated values from multiple simulations to compute
#' the contours for a donut VPC
#'
#' @param vine an object of class "vine_dist", from rvinecopulib
#' @param percentiles (a vector of) numeric value(s) representing the 
#' percentiles of the density distribution that appear in the VPC donut plot
#' @param B an integer indication number of simulations
#' @param pairs_matrix matrix with 2 column and each row containing a pair of 
#' covariate names from the copula. If set to NULL, every possible covariate 
#' pair is included 
#' @param nobs number of observations for each simulation, if set to NULL, the 
#' number of observations from the original data will be picked
#'
#' @return data.frame with contours from the different simulations over all 
#' percentiles 
#' @export
#'
#' @examples
#' # specify pair-copulas
#' bicop <- bicop_dist("bb1", 90, c(3, 2))
#' pcs <- list(
#'   list(bicop, bicop), # pair-copulas in first tree
#'   list(bicop) # pair-copulas in second tree
#' )
#' # specify R-vine matrix
#' mat <- matrix(c(1, 2, 3, 1, 2, 0, 1, 0, 0), 3, 3)
#' # set up vine copula model
#' vc <- vine_dist(list(distr = "norm"), pcs, mat)
#' 
#' sim_contour_data <- simulate_contours(vc, c(10, 50, 90), B = 10, nobs = 10)
#'


simulate_contours <- function(vine, percentiles, B, pairs_matrix = NULL, nobs = NULL) {
  sim_contours_list <- list()
  
  if (is.null(vine$names)) {
    vine$names <- paste0("V", 1:vine$copula$structure$d)
  }

  if (is.numeric(nobs)) {
    vine$nobs <- nobs
  } else if (is.null(vine$nobs)) {
    stop("nobs should be specified if the vine is not estimated from data",
         " and does not contain nobs in the vine object (vine$nobs)")
  }
  
  if (is.null(pairs_matrix)) {
    pairs_matrix <- t(combn(vine$names, 2))
  } else if (!all(c(pairs_matrix) %in% c(vine$names))) {
    #check for names not matching between vine and pairs_matrix
    stop("Covariate names in pairs_matrix differ from names in vine")
  }
  
  #simulate from vine
  full_sim_data <- as.data.frame(rvine(vine$nobs*B, vine))
  full_sim_data$b <- rep(1:B, each = vine$nobs)
  
  #calculate contours for every percentile and every covariate combination
  i <- 1
  for (b in 1:B) {
    cat("\r", "Calculate contours simulation:", b, "/", B)
    sim_data_b <- full_sim_data[full_sim_data$b == b, ]
    for (p in 1:nrow(pairs_matrix)) {
      #use ks for density computation
      kd_sim <- ks::kde(sim_data_b[, pairs_matrix[p, ]], compute.cont = TRUE, approx.cont = FALSE)
      contour_sim <- with(kd_sim, contourLines(x = eval.points[[1]], y = eval.points[[2]],
                                               z = estimate, levels = cont[paste0(percentiles, "%")]))
      #extract information
      sim_contours_list[[i]] <- extract_contour_df(contour_sim, kd_sim$cont, b, pairs_matrix[p, ])
      i <- i + 1
      
    }

  }
  cat("\n")
  sim_contours <- dplyr::bind_rows(sim_contours_list)
  
  return(sim_contours)
}

#' Create ggplot layer for contours
#'
#' Create the ggplot contour layers for a donut VPC
#'
#' @param sim_contours output data.frame from the simulate_contours function
#' @param conf_band numeric value indicating the empirical confidence level for 
#' the width of the bands

#' @param colors_bands colors of the confidence bands, string with hex or color
#' 
#' @return  
#' @export
#'
#' @examples
#' 

create_geom_donutVPC <- function(sim_contours, conf_band = 5, colors_bands = c("#99E0DC", "#E498B4"), return_polygons = FALSE) {

  percentiles <- as.numeric(gsub("%", "", unique(sim_contours$percentile), fixed = T))
  pairs_matrix <- unique(sim_contours[, c("var1", "var2")])
  
  colors_bands <- colors_bands[((1:length(percentiles))/2 == round((1:length(percentiles))/2)) + 1]
  names(colors_bands) <- percentiles
  
  contour_geoms <- list()
  conf_geom_data <- list()
  for (p in 1:nrow(pairs_matrix)) {
    var_pair <- paste0(pairs_matrix[p, 1], "-", pairs_matrix[p, 2])
    
    contour_geoms[[var_pair]] <- list()
    for (pr in percentiles) {
      sim_full_df <- sim_contours %>% 
        filter(var1 == pairs_matrix[p, 1], var2 == pairs_matrix[p, 2]) %>% 
        filter(percentile == paste0(pr, "%"))
      kd_sim_full <- ks::kde(sim_full_df[, c("x", "y")], compute.cont = TRUE, approx.cont = FALSE)
      contour_sim_full <- with(kd_sim_full, contourLines(x = eval.points[[1]], y = eval.points[[2]],
                                                         z = estimate, levels = cont[paste0(conf_band, "%")]))
      contour_data <- extract_contour_df(contour_sim_full, kd_sim_full$cont, pr, pairs_matrix[p, ])
      
      conf_geom_data[[paste0(var_pair, ": ", pr, "%")]] <- contour_data %>% 
        create_polygon()
      
      contour_geoms[[var_pair]][[paste0(pr, "%")]] <- geom_sf(data = conf_geom_data[[paste0(var_pair, ": ", pr, "%")]], 
                                                              color = colors_bands[as.character(pr)], fill = colors_bands[as.character(pr)])
    }
    
  }
  
  if (return_polygons) {
    return(conf_geom_data)
  } else {
    return(contour_geoms)
  }
}


#' Create ggplot(s) with geom_vpc
#' 
#' 
#' @param save_path string with path to save the resulting simulations as R 
#' object
#' 
#'
#'
#'
#work in progress
ggVPC_donut <- function(geom_vpc, obs_data, pairs_data = NULL) {

  if (is.null(pairs_data)) {
    pairs_data <- matrix(c(gsub("-([A-Z,a-z])\\w+", "", names(geom_vpc)),
                           gsub("([A-Z,a-z])\\w+-", "", names(geom_vpc))), ncol = 2)
  }
  
  percentiles <- names(geom_vpc[[1]])
  
  if (all(pairs_data %in% colnames(obs_data))) {
    obs_contours <- NULL
    for (p in 1:nrow(pairs_data)) {
      kd_obs <- ks::kde(obs_data %>% select(pairs_data[p, ]), compute.cont = TRUE)
      contour_obs <- with(kd_obs, contourLines(x = eval.points[[1]], y = eval.points[[2]],
                                               z = estimate, levels = cont[percentiles]))
      obs_contours <- rbind.data.frame(obs_contours, extract_contour_df(contour_obs, kd_obs$cont, 0, pairs_data[p, ]))
    }
  } else if (all(colnames(obs_data) %in% 
                 c("level", "x", "y", "circle", "percentile", 
                   "sim_nr", "var1", "var2"))) {
    obs_contours <- obs_data
  } else {
    stop("The obs_data does not seem to be a observed contour, or does not ",
         "contain the same names as the covariate names in pairs data or ", 
         "the geom_vpc object.")
  }
  
  plot_list <- list()
  for (i in 1:nrow(pairs_data)) {
    var_pair <- paste0(pairs_data[i, 1], "-", pairs_data[i, 2])
    plot_list[[var_pair]] <- obs_contours %>% 
      mutate(key = paste(percentile, var1, var2, sim_nr, circle)) %>% 
      filter(var1 == pairs_data[i, 1], var2 == pairs_data[i, 2]) %>% 
      ggplot() +
      geom_vpc[[i]] +
      geom_path(aes(x = x, y = y, color = percentile, group = key), color = "black") +
      labs(x = pairs_data[i, 1], y = pairs_data[i, 2]) +
      theme_bw() + theme(aspect.ratio = 1)
  }
  
  attr(plot_list, "obs_contours") <- obs_contours
  return(plot_list)
}

# - add possibility of combining plots (in matrix?) 
# - add possibility to plot without observed data