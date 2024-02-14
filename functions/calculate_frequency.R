# helper functions for the calculation of frequency ----

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 1 - calculate the race-ethnicity   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# new race categories
get_race <- function(data_set, columns = NULL) {
  
  long_format <- as.data.frame(matrix(nrow = 1, ncol = 5))
  names(long_format) <- c("Hispanic", "White","African American","Asian","Other race")
  long_format[1,1] = nrow(data_set[data_set$Race==1,]) 
  long_format[1,2] = nrow(data_set[data_set$Race==2,]) 
  long_format[1,3] = nrow(data_set[data_set$Race==3,]) 
  long_format[1,4] = nrow(data_set[data_set$Race==4,]) 
  long_format[1,5] = nrow(data_set[data_set$Race==5,]) 
  
  long_format <- reshape2::melt(long_format, na.rm = TRUE)
  return(long_format)
}

get_race_multiple <- function(data_set, m, columns = NULL) {
  data_set <- data_set[,c("Race","simulation_nr")]
  
  if (is.null(columns)) {
    columns <- setdiff(colnames(data_set), c("Race","simulation_nr"))
  }
  n_statistics <- 5
  full_results <- as.data.frame(matrix(nrow = m*n_statistics, ncol = 2))
  names(full_results) <- c("Race", "value")
  full_results$simulation_nr <- rep(1:m, each = n_statistics)
  
  
  
  for(i in unique(data_set$simulation_nr)) {
    sim_results <- get_race(data_set[data_set$simulation_nr == i, ], columns = columns)
    full_results[full_results$simulation_nr == i,c(1,2)] <- sim_results
  }
  
  return(full_results)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 2 - calculate the gender   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get_gender <- function(data_set, columns = NULL) {
#   
#   long_format <- as.data.frame(matrix(nrow = 1, ncol = 2))
#   names(long_format) <- c("Male", "Female")
#   long_format[1,1] = nrow(data_set[data_set$Gender=="Male",]) 
#   long_format[1,2] = nrow(data_set[data_set$Gender=="Female",]) 
#   
#   long_format <- reshape2::melt(long_format, na.rm = TRUE)
#   return(long_format)
# }

get_gender <- function(data_set, columns = NULL) {
  
  long_format <- as.data.frame(matrix(nrow = 1, ncol = 2))
  names(long_format) <- c("Male", "Female")
  long_format[1,1] = nrow(data_set[data_set$Gender==1,]) 
  long_format[1,2] = nrow(data_set[data_set$Gender==2,]) 
  
  long_format <- reshape2::melt(long_format, na.rm = TRUE)
  return(long_format)
}

get_gender_multiple <- function(data_set, m, columns = NULL) {
  data_set <- data_set[,c("Gender","simulation_nr")]
  # data_set <- data_set[, colSums(is.na(data_set)) != nrow(data_set)]
  
  if (is.null(columns)) {
    columns <- setdiff(colnames(data_set), c("Gender","simulation_nr"))
  }
  # n_statistics <- length(columns)*7  + 3*choose(length(columns), 2)
  n_statistics <- 2
  full_results <- as.data.frame(matrix(nrow = m*n_statistics, ncol = 2))
  names(full_results) <- c("Gender", "value")
  full_results$simulation_nr <- rep(1:m, each = n_statistics)
  
  
  
  for(i in unique(data_set$simulation_nr)) {
    sim_results <- get_gender(data_set[data_set$simulation_nr == i, ], columns = columns)
    full_results[full_results$simulation_nr == i,c(1,2)] <- sim_results
  }
  
  return(full_results)
}