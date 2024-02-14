########################################
# Evaluation of the pediction performance of fat mass data
## y.guo@lacdr.leidenuniv.nl - Jan 2024
########################################

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 0 - Protocols   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1. Generate a complete data set.
# 2. Apply a prespecified missing data mechanism to remove some observations (dataset A).
# 3. Use copula to simulate completed data sets (dataset B) with same amount of individuals in dataset A.
# 4. Compare the metrics of known missing data and simulated missing data.  

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 1 - Load libraries   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyr)
library(dplyr)
library(combinat)
library(ggplot2)
library(ggpubr)
library(rvinecopulib)
library(tibble)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 2 - Source files   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function
source("functions/calculate_metric.R") 
source("functions/calculate_metric Fat.R") 

# data
load("clean_data/nhanes_data_12d_5r.Rdata") # observation data, original scale
load("clean_data/sim_Lscale_back.Rdata") # simulation data, original scale

obs_data <- nhanes_data_5r
combi <- permn(1:5) # 120, separate them in 6 files
i = 19
obs_data[, 2] <- ordered(obs_data$Race, levels = combi[[i]])
levels(obs_data$Race)
obs_data[, 1] <- ordered(obs_data$Gender, levels = c(1,2))
levels(obs_data$Gender)
# assign NA for fat mass if age above 60 years old
sim_para_log_back$Fat[sim_para_log_back$Age >= 60] <- NA


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 3 - step 1: generate complete dataset  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
obs_data <- nhanes_data_5r %>% select(-Gender,-Race) # keep 10 continuous variables

#### statistics of obs_data
# report the missing percentage of variables
ag <- sapply(obs_data, function(x) list(ActualN=nrow(obs_data)-sum(is.na(x)),MissingP=100*(sum(is.na(x)))/nrow(obs_data),mean=round(mean(x,na.rm=T),2),
                                    sd=round(sd(x,na.rm=T),2),min=round(min(x,na.rm=T),2),max=round(max(x,na.rm=T),2)))
ag1 <- t(ag) %>%  
  data.frame() %>%
  mutate(meanSD=paste(mean,sd,sep="±"), 
         range=paste(min,max,sep="~")) %>%
  mutate(statistic = paste0(meanSD," ", "[",range,"]"))
ag1 # results: 56% of missing data in fat
# used for table 1 in main

#### generate complete dataset ----
data_cplt <- na.omit(obs_data)   # 11179 individuals, 10 variables

#### statistics of complete dataset 
ag <- sapply(data_cplt, function(x) list(ActualN=nrow(data_cplt)-sum(is.na(x)),MissingP=100*(sum(is.na(x)))/nrow(data_cplt),mean=round(mean(x,na.rm=T),2),
                                         sd=round(sd(x,na.rm=T),2),min=round(min(x,na.rm=T),2),max=round(max(x,na.rm=T),2)))
ag1 <- t(ag) %>%  
  data.frame() %>%
  mutate(meanSD=paste(mean,sd,sep="±"), 
         range=paste(min,max,sep="~")) %>%
  mutate(statistic = paste0(meanSD," ", "[",range,"]"))
ag1 # range of age: [18,59]

# we assume the missing data in fat is only related to age
# calculate the threshold (ending) age for different missing percentage of Fat data

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 4 - step 2: generate incomplete datasets ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### calculate the threshold ----
findEndAge <- as.data.frame(cbind(End_age=c(19:59),MissingP=rep(0,41)))
for (i in 19:59){
  data_Fat = data_cplt$Fat[data_cplt$Age >= i]
  MissingP=100*(length(data_Fat)/nrow(data_cplt))
  findEndAge$MissingP[findEndAge$End_age == i] = MissingP
}

# for the chosen missing percentage. corresponding cutoff age would be
# 30% : 47 year
# 50% : 38 year
# 60% : 34 year

#### Three missing data set ----
data_30 <- data_cplt %>% 
  mutate(Fat = ifelse(Age >= 47, NA, Fat))

data_50 <- data_cplt %>% 
  mutate(Fat = ifelse(Age >= 38, NA, Fat))

data_60 <- data_cplt %>% 
  mutate(Fat = ifelse(Age >= 34, NA, Fat))

data_cplt <- data_cplt %>% 
  mutate(type = "Complete dataset")
data_30 <- data_30 %>% 
  mutate(type ="30% missing dataset")
data_50 <- data_50 %>% 
  mutate(type ="50% missing dataset")
data_60 <- data_60 %>% 
  mutate(type ="60% missing dataset")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 5 - step 3: construct the copula for each dataset ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data_list <- list(data_cplt,data_30,data_50,data_60)
sim <- list()
data_set_names <- c("data_cplt", "data_30", "data_50", "data_60")
n_sim <- nrow(data_cplt)
m <- 100
seed_nr <- 12345
set.seed(seed_nr)

for (i in 1:4) {
vine <- vine(dat = data_list[[i]][,c(1:10)], # may not use this model for comparison
                   margins_controls = list(mult = NULL, xmin = NaN, xmax = NaN, bw = NA, deg = 2),
                   copula_controls = list(family_set = "parametric", 
                                          par_method = "mle",
                                          selcrit = "aic",
                                          keep_data = TRUE,  
                                          cores = 1,
                                          var_types = rep("c", 10)),
                   weights = numeric(),
                   keep_data = TRUE,
                   cores = 1)
set.seed(12345)
sim[[i]] <-  rvine(n_sim*m, vine) %>% as.data.frame() %>% 
  mutate(simulation_nr = rep(1:m, each = n_sim), data_set = data_set_names[[i]]) 

i = i +1
}
save(sim, file = "clean_data/sim_incplt_fat.Rdata")
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 6 - step 4: compare the marginal metric for fat mass data ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("clean_data/sim_incplt_fat.Rdata")
m <- 100
fatdata <- list()
fatdata_sim <- list()      # simulation from incomplete dataset model
fatdata_sim_cplt <- list() # simulation from complete dataset model

# extract the true fat data for different age cutoff values
ageCutoff <- c(47,38,34)
for (i in 1:3) {
  fatdata[[i]] <- data_cplt %>% 
    filter(Age >= ageCutoff[i]) %>% 
    dplyr::select(Fat) %>% 
    as.data.frame()
}

for (i in 1:3) {
  fatdata_sim_cplt[[i]] <- sim[[1]] %>%  
    filter(Age >= ageCutoff[i]) %>% 
    dplyr::select(Fat, simulation_nr) %>% 
    as.data.frame()
}

for (i in 2:4) {
  data_sim <- sim[[i]] %>% 
    filter(Age >=ageCutoff[i-1]) %>% 
    dplyr::select(Fat, simulation_nr) 
  fatdata_sim[[i-1]] <- as.data.frame(data_sim)
}

m = 100
missing_types <- c("30% fat missing", "50% fat missing", "60% fat missing")
plot_dat <- list()
for (i in 1:3) {
  plot_dat[[i]] <- get_statistics_multiple_f(fatdata_sim_cplt[[i]], m = m) %>% mutate(data_set = "complete") %>% 
    bind_rows(get_statistics_multiple_f(fatdata_sim[[i]], m = m) %>% mutate(data_set = "incomplete") ) %>%
    mutate(value = ifelse(abs(value) == Inf, NA, value)) %>% 
    left_join(get_statistics_f(fatdata[[i]]) %>% rename(observed = value)) %>% 
    mutate(rel_value = (value - observed)/observed) %>% 
    mutate(missing_type = missing_types[i])
    
}
plot_dat_joint <- do.call(rbind, plot_dat) 

save(plot_dat_joint, file = "results/imputation_locally.Rdata")

plot_dat_5d <- plot_dat_joint %>% 
  filter(statistic %in% c("Q5.5%","median","Q95.95%","mean","sd")) %>% 
  mutate(statistic = gsub("sd", "standard deviation", statistic, fixed = TRUE)) %>% 
  mutate(covariate = gsub("_", " - ", covariate, fixed = TRUE))

#### numerical analysis of the metrics ----
fat_missing_info <- plot_dat_5d %>% 
  group_by(statistic, missing_type) %>% 
  summarise(median_rel = median(rel_value))
save(fat_missing_info, file = "results/imputation_locally_analysis.Rdata")

#### plot ----
plot_dat_5d$statistic <- factor(plot_dat_5d$statistic,
                                levels = c("Q5.5%","median","Q95.95%","mean","standard deviation"),
                                labels = c("5th percentile","50th percentile","95th percentile","mean","standard deviation"))

Figure_S4 <- plot_dat_5d %>% 
  ggplot(aes(y = rel_value, x = data_set , color = missing_type)) +
  geom_vline(xintercept = seq(0.5, 10, by = 1), color = "grey95") +
  geom_boxplot(fill = "white") +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = 2, color = "grey65") +
  geom_hline(yintercept = 0, linetype = 1, color = "black") +
  scale_color_manual(values = c("#EFAFA8","#E5CC1A","#5499B2"), limits = force) +
  coord_cartesian(ylim = c(-0.4,0.4))+
  labs(x = "Data type", y = "Relative error", color = "Missing% of fat") +
  scale_x_discrete(guide = guide_axis(angle = 90), expand = expansion(mult = c(0.05, 0.05))) +
  facet_grid(missing_type ~ statistic) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text = element_text(size=6),
        legend.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.4, 'cm'))
Figure_S4
ggsave(Figure_S4, file = "figure/figure_S4.pdf",width = 6.7, height = 4, units = "in")
ggsave(Figure_S4, file = "figure/figure_S4.tiff",width = 6.7, height = 4, units = "in")

