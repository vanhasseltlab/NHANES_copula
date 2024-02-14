########################################
# Subgroup analysis - Race-ethnicity
## y.guo@lacdr.leidenuniv.nl - Jan 2024
########################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 0 - Load libraries   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## load the packages
library(tidyr)
library(dplyr)
library(combinat)
library(reshape2)
library(rvinecopulib)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 1 - Source files   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function
source("functions/calculate_metric.R") 
source("functions/calculate_frequency.R") 

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

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 2 - Construct subgroup copulas for different race-ethnic groups   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# dataset for modelling: obs_data
# info collection
unique(obs_data$Race)
sim_race_collection <-list()
name_race <- c("Hispanic", "White","African American","Asian","Other race")

for (i in 1:5) {
  start_time <- Sys.time()
  rv <- vine(dat = obs_data[obs_data$Race == i,] %>% select(-Race),
             margins_controls = list(mult = 1, xmin =NaN, xmax = NaN, bw = NA, deg = 2), # set the xmin to avoid negative values
             copula_controls = list(family_set = "parametric", 
                                    par_method = "mle",
                                    selcrit = "aic",
                                    keep_data = TRUE,  
                                    cores = 1,
                                    var_types = c("d",rep("c", 10))),
             weights = numeric(),
             keep_data = TRUE,
             cores = 1)
  end_time <- Sys.time()
  cat("time of copula fitting: ", end_time - start_time, "\n") 
  
  n_sim <- nrow(obs_data[obs_data$Race == i,])
  m <- 100
  set.seed(12345)
  sim <- rvine(n_sim*m, rv) %>%
    mutate(simulation_nr = rep(1:m, each = n_sim))
  
  sim <- sim %>%
    mutate(method = name_race[i])
  
  sim_race_collection[[i]] <- sim
  
}

save(sim_race_collection, file = "clean_data/subgroup_race_simulation.Rdata")
load("clean_data/subgroup_race_simulation.Rdata")

# data preparation
for (i in 1:5) {
  test <- sim_race_collection[[i]] 
  test$Fat[test$Age >= 60] <- NA
  sim_race_collection[[i]] <- test
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 3 - Marginal performance for full copula and subgroup copulas   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## subgroup copula X statistics ----
m <- 100
statistics_race_sub <- list()
# name_race <- c("Hispanic", "White","African American","Asian","Other race")
# rewrite it in R loop
for (i in 1:5) {
  statistics_race_sub[[i]] <- get_statistics_multiple(sim_race_collection[[i]][,-c(1,13)], m = m) %>% 
    mutate(population= name_race[i]) %>% 
    mutate(value = ifelse(abs(value) == Inf, NA, value)) %>% 
    left_join(get_statistics(obs_data[obs_data$Race == i,-c(1,2)] ) %>% 
                rename(observed = value)) 
}
statistics_race_sub <- do.call(rbind, statistics_race_sub) %>% 
  mutate(rel_value = (value - observed)/observed,
         model = "subgroup copula",
         cov1 = gsub("\\_.*", "", covariate),
         cov2 = gsub("\\w+.\\_", "", covariate))

statistics_race_sub[, c("covA", "covB")] <- t(apply(as.data.frame(statistics_race_sub[, c("cov1", "cov2")]), 1, sort))

## full copula X statistics ----
sim_para <- sim_para_log_back
statistics_race_full <- list()
for (i in 1:5) {
  statistics_race_full[[i]] <- get_statistics_multiple(sim_para[sim_para$Race == i,-c(1,2)], m = m) %>% 
    mutate(population= name_race[i]) %>% 
    mutate(value = ifelse(abs(value) == Inf, NA, value)) %>% 
    left_join(get_statistics(obs_data[obs_data$Race == i,-c(1,2)] ) %>% 
                rename(observed = value)) 
}
statistics_race_full <- do.call(rbind, statistics_race_full) %>% 
  mutate(rel_value = (value - observed)/observed,
         model = "full copula",
         cov1 = gsub("\\_.*", "", covariate),
         cov2 = gsub("\\w+.\\_", "", covariate))
statistics_race_full[, c("covA", "covB")] <- t(apply(as.data.frame(statistics_race_full[, c("cov1", "cov2")]), 1, sort))

statistics_race <- rbind(statistics_race_sub,statistics_race_full)
save(statistics_race, file = "results/subgroup_race_marginal_data.Rdata")
load("results/subgroup_race_marginal_data.Rdata")

# plot for error
statistics_race <- statistics_race %>% 
  mutate(error_value = (value - observed)) 

plot_data_race_5d <- statistics_race %>% 
  filter(statistic %in% c("Q5.5%","median","Q95.95%","mean","sd")) %>% 
  mutate(statistic = gsub("sd", "standard deviation", statistic, fixed = TRUE)) %>% 
  mutate(covariate = gsub("_", " - ", covariate, fixed = TRUE))

### numerical analysis ----
num_als_race <- plot_data_race_5d %>% 
  group_by(population, model, covariate,statistic) %>% 
  summarise(mean= mean(rel_value), sd = sd(rel_value), min = min(rel_value), max= max(rel_value), median = median(rel_value))

data_report_Race <- num_als_race %>% 
  group_by(population,model) %>% 
  summarise(min = min(median), max= max(median))
save(data_report_Race, file = "results/subgroup_race_marginal_analysis.Rdata")


plot_data_race_5d$statistic <- factor(plot_data_race_5d$statistic,
                                              levels = c("Q5.5%","median","Q95.95%","mean","standard deviation"),
                                              labels = c("5th percentile","50th percentile","95th percentile","mean","standard deviation"))
plot_data_race_5d$population <- factor(plot_data_race_5d$population,
                                               levels = c("Hispanic", "White","African American","Asian","Other race"),
                                               labels = c("Hispanic", "White","African American","Asian","Other race"))
save(plot_data_race_5d,file="results/subgroup_race_marginal_plotdata.Rdata")

# visualize the relative error of marginal metrics in observed and simulated populations
marginal_race <- plot_data_race_5d %>% 
  ggplot(aes(y = rel_value, x = covariate, color = model)) +
  geom_vline(xintercept = seq(0.5, 10, by = 1), color = "grey95") +
  geom_boxplot(fill = "white",outlier.shape = NA) +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = 2, color = "grey65") +
  geom_hline(yintercept = 0, linetype = 1, color = "black") +
  scale_color_manual(name = "Model",values = c("#C07A92","#8FE2FF"), limits = force) +
  labs(x = "Covariates", y = "Relative error", color = "Model") +
  scale_x_discrete(guide = guide_axis(angle = 90), expand = expansion(mult = c(0.1, 0.1))) +
  facet_grid(statistic ~ population , scales = "free")+
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text = element_text(size=6),
        axis.text=element_text(size=5.5),
        axis.title=element_text(size=10),
        legend.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.4, 'cm'))
marginal_race
ggsave(marginal_race, file = "figure/figure_S2.pdf",width = 6.7, height = 5.5, units = "in")
ggsave(marginal_race, file = "figure/figure_S2.tiff",width = 6.7, height = 5.5, units = "in")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 4 - Dependent performance for full copula and subgroup copulas   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## subgroup copula X overlap ----
m <- 100
overlap_race_sub <- list()

for (i in 1:5) {
  overlap_race_sub[[i]] <- get_overlap_multiple(obs_data[obs_data$Race == i,c(-1,-2)],
                                                sim_race_collection[[i]][,-c(1,13)],
                                                m = m) %>%
    mutate(population= name_race[i]) %>% 
    mutate(value = ifelse(abs(value) == Inf, NA, value)) %>% 
    left_join(get_statistics(obs_data[obs_data$Race == i,-c(1,2)] ) %>% 
                rename(observed = value)) 
}
overlap_race_sub <- do.call(rbind, overlap_race_sub) %>% 
  mutate(rel_value = (value - observed)/observed,
         observed = ifelse(statistic=="overlap",100,observed),
         model = "subgroup copula",
         cov1 = gsub("\\_.*", "", covariate),
         cov2 = gsub("\\w+.\\_", "", covariate))
overlap_race_sub[, c("covA", "covB")] <- t(apply(as.data.frame(overlap_race_sub[, c("cov1", "cov2")]), 1, sort))

## full copula X overlap ----
overlap_race_full <- list()
for (i in 1:5) {
  overlap_race_full[[i]] <- get_overlap_multiple(obs_data[obs_data$Race == i,c(-1,-2)],
                                                sim_para[sim_para$Race == i,-c(1,2)],
                                                m = m) %>%
    mutate(population= name_race[i]) %>% 
    mutate(value = ifelse(abs(value) == Inf, NA, value)) %>% 
    left_join(get_statistics(obs_data[obs_data$Race == i,-c(1,2)] ) %>% 
                rename(observed = value)) 
}
overlap_race_full <- do.call(rbind, overlap_race_full) %>% 
  mutate(rel_value = (value - observed)/observed,
         observed = ifelse(statistic=="overlap",100,observed),
         model = "full copula",
         cov1 = gsub("\\_.*", "", covariate),
         cov2 = gsub("\\w+.\\_", "", covariate))
overlap_race_full[, c("covA", "covB")] <- t(apply(as.data.frame(overlap_race_full[, c("cov1", "cov2")]), 1, sort))

overlap_race <- rbind(overlap_race_sub,overlap_race_full)
save(overlap_race, file = "results/subgroup_race_overlap.Rdata")

### numerical analysis ----
load("results/subgroup_race_overlap.Rdata")
data_report_Race_ovlp <- overlap_race %>% 
  group_by(population, model) %>% 
  summarise(median = median(value))
save(data_report_Race_ovlp, file = "results/subgroup_race_olp_analysis.Rdata")

### plot ----
load("results/subgroup_race_overlap.Rdata")
plot_overlap_race <- overlap_race %>%
  filter(statistic == "overlap") %>%
  mutate(covariate = gsub("_", " - ", covariate, fixed = TRUE)) 
plot_overlap_race$population <- factor(plot_overlap_race$population,
                                              levels = c("Hispanic", "White","African American","Asian","Other race"),
                                              labels = c("Hispanic", "White","African American","Asian","Other race"))


olp_race_plot <- plot_overlap_race %>% 
  ggplot() +
  geom_vline(xintercept = seq(0.5, 46, by = 1), color = "grey95") + 
  geom_hline(yintercept = 100, linetype = 2, color = "grey65") +
  geom_boxplot(aes(y = value, x = covariate, color = model), fill = "white", outlier.shape = NA) +
  labs(x = "Covariate combinations", y = "Overlap (%)") +
  coord_cartesian(ylim = c(0,100))+
  scale_x_discrete(guide = guide_axis(angle = 90), expand = expansion(mult = c(0.05, 0.05))) +
  scale_color_manual(name = "Model",values = c("#C07A92","#8FE2FF"), limits = force) +
  scale_shape_manual(name = "Observation",values = 23) +
  facet_grid(population ~ ., scales = "free")+
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text = element_text(size=6),
        axis.text=element_text(size=5.5),
        axis.title=element_text(size=10),
        legend.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.5, 'cm'))
olp_race_plot # figure 4A
