########################################
# Subgroup analysis - Sex
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
##  ~ 2 - Construct subgroup copulas for different sex groups   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# dataset for modelling: obs_data
# info collection
unique(obs_data$Gender)
sim_gender_collection <-list()
name_gender <- c("Male", "Female")

for (i in 1:2) {
  start_time <- Sys.time()
  rv <- vine(dat = obs_data[obs_data$Gender == i,] %>% select(-Gender),
             margins_controls = list(mult = 1, xmin = NaN, xmax = NaN, bw = NA, deg = 2), # set the xmin to avoid negative values
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
  
  n_sim <- nrow(obs_data[obs_data$Gender == i,])
  m <- 100
  set.seed(12345)
  sim <- rvine(n_sim*m, rv) %>%
    mutate(simulation_nr = rep(1:m, each = n_sim))
  
  sim <- sim %>%
    mutate(method = name_gender[i])
  
  sim_gender_collection[[i]] <- sim
  
}

save(sim_gender_collection, file = "clean_data/subgroup_gender_simulation.Rdata")

load("clean_data/subgroup_gender_simulation.Rdata")

# data preparation
for (i in 1:2) {
  test <- sim_gender_collection[[i]] 
  test$Fat[test$Age >= 60] <- NA
  sim_gender_collection[[i]] <- test
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 3 - Marginal performance for full copula and subgroup copulas   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### subgroup copula X statistics ----
m <- 100
statistics_gender_sub <- get_statistics_multiple(sim_gender_collection[[1]][,-c(1,13)], m = m) %>% 
  mutate(population= "Male",
         value = ifelse(abs(value) == Inf, NA, value)) %>% 
  left_join(get_statistics(obs_data[obs_data$Gender == 1,-c(1,2)] ) %>% 
              rename(observed = value)) %>% 
  bind_rows(get_statistics_multiple(sim_gender_collection[[2]][,-c(1,13)], m = m) 
            %>% mutate(population = "Female",
                       value = ifelse(abs(value) == Inf, NA, value)) %>% 
              left_join(get_statistics(obs_data[obs_data$Gender == 2,-c(1,2)] ) %>% 
                          rename(observed = value))) %>% 
  mutate(rel_value = (value - observed)/observed,
         cov1 = gsub("\\_.*", "", covariate),
         cov2 = gsub("\\w+.\\_", "", covariate),
         model = "subgroup copula")
statistics_gender_sub[, c("covA", "covB")] <- t(apply(as.data.frame(statistics_gender_sub[, c("cov1", "cov2")]), 1, sort))

### full copula X statistics ----
sim_para <- sim_para_log_back
statistics_gender_full <- get_statistics_multiple(sim_para[sim_para$Gender == 1,-c(1,2)], m = m) %>% 
  mutate(population= "Male",
         value = ifelse(abs(value) == Inf, NA, value)) %>% 
  left_join(get_statistics(obs_data[obs_data$Gender == 1,-c(1,2)] ) %>% 
              rename(observed = value)) %>% 
  bind_rows(get_statistics_multiple(sim_para[sim_para$Gender == 2,-c(1,2)], m = m) %>% 
              mutate(population = "Female",
                     value = ifelse(abs(value) == Inf, NA, value)) %>% 
              left_join(get_statistics(obs_data[obs_data$Gender == 2,-c(1,2)] ) %>% 
                          rename(observed = value))) %>% 
  mutate(rel_value = (value - observed)/observed,
         cov1 = gsub("\\_.*", "", covariate),
         cov2 = gsub("\\w+.\\_", "", covariate),
         model = "full copula")
statistics_gender_full[, c("covA", "covB")] <- t(apply(as.data.frame(statistics_gender_full[, c("cov1", "cov2")]), 1, sort))

statistics_gender <- rbind(statistics_gender_sub,statistics_gender_full)
save(statistics_gender, file = "results/subgroup_gender_marginal_data.Rdata")
load("results/subgroup_gender_marginal_data.Rdata")

# plot for error
statistics_gender <- statistics_gender %>% 
  mutate(error_value = (value - observed)) 

plot_data_gender_5d <- statistics_gender %>% 
  filter(statistic %in% c("Q5.5%","median","Q95.95%","mean","sd")) %>% 
  mutate(statistic = gsub("sd", "standard deviation", statistic, fixed = TRUE)) %>% 
  mutate(covariate = gsub("_", " - ", covariate, fixed = TRUE))

### numerical analysis ----
num_als_sex <- plot_data_gender_5d %>% 
group_by(population, model, covariate,statistic) %>% 
summarise(median = median(rel_value))

data_report_sex <- num_als_sex %>% 
  group_by(population,model) %>% 
  summarise(min = min(median), max= max(median))
save(data_report_sex, file = "results/subgroup_gender_marginal_analysis.Rdata")

### plot ----
plot_data_gender_5d$statistic <- factor(plot_data_gender_5d$statistic,
                                      levels = c("Q5.5%","median","Q95.95%","mean","standard deviation"),
                                      labels = c("5th percentile","50th percentile","95th percentile","mean","standard deviation"))
plot_data_gender_5d$population <- factor(plot_data_gender_5d$population,
                                         levels = c("Male", "Female"),
                                         labels = c("Male", "Female"))
save(plot_data_gender_5d,file="results/subgroup_gender_marginal_plotdata.Rdata")

marginal_gender <- plot_data_gender_5d %>% 
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
marginal_gender
ggsave(marginal_gender, file = "figure/figure_S3.pdf",width = 4, height = 5, units = "in")
ggsave(marginal_gender, file = "figure/figure_S3.tiff",width = 4, height = 5, units = "in")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 4 - Dependent performance for full copula and subgroup copulas   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### subgroup copula X overlap ----
m <- 100
overlap_gender_sub <- get_overlap_multiple(obs_data[obs_data$Gender == 1,c(-1,-2)],sim_gender_collection[[1]][,-c(1,13)], m = m) %>% 
  mutate(population= "Male",
         value = ifelse(abs(value) == Inf, NA, value)) %>% 
  bind_rows(get_overlap_multiple(obs_data[obs_data$Gender == 2,c(-1,-2)], sim_gender_collection[[2]][,-c(1,13)], m = m) %>% 
              mutate(population = "Female",
                     value = ifelse(abs(value) == Inf, NA, value))) %>% 
  mutate(observed = ifelse(statistic=="overlap",100,observed),
         rel_value = (value - observed)/observed,
         cov1 = gsub("\\_.*", "", covariate),
         cov2 = gsub("\\w+.\\_", "", covariate),
         model = "subgroup copula")
overlap_gender_sub[, c("covA", "covB")] <- t(apply(as.data.frame(overlap_gender_sub[, c("cov1", "cov2")]), 1, sort))

### full copula X overlap ----
overlap_gender_full <- get_overlap_multiple(obs_data[obs_data$Gender == 1,c(-1,-2)],sim_para[sim_para$Gender == 1,-c(1,2)], m = m) %>% 
  mutate(population= "Male",
         value = ifelse(abs(value) == Inf, NA, value)) %>% 
  bind_rows(get_overlap_multiple(obs_data[obs_data$Gender == 2,c(-1,-2)],sim_para[sim_para$Gender == 2,-c(1,2)], m = m) %>% 
              mutate(population = "Female",
                     value = ifelse(abs(value) == Inf, NA, value))) %>% 
  mutate(observed = ifelse(statistic=="overlap",100,observed),
         rel_value = (value - observed)/observed,
         cov1 = gsub("\\_.*", "", covariate),
         cov2 = gsub("\\w+.\\_", "", covariate),
         model = "full copula")
overlap_gender_full[, c("covA", "covB")] <- t(apply(as.data.frame(overlap_gender_full[, c("cov1", "cov2")]), 1, sort))

overlap_gender <- rbind(overlap_gender_sub,overlap_gender_full)
save(overlap_gender, file = "results/subgroup_gender_overlap.Rdata")

### numerical analysis ----
load("results/subgroup_gender_overlap.Rdata")
data_report_gender_ovlp <- overlap_gender %>% 
  group_by(population, model) %>% 
  summarise(median = median(value))
save(data_report_gender_ovlp, file = "results/subgroup_gender_olp_analysis.Rdata")

### plot ----
load("results/subgroup_gender_overlap.Rdata")
plot_overlap_gender <- overlap_gender %>%
  mutate(covariate = gsub("_", " - ", covariate, fixed = TRUE)) 
plot_overlap_gender$population <- factor(plot_overlap_gender$population,
                                         levels = c("Male", "Female"),
                                         labels = c("Male", "Female"))

olp_gender_plot <- plot_overlap_gender %>% 
  ggplot() +
  geom_vline(xintercept = seq(0.5, 46, by = 1), color = "grey95") + # grid background
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
olp_gender_plot # figure 4B

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 5 - combine the gender overlap and race overlap   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Figure_4 <- ggarrange(olp_race_plot,olp_gender_plot, ncol = 1, nrow = 2, heights=c(1,0.5),labels = c("A","B"))
ggsave(Figure_4, file = "figure/figure_4.pdf",width = 6.7, height = 7, units = "in") 
ggsave(Figure_4, file = "figure/figure_4.tiff",width = 6.7, height = 7, units = "in", dpi = 600) 
