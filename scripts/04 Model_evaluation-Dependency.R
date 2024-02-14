########################################
# Simulation-based model evaluation - Dependency ----
## y.guo@lacdr.leidenuniv.nl - Jan 2024
########################################

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 0 - Load libraries   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyr)
library(ks)
library(sf)
library(dplyr)
library(tibble)
library(combinat)
library(reshape2)
library(ggplot2)
library(ggpubr)

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

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 1 - Quantitative evaluation: Dependency comparison   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate the correlation between covariates in observed and simulated populations
NHANES_sim_con <- sim_para_log_back[,c(-1,-2)] 
m <- 100
statistics <- get_statistics_multiple(NHANES_sim_con, m = m) %>% 
  mutate(value = ifelse(abs(value) == Inf, NA, value)) %>% 
  left_join(get_statistics(nhanes_data_5r[,-c(1,2)] ) %>% rename(observed = value)) %>% 
  mutate(rel_value = (value - observed)/observed) %>%  # calculate the relative error
  mutate(cov1 = gsub("\\_.*", "", covariate),
         cov2 = gsub("\\w+.\\_", "", covariate))
statistics[, c("covA", "covB")] <- t(apply(as.data.frame(statistics[, c("cov1", "cov2")]), 1, sort))
statistics_cor <- statistics %>% 
  filter(statistic == "correlation") %>% 
  mutate(covariate = gsub("_", " - ", covariate, fixed = TRUE)) %>% 
  select(-c(Var1,Var2))

# Calculate the overlap metric of covariate pairs of observed and simulated populations
# Overlap metric: overlap of the density contours in observed and simulated data. 
# For each covariate pair, 95th percentile density contours were calculated for observed and simulated populations. 
# The overlap metric was computed as the Jaccard index: the ratio between the intersection area and union area 
# [Calculation takes a long time!!!]
data_obs <- nhanes_data_5r[,-c(1,2)]
overlap <- get_overlap_multiple(data_obs, NHANES_sim_con, m = m) %>% 
  mutate(value = ifelse(abs(value) == Inf, NA, value)) %>% 
  mutate(observed = ifelse(statistic=="overlap",100,observed)) %>% 
  mutate(rel_value = (value - observed)/observed) %>% 
  mutate(cov1 = gsub("\\_.*", "", covariate),
         cov2 = gsub("\\w+.\\_", "", covariate))
overlap[, c("covA", "covB")] <- t(apply(as.data.frame(overlap[, c("cov1", "cov2")]), 1, sort))
save(statistics_cor, file ="results/dependency_metrics_cor.Rdata")
save(overlap, file ="results/dependency_metrics_olp.Rdata")


### numerical analysis ----
# see if there is any similarity between the bad fit if we look at the correlation and overlap metric
corr_info_covariate <- statistics_cor %>% 
  group_by(covariate) %>% 
  summarise(median_err = round(median(value-observed),3)) %>% 
  mutate(abs = abs(median_err))
median(corr_info_covariate$median_err) # median error
corr_info_covariate[order(corr_info_covariate$abs, decreasing = TRUE),] # find the top2 largest error 

overlap_info_covariate <- overlap %>% 
  group_by(covariate) %>% 
  summarise(median = round(median(value),3))
median(overlap_info_covariate$median)

amount_a85 <- nrow(overlap_info_covariate[overlap_info_covariate$median >= 85,])
amount_all <- nrow(overlap_info_covariate)
proportion <- amount_a85/amount_all
overlap_info_covariate[order(overlap_info_covariate$median, decreasing = FALSE),] # find the top2 largest error 


# plot ----
# visualize the correlation and overlap measure of each covariate pair in observed and simulated populations
load("results/dependency_metrics_cor.Rdata")
p_combo_cor <- ggplot(statistics_cor) +
  geom_vline(xintercept = seq(0.5, 46, by = 1), color = "grey95") + 
  geom_hline(aes(yintercept = 0), color = "#515151",linetype = "dashed") +
  geom_boxplot(aes(y = value, x = covariate), fill = "white", outlier.shape = NA) +
  geom_point(aes(y = observed, x = covariate, color = statistic),size = 1.7, shape = 18, alpha = 0.3) + # 2.2 size for manuscript; 1.7 for poster
  labs(y = "Correlation") +
  scale_x_discrete() +
  scale_color_manual(values = c("#EE9964","white"), limits = force) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),legend.position="none",
        strip.text.y = element_text(size=10),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())
p_combo_cor
load("results/dependency_metrics_olp.Rdata")
p_combo_olp <- ggplot(overlap) +
  geom_vline(xintercept = seq(0.5, 46, by = 1), color = "grey95") + 
  geom_hline(aes(yintercept = 100), color = "#515151",linetype = "dashed") +
  geom_hline(aes(yintercept = 85), color = "#515151",linetype = "dashed") +
  geom_boxplot(aes(y = value, x = covariate), fill = "white", outlier.shape = NA) +
  labs(x = "Covariate combinations", y = "Ovarlap (%)") +
  coord_cartesian(ylim = c(0,100)) +
  scale_x_discrete(guide = guide_axis(angle = 90), expand = expansion(mult = c(0.05, 0.05))) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),legend.position="none",
        strip.text.y = element_text(size=10))
p_combo_olp
# combine plots 
Figure_3 <- ggarrange(p_combo_cor, p_combo_olp, ncol = 1, nrow = 2,heights=c(3,5),labels = c("A","B"))
Figure_3

ggsave(Figure_3, file = "figure/figure_3.pdf",width = 6.7, height = 4.5, units = "in") # for manuscript
ggsave(Figure_3, file = "figure/figure_3.tiff",width = 6.7, height = 4.5, units = "in", dpi = 600) # for manuscript
