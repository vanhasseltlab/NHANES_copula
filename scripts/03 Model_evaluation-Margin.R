########################################
# Simulation-based model evaluation -Margin ----
## y.guo@lacdr.leidenuniv.nl - Jan 2024
########################################

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 0 - Load libraries   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyr)
library(dplyr)
library(combinat)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggside)
library(rvinecopulib)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 1 - Source files   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function
source("functions/calculate_metric.R") 
source("functions/calculate_frequency.R") # functions for the batch calculation of frequency of discrete var

# data
load("results/nhanes_copula_Lscale.Rdata")
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

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 2 - Tree structure and density contours ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### plot the first tree ####
name_new <- nhanes_copula_log$copula$names
name_new[1] <- "Sex"
name_new[2] <- "Race-ethnicity"
nhanes_copula_log$copula$names <- name_new
tree_1 <- plot(nhanes_copula_log$copula, 
               tree = 1, 
               var_names = "use", 
               edge_labels = "family")
tree <- tree_1
tree$labels$title <- ""
# visualize the first tree structure of NHANES vine copula with 12 nodes, and 11 edges
print(tree)

### plot the density contours ####
# visualize 6 examples
cat("Visualize set of 6 densities\n")
sim_data <- sim_para_log_back %>% filter(simulation_nr == 1) %>% select(c(-Gender,-Race,-simulation_nr))
obs_data <- obs_data %>%select(c(-Gender,-Race))
all_NHANES <- sim_data %>% mutate(type = "simulated") %>% 
  bind_rows(obs_data %>% mutate(type = "observed")) 

# selected covariate pairs
sets_of_interest <- matrix(c(
  c("SCR", "BR"),
  c("Weight", "Height"),
  c("ALT", "AST"),
  c("Fat", "Albumin"),
  c("Age", "ALT"),
  c("Weight","Fat")),
  ncol = 2, byrow = T)

# plots
p <- list()
for(i in 1 : 6){
  plot_dat_ind <- !is.na(all_NHANES[, sets_of_interest[i, 1]]) & !is.na(all_NHANES[, sets_of_interest[i, 2]])
  data <- all_NHANES[plot_dat_ind, ]
  p[[i]] <- ggplot(data, aes_string(x = sets_of_interest[i, 1], y = sets_of_interest[i, 2], color = "type", linetype = "type")) +
    geom_density_2d(bins = 10,linewidth = 1,show.legend = F) + 
    scale_linetype_manual(values = c(5, 1)) +
    scale_color_manual(values = c("#EE9964", "gray40"), limits = force) +
    scale_x_continuous(expand = expansion(mult = c(0, 0)),
                       limits = quantile(data[,sets_of_interest[i, 1]], probs = c(0.01, 0.95), na.rm = TRUE)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)),
                       limits = quantile(data[,sets_of_interest[i, 2]], probs = c(0.01, 0.95), na.rm = TRUE)) +
    geom_xsidedensity(show.legend = F,linewidth = 1) +
    geom_ysidedensity(show.legend = F,linewidth = 1) +
    scale_ysidex_continuous(minor_breaks = NULL, limits = c(0,NA), expand = expansion(mult = c(0.001, 0.001))) +
    scale_xsidey_continuous(minor_breaks = NULL, limits = c(0,NA), expand = expansion(mult = c(0.001, 0.001))) +
    theme_bw() +
    theme(aspect.ratio = 1, ggside.panel.grid = element_blank(), ggside.axis.line = element_line(color = "white"),
          ggside.axis.text = element_blank(), ggside.axis.ticks = element_blank(), 
          ggside.panel.border = element_rect(colour = "white"), ggside.panel.scale = .2,
          axis.title.x = element_text(color='black',size=12,hjust = 0.4),
          axis.title.y = element_text(color='black',size=12,hjust = 0.4))
  i = i +1
}
p_densityContour <- ggarrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]], ncol = 3, nrow = 2)
p_densityContour 
# combine the tree and density contours

#### ----
library(png)
library(grid)
library(gridExtra)
image <- readPNG("figure/Picture1.png")
raster_grob <- rasterGrob(image, interpolate = TRUE)

Figure_1 <- ggarrange(raster_grob,p_densityContour, ncol = 2, nrow = 1,widths=c(1.5,3),labels=c("A","B"),
                      font.label = list(size = 24))
Figure_1

# revision
ggsave(Figure_1,  file = "results/figure_1.pdf", width = 12, height = 5.5) 
ggsave(Figure_1,  file = "results/figure_1.tiff", width = 12, height = 4.5, dpi = 600) 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 3 - Marginal comparison - discrete var   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### investigate the metrics for discrete variables ####
#### Race-ethnicity ####
m = 100
# Calculate the frequency of each race-ethnicity subgroup
freq_race <- get_race_multiple(sim_para_log_back, m = m) %>% 
  mutate(Race = factor(Race, labels = c("Hispanic", "White","African American","Asian","Other race"))) 
target <- c("Hispanic", "White","African American","Asian","Other race")
freq_race$Race <- factor(freq_race$Race,
                                levels = target,
                                labels = target)
freq_race <- freq_race %>% 
  group_by(Race) %>% 
  summarise(mean= mean(value), sd = sd(value)) %>% 
  mutate(type ="Virtual population")

freq_race_obs <- get_race(nhanes_data_5r) %>%
  rename(Race = variable) %>% 
  dplyr::rename(mean = value) %>% 
  mutate(sd = 0) %>% 
  mutate(type ="Observed population")

freq_race_combo <- rbind(freq_race_obs,freq_race)
# visualize the frequency of each race-ethnicity subgroup in observed and simulated populations
race_frequency <-  ggplot() +
  geom_errorbar(data = freq_race_combo, aes(x = Race,y = mean,ymin=mean-sd, max=mean+sd, color = type), width=.3, position=position_dodge(width=0.9)) + 
  geom_bar(data = freq_race_combo, aes(x = Race,y = mean, fill = type), color = "black", alpha = .4,position = "dodge", stat="identity") +
  scale_fill_manual(name = NULL, values=c("#FF7F00","gray"),limits = force) +
  scale_color_manual(name = NULL, values=c("white","black"),limits = force) +
  labs(x = "Race-ethnicity", y = "Frequency") +
  guides(fill="none") +
  guides(color="none") +
  scale_x_discrete(guide = guide_axis(angle = 45), expand = expansion(mult = c(0.15, 0.15))) +
  theme_bw() +
  coord_cartesian(ylim = c(0,10000))+
  theme(strip.background = element_rect(fill = "white"), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size=8))

#### Gender ####
target <- c("Male", "Female")
# Calculate the frequency of each sex subgroup
freq_gender <- get_gender_multiple(sim_para_log_back, m = m) %>% 
  mutate(Gender = factor(Gender, labels = c("Male", "Female"))) 
freq_gender$Gender <- factor(freq_gender$Gender,levels = target,labels = target)
freq_gender <- freq_gender %>% 
  group_by(Gender) %>% 
  summarise(mean= mean(value), sd = sd(value)) %>% 
  mutate(type ="Virtual population")

freq_gender_obs <- get_gender(nhanes_data_5r) %>%
  rename(Gender = variable) %>% 
  dplyr::rename(mean = value) %>% 
  mutate(sd = 0) %>% 
  mutate(type ="Observed population")

freq_gender_combo <- rbind(freq_gender_obs,freq_gender)
# visualize the frequency of each sex subgroup in observed and simulated populations
gender_frequency <- ggplot() +
  geom_errorbar(data = freq_gender_combo, aes(x = Gender,y = mean,ymin=mean-sd, max=mean+sd, color = type), width=.4, position=position_dodge(width=0.9)) +
  geom_bar(data = freq_gender_combo, aes(x = Gender,y = mean, fill = factor(type)), color = "black", stat="identity", position = "dodge",alpha = .4) +
  scale_fill_manual(name = NULL, values=c("#FF7F00","gray"), breaks= c("Observed population","Virtual population"),limits = force) +
  scale_color_manual(name = NULL, values=c("white","black"),limits = force) +
  labs(x = "Sex", y = "Frequency") +
  guides(fill = guide_legend(byrow = TRUE)) +
  guides(color="none") +
  coord_cartesian(ylim = c(0,18000))+
  scale_x_discrete(guide = guide_axis(angle = 45), expand = expansion(mult = c(0.75,0.75))) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.title = element_blank(),
        legend.spacing.y = unit(0.5, 'cm'),
        legend.direction = "vertical",
        axis.text.x = element_text(size=8),
        axis.title.x = element_text(margin = unit(c(9, 0, 0, 0), "mm"))) 
Figure_2_A <- ggarrange(race_frequency, gender_frequency, ncol = 2, nrow = 1, widths=c(1.75,1.7))
Figure_2_A 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 4 - Marginal comparison - continuous var   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### Density curves ####
long_sim <- melt(sim_para_log_back[,c(-1,-2)], id.vars = "simulation_nr") %>% 
  mutate(data = "Virtual population")

long_obs <- melt(nhanes_data_5r[,c(-1,-2)]) %>% 
  mutate(simulation_nr = 101,
         data = "Observed population")

all_line <- rbind(long_sim, long_obs) 
density_curve_1D <- ggplot(all_line) +
  geom_density(aes(value, color = factor(data, levels = c("Virtual population", "Observed population")) , fill = factor(simulation_nr), linewidth = factor(data)),alpha = 1/10) +
  scale_fill_manual(name = NULL, values=c(rep("grey80",100),"white"),limits = force) +
  scale_colour_manual(name="Data type", values=c("#EE9964", "gray40"), breaks= c("Observed population","Virtual population"), limits = force) + 
  scale_linewidth_manual(name = NULL, values = c(0.3,0.5)) +
  facet_wrap(variable~.,scales="free",nrow = 4) +
  labs(x = "", y = "Density") +
  guides(fill = "none", linewidth = "none") +
  guides(color = guide_legend(byrow = TRUE)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(.55, .10),
        legend.direction = "vertical",
        legend.spacing.y = unit(0.5, 'cm'),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        strip.text.x = element_text(size=6.5),
        axis.text.x = element_text(size=5.5),
        axis.title=element_text(size=10))
density_curve_1D

#### Marginal metrics ####
NHANES_sim_con <- sim_para_log_back[,c(-1,-2)] 
m <- 100
statistics <- get_statistics_multiple(NHANES_sim_con, m = m) %>% 
  mutate(value = ifelse(abs(value) == Inf, NA, value)) %>% 
  left_join(get_statistics(nhanes_data_5r[,-c(1,2)]) %>% rename(observed = value)) %>% 
  mutate(rel_value = (value - observed)/observed) %>%  # calculate the relative error
  mutate(cov1 = gsub("\\_.*", "", covariate),
         cov2 = gsub("\\w+.\\_", "", covariate))
statistics[, c("covA", "covB")] <- t(apply(as.data.frame(statistics[, c("cov1", "cov2")]), 1, sort))

statistics_5d <- statistics %>% 
  filter(statistic %in% c("Q5.5%","median","Q95.95%","mean","sd")) %>% 
  mutate(statistic = gsub("sd", "standard deviation", statistic, fixed = TRUE)) %>% 
  mutate(covariate = gsub("_", " - ", covariate, fixed = TRUE))

statistics_5d$statistic <- factor(statistics_5d$statistic,
                                  levels = c("Q5.5%","median","Q95.95%","mean","standard deviation"),
                                  labels = c("5th P","median","95th P","mean","SD"))
save(statistics_5d, file = "results/overall_marginal.Rdata")
### numerical analysis ----
data_report_statistics_5d <- statistics_5d %>% 
  group_by(statistic, covariate) %>% 
  summarise(avg = mean(value),
            sdvalue = sd(rel_value),
            cv = sd(value)/mean(value))

data_report_statistics_5d_1 <- data_report_statistics_5d %>% 
  filter(statistic %in% c("5th P", "median", "9th P", "mean"))
print(max(data_report_statistics_5d_1$cv))

data_report_statistics_5d_2 <- data_report_statistics_5d %>% 
  filter(!statistic %in% c("5th P", "median", "9th P", "mean"))
print(max(data_report_statistics_5d_2$cv))

### plot ----
# visualize the relative error of marginal metrics in observed and simulated populations
plot_statistics_5d <- statistics_5d %>% 
  ggplot(aes(y = rel_value, x = covariate)) +
  geom_vline(xintercept = seq(0.5, 10, by = 1), color = "grey95") +
  geom_boxplot(fill = "white", color = "black",outlier.shape = NA) +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = 2, color = "grey65") +
  geom_hline(yintercept = 0, linetype = 1, color = "black") +
  labs(x = "Covariates", y = "Relative error", color = "Method") +
  scale_x_discrete(guide = guide_axis(angle = 90), expand = expansion(mult = c(0.1, 0.1))) +
  coord_cartesian(ylim = c(-0.22,0.22))+ 
  facet_grid(statistic~.)+
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text=element_text(size=6.5),
        axis.title=element_text(size=10),
        strip.text.x = element_text(size=7))
plot_statistics_5d

Figure_2_B <- ggarrange(density_curve_1D, plot_statistics_5d, ncol = 2, nrow = 1, widths=c(1.4,1),labels = c("B","C"))
Figure_2 <- ggarrange(Figure_2_A,Figure_2_B, ncol = 1, nrow = 2, heights=c(0.8,1),labels = c("A",""))
ggsave(Figure_2, file = "figure/figure_2.pdf",width = 8, height = 8, units = "in") 
ggsave(Figure_2, file = "figure/figure_2.tiff",width = 8, height = 8, units = "in", dpi = 600) 
