######################################################################
# Development of copula with unordered-categorical covariate variable  ----
## & copula based simulations
## y.guo@lacdr.leidenuniv.nl - Jan 2024
######################################################################

##~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 0 - Load libraries   ----
##~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyr)
library(combinat)
library(rvinecopulib)
library(stringr) 
library(dplyr)

##~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 1 - Source files   ----
##~~~~~~~~~~~~~~~~~~~~~~~~
# data
load("clean_data/nhanes_data_12d_log.Rdata")

##~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 2 - model development   ----
##~~~~~~~~~~~~~~~~~~~~~~~~

# all the possible order combinations of "race-ethnicity"
combi <- permn(1:5) # 120

# create a form to store the information of each run
AICvalue <- rep(NA,length(combi))
Num <- 1:length(combi)# the order of combination
RunningTime <- rep(NA,length(combi))
AICform <- cbind(Num, AICvalue,RunningTime) %>% as.data.frame() 

# set the order of gender
data <- nhanes_data_log
data[, 1] <- ordered(data$Gender, levels = c(1,2))
levels(data$Gender)

# run in loop
for (i in 1 : length(combi)) {
  start_time <- Sys.time()
  data$Race <- ordered(data$Race, levels = combi[[i]])
  vineTempo <- vine(dat = data,
                         margins_controls = list(mult = 1, xmin =NaN, xmax = NaN, bw = NA, deg = 2),
                         copula_controls = list(family_set = "parametric", 
                                                par_method = "mle",
                                                selcrit = "aic",
                                                keep_data = TRUE,  
                                                cores = 1,
                                                var_types = c("d","d",rep("c", 10))),
                         weights = numeric(),
                         keep_data = TRUE,
                         cores = 1) 
  
  string_vine <- capture.output(vineTempo$copula) # two element, AIC in the second element
  getAIC <- word(string_vine[2],start=18,end=18,sep=fixed(" ")) %>% as.numeric()
  end_time <- Sys.time()
  AICform[i,2] <- getAIC
  AICform[i,3] <- end_time - start_time
  cat("time of copula fitting: ", end_time - start_time, "\n",i, "run", "\n") # ~10 min for each run
  save(AICform, file = "results/AICform.Rdata")
}

##~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 3 - model selection : AIC   ----
##~~~~~~~~~~~~~~~~~~~~~~~~
# find the lowest AIC
load("results/AICform.Rdata")
min(AICform$AICvalue) # -112347.4
AICform[AICform$AICvalue == -112347.4,] # best order: No.19

##~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 4 - model rerun   ----
##~~~~~~~~~~~~~~~~~~~~~~~~
# best model with 19th rank combination for race
i = 19
nhanes_data_log$Race <- ordered(nhanes_data_log$Race, levels = combi[[i]])
levels(nhanes_data_log$Race)
nhanes_data_log$Gender <- ordered(nhanes_data_log$Gender, levels = c(1,2))
levels(nhanes_data_log$Gender)
nhanes_copula_log <- vine(dat = nhanes_data_log,
                          margins_controls = list(mult = 1, xmin =NaN, xmax = NaN, bw = NA, deg = 2),
                          copula_controls = list(family_set = "parametric", 
                                                 par_method = "mle",
                                                 selcrit = "aic",
                                                 keep_data = TRUE,  
                                                 cores = 1,
                                                 var_types = c("d","d",rep("c", 10))),
                          weights = numeric(),
                          keep_data = TRUE,
                          cores = 1) 

save(nhanes_copula_log, file = "results/nhanes_copula_Lscale.Rdata")

##~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 5 - simulation   ----
##~~~~~~~~~~~~~~~~~~~~~~~~
load("results/nhanes_copula_Lscale.Rdata")
n_sim <- nrow(nhanes_data_log)
m <- 100
seed_nr <- 12345
set.seed(seed_nr)
sim_para_log <- rvine(n_sim*m, nhanes_copula_log) %>%
  mutate(simulation_nr = rep(1:m, each = n_sim)) 
# 2700800 obs, 13 var
save(sim_para_log, file = "clean_data/sim_Lscale.Rdata")

# data back-transform
sim_para_log_back <- sim_para_log %>% 
  mutate(Fat = exp(logFat)) %>%
  mutate(SCR = exp(logSCR)) %>%
  mutate(ALT = exp(logALT)) %>%
  mutate(AST = exp(logAST)) %>%
  mutate(ALP = exp(logALP)) %>%
  mutate(Albumin = exp(logAlbumin)) %>% 
  mutate(BR = exp(logBR) - 0.01) %>% 
  select(-c(logFat,logSCR,logALT,logAST,logALP,logAlbumin,logBR)) 
# 2700800 obs, 13 var
save(sim_para_log_back, file = "clean_data/sim_Lscale_back.Rdata")