#############################################################
## Processing NHANES data and generate weighted copula
## y.guo@lacdr.leidenuniv.nl - April 2024
#############################################################

rm(list=ls())

##~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 0 - Load libraries   ----
##~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
library(survey)
library(ggplot2)
library(ggpubr)  
library(reshape2) 
library(tidyr)  
library(dplyr)
library(combinat)
library(rvinecopulib)
source("functions/calculate_metric.R") 
source("functions/calculate_frequency.R") 

##~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 1 - Data list   ----
##~~~~~~~~~~~~~~~~~~~~~~~~

# 1.SEQN     : id of participants
# 2.RIAGENDR : gender 
# 3.RIDAGEYR : age  (year)
# 4.RIDRETH3; RIDRETH1 : race(with/without Asian) 
# 5.BMXWT    : weight (kg)
# 6.BMXHT    : height (cm)
# 7.DXDTOFAT : total body fat (g)
# 8.LBXSCR   : serum creatinine (mg/dL)
# 9.LBXSATSI : ALT (U/L)
# 10.LBXSASSI: AST (U/L)
# 11.LBXSAPSI: ALP (IU/L)
# 12.LBXSAL  : albumin (g/dL)
# 13.LBXSTB  ：total bilirubin (mg/dL)
# 14.WTMEC2YR: full Sample 2 Year MEC Exam Weight, https://wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx

##~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 2 - Data downloading   ----
##~~~~~~~~~~~~~~~~~~~~~~~~
# 2.1 Demographic (DEMO)
# 2017-2018
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/DEMO_J.XPT", tf <- tempfile(), mode="wb")
DEMO_2017_2018 <- foreign::read.xport(tf)[,c("SEQN","RIAGENDR","RIDAGEYR","RIDRETH1","RIDRETH3", "WTMEC2YR")]
# 2015-2016
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/DEMO_I.XPT", tf <- tempfile(), mode="wb")
DEMO_2015_2016 <- foreign::read.xport(tf)[,c("SEQN","RIAGENDR","RIDAGEYR","RIDRETH1","RIDRETH3", "WTMEC2YR")]
# 2013-2014
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/DEMO_H.XPT", tf <- tempfile(), mode="wb")
DEMO_2013_2014 <- foreign::read.xport(tf)[,c("SEQN","RIAGENDR","RIDAGEYR","RIDRETH1","RIDRETH3", "WTMEC2YR")]
# 2011-2012
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/DEMO_G.XPT", tf <- tempfile(), mode="wb")
DEMO_2011_2012 <- foreign::read.xport(tf)[,c("SEQN","RIAGENDR","RIDAGEYR","RIDRETH1","RIDRETH3", "WTMEC2YR")]
# 2009-2010
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2009-2010/DEMO_F.XPT", tf <- tempfile(), mode="wb")
DEMO_2009_2010 <- foreign::read.xport(tf)[,c("SEQN","RIAGENDR","RIDAGEYR","RIDRETH1", "WTMEC2YR")] # without "RIDRETH3" 

# label the year of release
DEMO_2017_2018$release <- 1
DEMO_2015_2016$release <- 2
DEMO_2013_2014$release <- 3
DEMO_2011_2012$release <- 4
DEMO_2009_2010$release <- 5

# 2.2 Body measures (BMX) 
# 2017-2018
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/BMX_J.XPT", tf <- tempfile(), mode="wb")
BMX_2017_2018 <- foreign::read.xport(tf)[,c("SEQN","BMXWT","BMXHT")]
# 2015-2016
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/BMX_I.XPT", tf <- tempfile(), mode="wb")
BMX_2015_2016 <- foreign::read.xport(tf)[,c("SEQN","BMXWT","BMXHT")]
# 2013-2014
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/BMX_H.XPT", tf <- tempfile(), mode="wb")
BMX_2013_2014 <- foreign::read.xport(tf)[,c("SEQN","BMXWT","BMXHT")]
# 2011-2012
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/BMX_G.XPT", tf <- tempfile(), mode="wb")
BMX_2011_2012 <- foreign::read.xport(tf)[,c("SEQN","BMXWT","BMXHT")]
# 2009-2010
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2009-2010/BMX_F.XPT", tf <- tempfile(), mode="wb")
BMX_2009_2010 <- foreign::read.xport(tf)[,c("SEQN","BMXWT","BMXHT")]

# 2.3 Body absorptiometry (DXX) 
# 2017-2018
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/DXX_J.XPT", tf <- tempfile(), mode="wb")
DXX_2017_2018 <- foreign::read.xport(tf)[,c("SEQN","DXDTOFAT")]
# 2015-2016
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/DXX_I.XPT", tf <- tempfile(), mode="wb")
DXX_2015_2016 <- foreign::read.xport(tf)[,c("SEQN","DXDTOFAT")]
# 2013-2014
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/DXX_H.XPT", tf <- tempfile(), mode="wb")
DXX_2013_2014 <- foreign::read.xport(tf)[,c("SEQN","DXDTOFAT")]
# 2011-2012
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/DXX_G.XPT", tf <- tempfile(), mode="wb")
DXX_2011_2012 <- foreign::read.xport(tf)[,c("SEQN","DXDTOFAT")]
# 2009-2010 Not available

# 2.4 Biochemistry Profile (BIO) 
# 2017-2018
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/BIOPRO_J.XPT", tf <- tempfile(), mode="wb")
BIO_2017_2018 <- foreign::read.xport(tf)[,c("SEQN","LBXSCR","LBXSATSI","LBXSASSI","LBXSAPSI","LBXSAL","LBXSTB")]
# 2015-2016
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/BIOPRO_I.XPT", tf <- tempfile(), mode="wb")
BIO_2015_2016 <- foreign::read.xport(tf)[,c("SEQN","LBXSCR","LBXSATSI","LBXSASSI","LBXSAPSI","LBXSAL","LBXSTB")]
# 2013-2014
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/BIOPRO_H.XPT", tf <- tempfile(), mode="wb")
BIO_2013_2014 <- foreign::read.xport(tf)[,c("SEQN","LBXSCR","LBXSATSI","LBXSASSI","LBXSAPSI","LBXSAL","LBXSTB")]
# 2011-2012
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/BIOPRO_G.XPT", tf <- tempfile(), mode="wb")
BIO_2011_2012 <- foreign::read.xport(tf)[,c("SEQN","LBXSCR","LBXSATSI","LBXSASSI","LBXSAPSI","LBXSAL","LBXSTB")]
# 2009-2010
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2009-2010/BIOPRO_F.XPT", tf <- tempfile(), mode="wb")
BIO_2009_2010 <- foreign::read.xport(tf)[,c("SEQN","LBXSCR","LBXSATSI","LBXSASSI","LBXSAPSI","LBXSAL","LBXSTB")]


##~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 3 - Data configuration and adjustment   ----
##~~~~~~~~~~~~~~~~~~~~~~~~

# Append Files
DEMO <- bind_rows(DEMO_2009_2010, DEMO_2011_2012,DEMO_2013_2014,DEMO_2015_2016,DEMO_2017_2018)
BMX <- bind_rows(BMX_2009_2010, BMX_2011_2012,BMX_2013_2014,BMX_2015_2016,BMX_2017_2018)
DXX <- bind_rows(DXX_2011_2012,DXX_2013_2014,DXX_2015_2016,DXX_2017_2018) # there is no DXX report in 2009-2010 release

# Adjustment of biochemistry biofiles
BIO_2017_2018_2 <- BIO_2017_2018 %>%
  mutate(LBXSCR=1.051*LBXSCR-0.06945,                  # serum creatinine
         LBXSATSI=1.013*LBXSATSI+2.688,                # ALT
         LBXSASSI=1.018*LBXSASSI+ 3.762,               # AST
         LBXSAPSI=10^(1.001*log10(LBXSAPSI)-0.04294),  # ALP
         LBXSAL=1.044*LBXSAL+0.01128) %>%              # Albumin
  mutate(LBXSCR=round(LBXSCR,2),                  # serum creatinine
         LBXSATSI=round(LBXSATSI,0),              # ALT
         LBXSASSI=round(LBXSASSI,0),              # AST
         LBXSAPSI=round(LBXSAPSI,0),              # ALP
         LBXSAL=round(LBXSAL,1))                  # Albumin

BIO <- bind_rows(BIO_2009_2010, BIO_2011_2012,BIO_2013_2014,BIO_2015_2016,BIO_2017_2018_2)


nhanes_data <- left_join(DEMO, BMX, by="SEQN") # +2 covariates ( 8 in total)
nhanes_data <- left_join(nhanes_data, DXX, by="SEQN")  # +1 covariates ( 9 in total)
nhanes_data <- left_join(nhanes_data, BIO, by="SEQN")  # +6 covariates (15 in total)


# set the age group and change the class of variables
ageThresholdL <- 18
ageThresholdH <- 80
colnames(nhanes_data)
nhanes_data <- nhanes_data %>%
  filter(RIDAGEYR >= ageThresholdL & RIDAGEYR < ageThresholdH) 

# handle the different records of race: RIDRETH1,-RIDRETH3 --- 6 races in total
nhanes_data <- nhanes_data %>%
  mutate(RIDRETH = case_when(release < 5  ~ RIDRETH3, 
                             release == 5 & RIDRETH1 <= 4 ~ RIDRETH1,
                             release == 5 & RIDRETH1 > 4 ~ 8,)) %>% # adjust of race records (add non-hispanic Asian population)
  filter(RIDRETH!=8) %>% # remove the "Other Race - Including Multi-Racial" category recorded by RIDRETH1
  select(-RIDRETH1,-RIDRETH3) %>% # only keep RIDRETH to indicate race
  mutate_at(.vars = vars("SEQN", "RIAGENDR","RIDRETH"), factor) %>% # change the class of variables
  mutate(RIAGENDR = factor(RIAGENDR, labels  = c("male", "female"))) %>%
  mutate(RIDRETH = factor(RIDRETH, labels = c("Mexican American", "Other Hispanic","Non-Hispanic White","Non-Hispanic Black","Non-Hispanic Asian","Other Race - Including Multi-Racial")))  %>%
  mutate(DXDTOFAT=DXDTOFAT/1000) # convert the unit of total body fat from "g" to "kg"
str(nhanes_data)

##~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 4 - Reshaping the data  ----
##~~~~~~~~~~~~~~~~~~~~~~~~
data <- nhanes_data %>%
  filter(is.na(BMXHT) + is.na(BMXWT)+ is.na(DXDTOFAT) +
           is.na(LBXSCR) +is.na(LBXSATSI) + is.na(LBXSASSI) +
           is.na(LBXSAPSI) + is.na(LBXSAL) + is.na(LBXSAL) != 9) %>%   # generate the 27008 data file
  dplyr::select(-SEQN,-release) %>%
  rename("Age"="RIDAGEYR",
         "Weight"="BMXWT",
         "Height"="BMXHT",
         "Fat"="DXDTOFAT",
         "SCR"="LBXSCR",
         "ALT"="LBXSATSI",
         "AST"="LBXSASSI",
         "ALP"="LBXSAPSI",
         "Albumin"="LBXSAL",
         "BR"="LBXSTB",
         "Race"= "RIDRETH",
         "Gender" = "RIAGENDR") 

# reduce the races from 6 into 5 categories 
data$Race <- factor(data$Race,
                    levels = c("Mexican American", "Other Hispanic", "Non-Hispanic White","Non-Hispanic Black","Non-Hispanic Asian","Other Race - Including Multi-Racial"),
                    labels = c("Hispanic", "Hispanic", "White","African American","Asian","Other race"))
data$Race <- factor(data$Race,
                    levels = c("Hispanic", "White","African American","Asian","Other race"),
                    labels = c("1", "2", "3","4","5"))
data$Gender <- factor(data$Gender,
                      levels = c("male", "female"),
                      labels = c("1", "2")) 
nhanes_data_5r <- data
nhanes_data_5r <- nhanes_data_5r[,c(1,13,2,4,5,6,7,8,9,10,11,12,3)]
# 27008 obs, 12 var
nhanes_data_5r_sw <- nhanes_data_5r
save(nhanes_data_5r_sw, file = "clean_data/nhanes_data_12d_5r_sw.Rdata") 

## lable for race
# 1 Hispanic
# 2 White
# 3 African American
# 4 Asian
# 5 Other race

nhanes_data_sw_log <- nhanes_data_5r_sw %>% 
  mutate(logFat = log(Fat)) %>%
  mutate(logSCR = log(SCR)) %>%
  mutate(logALT = log(ALT)) %>%
  mutate(logAST = log(AST)) %>%
  mutate(logALP = log(ALP)) %>%
  mutate(logAlbumin = log(Albumin)) %>% 
  mutate(logBR = log(BR+0.01)) %>% 
  select(-c(Fat,SCR,ALT,AST,ALP,Albumin,BR))
# 27008 obs, 12 var
save(nhanes_data_sw_log, file = "clean_data/nhanes_data_12d_log_sw.Rdata")

##~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 5 - Clean the environment  ----
##~~~~~~~~~~~~~~~~~~~~~~~~
rm(DEMO_2009_2010, DEMO_2011_2012,DEMO_2013_2014,DEMO_2015_2016,DEMO_2017_2018,
   BMX_2009_2010, BMX_2011_2012,BMX_2013_2014,BMX_2015_2016,BMX_2017_2018,
   DXX_2011_2012,DXX_2013_2014,DXX_2015_2016,DXX_2017_2018,
   BIO_2009_2010, BIO_2011_2012,BIO_2013_2014,BIO_2015_2016,BIO_2017_2018,BIO_2017_2018_2,
   DEMO, BMX, DXX, BIO,
   nhanes_data, data, nhanes_data_5r)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 6 - Modelling copula with weight  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("clean_data/nhanes_data_12d_log_sw.Rdata")
# all the possible order combinations of "race-ethnicity"
combi <- permn(1:5) # 120
i = 19
nhanes_data_sw_log$Race <- ordered(nhanes_data_sw_log$Race, levels = combi[[i]])
levels(nhanes_data_sw_log$Race)
nhanes_data_sw_log$Gender <- ordered(nhanes_data_sw_log$Gender, levels = c(1,2))
levels(nhanes_data_sw_log$Gender)

obs_data <- nhanes_data_sw_log %>% 
  select(-WTMEC2YR)
sample_w <- nhanes_data_sw_log$WTMEC2YR/5 # five releases
start_time <- Sys.time()
nhanes_copula_sw <- vine(dat = obs_data,
                         margins_controls = list(mult = 1, xmin =NaN, xmax = NaN, bw = NA, deg = 2),
                         copula_controls = list(family_set = "parametric", 
                                                par_method = "mle",
                                                selcrit = "aic",
                                                keep_data = TRUE,  
                                                cores = 1,
                                                var_types = c("d","d",rep("c", 10))),
                         weights = sample_w,
                         keep_data = TRUE,
                         cores = 1) 
end_time <- Sys.time()
cat("time of copula fitting: ", end_time - start_time, "weighted") # 30.0533 weighted
save(nhanes_copula_sw, file = "results/nhanes_copula_Lscale_weighted.Rdata")

##~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 7 - simulation   ----
##~~~~~~~~~~~~~~~~~~~~~~~~
## weighted simulation
load("results/nhanes_copula_Lscale_weighted.Rdata")
n_sim <- nrow(nhanes_data_sw_log)
m <- 1
seed_nr <- 12345
set.seed(seed_nr)
sim_para_sw <- rvine(n_sim*m, nhanes_copula_sw)
# 2700800 obs, 13 var
save(sim_para_sw, file = "clean_data/sim_Lscale_sw.Rdata")

# data back-transform
sim_para_sw_back <- sim_para_sw %>% 
  mutate(Fat = exp(logFat)) %>%
  mutate(SCR = exp(logSCR)) %>%
  mutate(ALT = exp(logALT)) %>%
  mutate(AST = exp(logAST)) %>%
  mutate(ALP = exp(logALP)) %>%
  mutate(Albumin = exp(logAlbumin)) %>% 
  mutate(BR = exp(logBR) - 0.01) %>% 
  select(-c(logFat,logSCR,logALT,logAST,logALP,logAlbumin,logBR)) 
# 2700800 obs, 13 var
save(sim_para_sw_back, file = "clean_data/sim_Lscale_back_weighted.Rdata")

## original simulation
load("clean_data/sim_Lscale_back.Rdata") # # simulation data, original scale
sim_para_back <- sim_para_log_back %>% 
  filter(simulation_nr == 1) %>% 
  select(-simulation_nr)
rm(sim_para_log_back)

##~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 8 - Comparision   ----
##~~~~~~~~~~~~~~~~~~~~~~~~
load("clean_data/sim_Lscale_back_weighted.Rdata")
### 8.1 comparision of marginal distribution ----

#### 8.1.1 discrete var ----
##### 8.1.1.1 race ----
target <- c("Hispanic", "White","African American","Asian","Other race")

freq_race_sw <- get_race(sim_para_sw_back) %>%
  rename(Race = variable) %>% 
  dplyr::rename(mean = value) %>% 
  mutate(sd = 0) %>% 
  mutate(type ="weighted VP")

freq_race_uw <- get_race(sim_para_back) %>%
  rename(Race = variable) %>% 
  dplyr::rename(mean = value) %>% 
  mutate(sd = 0) %>% 
  mutate(type ="unweighted VP")

freq_race_combo <- rbind(freq_race_sw,freq_race_uw)

freq_race_combo$Race <- factor(freq_race_combo$Race,
                               levels = target,
                               labels = target)
type_order <- c("weighted VP","unweighted VP")
freq_race_combo$type <- factor(freq_race_combo$type,
                               levels = type_order,
                               labels = type_order)
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
  coord_cartesian(ylim = c(0,20000))+
  theme(strip.background = element_rect(fill = "white"), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size=8))
race_frequency

# comparison to the reported distribution of race origin 2015-2016 
# https://wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx

freq_race_comp <- freq_race_combo %>% 
  mutate(percent = 100*mean/27008)

freq_race_comp$Race <- as.character(freq_race_comp$Race)
freq_race_comp[2,"percent"] = freq_race_comp[2,"percent"] + freq_race_comp[5,"percent"] 
freq_race_comp[2,"Race"] = "White and other" 
freq_race_comp[7,"percent"] = freq_race_comp[7,"percent"] + freq_race_comp[10,"percent"] 
freq_race_comp[7,"Race"] = "White and other" 
freq_race_comp <- freq_race_comp %>% 
  filter(Race != "Other race") %>% 
  mutate(percent = round(x = percent,1))

# visualize the frequency of each race-ethnicity subgroup in observed and simulated populations
race_frequency_comp <-  ggplot(data = freq_race_comp) +
  geom_col(aes(x = Race,y = percent, fill = type, group = factor(type)), position = position_dodge(width = 0.8), width = 0.7, inherit.aes = TRUE, size = 0) +
  geom_text(mapping = aes(label = percent, x = Race, y = percent,group = factor(type)), vjust = "bottom", position = position_dodge(width = 0.8), inherit.aes = TRUE, size = 3) +
  scale_fill_manual(name = NULL, values=c("#FF7F00","gray"),limits = force) +
  labs(x = "Race-ethnicity", y = "Frequency (%)") +
  guides(fill = guide_legend(byrow = TRUE)) +
  guides(color="none") +
  scale_x_discrete(guide = guide_axis(angle = 45), expand = expansion(mult = c(0.15, 0.15))) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size=8))
race_frequency_comp
ggsave(race_frequency_comp, file = "figure/figure_S_race_frequency_comp_USoffical.pdf",width = 5, height = 4, units = "in") 

##### 8.1.1.2 Gender ----
target <- c("Male", "Female")

freq_gender_sw <- get_gender(sim_para_sw_back) %>%
  rename(Gender = variable) %>% 
  dplyr::rename(mean = value) %>% 
  mutate(sd = 0) %>% 
  mutate(type ="weighted VP")

freq_gender_uw <- get_gender(sim_para_back) %>%
  rename(Gender = variable) %>% 
  dplyr::rename(mean = value) %>% 
  mutate(sd = 0) %>% 
  mutate(type ="unweighted VP")
freq_gender_combo <- rbind(freq_gender_sw,freq_gender_uw)
freq_gender_combo$Gender <- factor(freq_gender_combo$Gender,levels = target,labels = target)
type_order <- c("weighted VP","unweighted VP")
freq_gender_combo$type <- factor(freq_gender_combo$type,
                                 levels = type_order,
                                 labels = type_order)

# visualize the frequency of each sex subgroup in observed and simulated populations
gender_frequency <- ggplot() +
  geom_errorbar(data = freq_gender_combo, aes(x = Gender,y = mean,ymin=mean-sd, max=mean+sd, color = type), width=.4, position=position_dodge(width=0.9)) +
  geom_bar(data = freq_gender_combo, aes(x = Gender,y = mean, fill = factor(type)), color = "black", stat="identity", position = "dodge",alpha = .4) +
  scale_fill_manual(name = NULL, values=c("#FF7F00","gray"), breaks= c("weighted VP","unweighted VP"),limits = force) +
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
gender_frequency
Figure_S6 <- ggarrange(race_frequency, gender_frequency, ncol = 2, nrow = 1, widths=c(1.75,1.7))
Figure_S6 

ggsave(Figure_S6, file = "figure/figure_S6.pdf",width = 8, height = 4, units = "in") 

#### 8.1.2 continuous var ----

##### 8.1.2.1 plot for all ----
long_sw <- melt(sim_para_sw_back[,c(-1,-2)]) %>% 
  mutate(data = "weighted VP")

long_uw <- melt(sim_para_back[,c(-1,-2)]) %>% 
  mutate(data = "unweighted VP")

all_line <- rbind(long_sw, long_uw) 
type_order <- c("weighted VP","unweighted VP")
all_line$type <- factor(all_line$data,
                        levels = type_order,
                        labels = type_order)
density_curve_1D <- ggplot(all_line) +
  geom_density(aes(value, color = factor(data, levels = c("weighted VP","unweighted VP")) , linewidth = factor(data)),alpha = 1/10) +
  scale_fill_manual(name = NULL, values=c(rep("grey80",100),"white"),limits = force) +
  scale_colour_manual(name="Data type", values=c("#EE9964", "gray40"), breaks= c("weighted VP","unweighted VP"), limits = force) + 
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
ggsave(density_curve_1D, file = "figure/figure_S6_con.pdf",width = 5, height = 5.7, units = "in") # for manuscript
ggsave(density_curve_1D, file = "figure/igure_S6_con.tiff",width = 6.7, height = 4.5, units = "in", dpi = 600) # for manuscript

##### 8.1.2.2 plot for age ----
long_sw <- cbind(Age = sim_para_sw_back[,"Age"], data = rep("weighted VP",27008)) %>% as.data.frame()
long_uw <- cbind(Age = sim_para_back[,"Age"], data = rep("unweighted VP",27008)) %>% as.data.frame()

all_line <- rbind(long_sw, long_uw) 
all_line$Age <- as.numeric(all_line$Age)
type_order <- c("weighted VP","unweighted VP")
all_line$data <- factor(all_line$data,
                        levels = type_order,
                        labels = type_order)

histogram <- ggplot(all_line, aes(Age, fill = data, color = data)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  labs(x = "Age", y = "Frequency") +
  scale_fill_manual(values = c("weighted VP" = "#EE9964", "unweighted VP" = "gray")) +
  scale_color_manual(values = c("weighted VP" = "#EE9964", "unweighted VP" = "gray10")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0.5, 'cm'),
        legend.direction = "vertical",
        axis.text.x = element_text(size=8),
        axis.title.x = element_text(margin = unit(c(9, 0, 0, 0), "mm"))) 
histogram

Figure_S7 <- ggarrange(Figure_S6, histogram, ncol = 2, nrow = 1, widths=c(1.3, 0.9))
Figure_S7
ggsave(Figure_S7, file = "figure/figure_S_weightedCheck.tiff",width = 11, height = 3.5, units = "in", dpi = 300) # for manuscript


### 8.2 comparison of density contours/tree ----
load("results/nhanes_copula_Lscale_weighted.Rdata") # nhanes_copula_sw
load("results/nhanes_copula_Lscale.Rdata") # nhanes_copula_log

tree_sw <- plot(nhanes_copula_sw$copula, 
                tree = 1, 
                var_names = "use", 
                edge_labels = "family")
tree_sw$labels$title <- "tree_sw"
print(tree_sw)
tree_uw <- plot(nhanes_copula_log$copula, 
                tree = 1, 
                var_names = "use", 
                edge_labels = "family")
tree_uw$labels$title <- "tree_uw"
print(tree_uw)
