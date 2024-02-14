#############################################################
## Processing NHANES data
## y.guo@lacdr.leidenuniv.nl - Jan 2024
#############################################################

##~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 0 - Load libraries   ----
##~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
library(survey)
library(ggplot2)
library(GGally)   
library(ggpubr)  
library(reshape2) 
library(tidyr)    

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

##~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 2 - Data downloading   ----
##~~~~~~~~~~~~~~~~~~~~~~~~
# 2.1 Demographic (DEMO)
# 2017-2018
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/DEMO_J.XPT", tf <- tempfile(), mode="wb")
DEMO_2017_2018 <- foreign::read.xport(tf)[,c("SEQN","RIAGENDR","RIDAGEYR","RIDRETH1","RIDRETH3")]
# 2015-2016
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/DEMO_I.XPT", tf <- tempfile(), mode="wb")
DEMO_2015_2016 <- foreign::read.xport(tf)[,c("SEQN","RIAGENDR","RIDAGEYR","RIDRETH1","RIDRETH3")]
# 2013-2014
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/DEMO_H.XPT", tf <- tempfile(), mode="wb")
DEMO_2013_2014 <- foreign::read.xport(tf)[,c("SEQN","RIAGENDR","RIDAGEYR","RIDRETH1","RIDRETH3")]
# 2011-2012
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/DEMO_G.XPT", tf <- tempfile(), mode="wb")
DEMO_2011_2012 <- foreign::read.xport(tf)[,c("SEQN","RIAGENDR","RIDAGEYR","RIDRETH1","RIDRETH3")]
# 2009-2010
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2009-2010/DEMO_F.XPT", tf <- tempfile(), mode="wb")
DEMO_2009_2010 <- foreign::read.xport(tf)[,c("SEQN","RIAGENDR","RIDAGEYR","RIDRETH1")] # without "RIDRETH3" 

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

# save the data
cat("saving nhanes_data.Rdata")
save(nhanes_data, file = "clean_data/nhanes_data.Rdata")

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
nhanes_data_5r <- nhanes_data_5r[,c(1,12,2,3,4,5,6,7,8,9,10,11)]
# 27008 obs, 12 var
save(nhanes_data_5r, file = "clean_data/nhanes_data_12d_5r.Rdata") 

## lable for race
# 1 Hispanic
# 2 White
# 3 African American
# 4 Asian
# 5 Other race

nhanes_data_log <- nhanes_data_5r %>% 
  mutate(logFat = log(Fat)) %>%
  mutate(logSCR = log(SCR)) %>%
  mutate(logALT = log(ALT)) %>%
  mutate(logAST = log(AST)) %>%
  mutate(logALP = log(ALP)) %>%
  mutate(logAlbumin = log(Albumin)) %>% 
  mutate(logBR = log(BR+0.01)) %>% 
  select(-c(Fat,SCR,ALT,AST,ALP,Albumin,BR))
# 27008 obs, 12 var
save(nhanes_data_log, file = "clean_data/nhanes_data_12d_log.Rdata")


##~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ 5 - Clean the environment  ----
##~~~~~~~~~~~~~~~~~~~~~~~~
rm(DEMO_2009_2010, DEMO_2011_2012,DEMO_2013_2014,DEMO_2015_2016,DEMO_2017_2018,
   BMX_2009_2010, BMX_2011_2012,BMX_2013_2014,BMX_2015_2016,BMX_2017_2018,
   DXX_2011_2012,DXX_2013_2014,DXX_2015_2016,DXX_2017_2018,
   BIO_2009_2010, BIO_2011_2012,BIO_2013_2014,BIO_2015_2016,BIO_2017_2018,BIO_2017_2018_2,
   DEMO, BMX, DXX, BIO,
   nhanes_data, data, nhanes_data_5r)
