########################################

# Important !!!!!

# environment 
# 0 = PE
# 1 = EE

# treatment
# 0 = SAL
# 1 = CVAD

# Stage
# 1 = training
# 2 = test

########################################

#####################
## Load in packages##
#####################
pacman::p_load(tidyverse, ggplot2, here, emmeans, logspline, cmdstanr, posterior, lme4,lmerTest, HDInterval, patchwork, effectsize, readxl,r2glmm)
path = here("Data Analysis", "Analysis tests", "MWM", "fits")
source(here("Data Analysis", "Utility.R"))

########################################
## Load in data and mutate conditions ##
########################################
data= read_excel(here("Data", "Statistics-Sabina T maze.xlsx"))

names(data)[2] = "Animal"
names(data)[3] = "Treatment"
names(data)[5] = "Distance"

meta_data = read.csv(here("Data", "data.csv"))

data = data %>% filter(Animal %in% meta_data$mice_ID | Animal == "M-15.2_R") %>% rowwise() %>% 
  mutate(treatment = factor(ifelse(str_detect(Treatment, "SAL"), 0, 1))) %>% 
  mutate(env = factor(ifelse(str_detect(Treatment, "PE"), 0,1))) %>% mutate(Distance = as.numeric(Distance)) %>% 
  mutate(Stage = ifelse(is.na(`Blocked arm`), 2,1))


sum(as.numeric(data$treatment)-1) # Check how many animals CVAD (should be 12)
sum(as.numeric(data$env)-1)  # Check how many animals CVAD (should be 18)

data = data %>% filter(!Animal %in% c("F-13.2_B")) # exclude the mouse that died 

######################
## Pathlength       ##
######################
mod = lmer(Distance ~ Stage + env * treatment +(1|Animal), data = data)

summary(mod)
r2beta(mod, method = "nsj")
