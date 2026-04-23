########################################

# Important !!!!!

# environment 
# 0 = PE
# 1 = EE

# treatment
# 0 = SAL
# 1 = CVAD

########################################

#####################
## Load in packages##
#####################
pacman::p_load(tidyverse, ggplot2, here, emmeans, logspline, cmdstanr, posterior, rstanarm, HDInterval, patchwork, effectsize, readxl)
path = here("Data Analysis", "Analysis tests", "MWM", "fits")
source(here("Data Analysis", "Utility.R"))

########################################
## Load in data and mutate conditions ##
########################################
data = read_excel(here("Data", "Statistics-Sabina MWM week 1.xlsx"))
data = data[-c(1:3),]

names(data)[2] = "Animal"
names(data)[3] = "Treatment"
names(data)[6] = "Distance"
names(data)[9] = "swim time"
names(data)[17] = "Dis_TQ"

meta_data = read.csv(here("Data", "data.csv"))

data = data %>% filter(Animal %in% meta_data$mice_ID | Animal == "M-15.2_R") %>% rowwise() %>% 
  mutate(treatment = factor(ifelse(str_detect(Treatment, "SAL"), 0, 1))) %>% 
  mutate(env = factor(ifelse(str_detect(Treatment, "PE"), 0,1))) %>% mutate(Distance = as.numeric(Distance)) %>% 
  mutate(Day = as.numeric(as.character(Day))) %>% mutate(Dis_TQ = as.numeric(Dis_TQ))


sum(as.numeric(data$treatment)-1) # Check how many animals CVAD (should be 12)
sum(as.numeric(data$env)-1)  # Check how many animals CVAD (should be 18)

data = data %>% filter(Animal != "F-13.2_B") # exclude the mouse that died 

######################
## Pathlength       ##
######################
mod = lmer(Distance ~Day * env *treatment +(1|Animal), data = data)

summary(mod)