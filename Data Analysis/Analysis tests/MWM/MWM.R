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
pacman::p_load(tidyverse, ggplot2, here, emmeans, logspline, cmdstanr, posterior, lme4,lmerTest, HDInterval, patchwork, effectsize, readxl)
path = here("Data Analysis", "Analysis tests", "MWM", "fits")
source(here("Data Analysis", "Utility.R"))

########################################
## Load in data and mutate conditions ##
########################################
data_w1 = read_excel(here("Data", "Statistics-Sabina MWM week 1.xlsx"))
data_w2 = read_excel(here("Data", "Statistics-Sabina MWM week 2.xlsx"))
data_w1 = data_w1[-c(1:3),]
data_w2 = data_w2[-c(1:3),]

names(data_w1)[2] = "Animal"
names(data_w1)[3] = "Treatment"
names(data_w1)[6] = "Distance"
names(data_w1)[9] = "swim time"
names(data_w1)[17] = "Dis_TQ"

names(data_w2)[2] = "Animal"
names(data_w2)[3] = "Treatment"
names(data_w2)[6] = "Distance"
names(data_w2)[9] = "swim time"
names(data_w2)[17] = "Dis_TQ"


meta_data = read.csv(here("Data", "data.csv"))

data = bind_rows(data_w1, data_w2) %>% filter(Animal %in% meta_data$mice_ID | Animal == "M-15.2_R") %>% rowwise() %>% 
  mutate(treatment = factor(ifelse(str_detect(Treatment, "SAL"), 0, 1))) %>% 
  mutate(env = factor(ifelse(str_detect(Treatment, "PE"), 0,1))) %>% mutate(Distance = as.numeric(Distance)) %>% 
  mutate(Day = as.numeric(as.character(Day))) %>% mutate(Dis_TQ = as.numeric(Dis_TQ))


sum(as.numeric(data$treatment)-1) # Check how many animals CVAD (should be 12)
sum(as.numeric(data$env)-1)  # Check how many animals CVAD (should be 18)

data = data %>% filter(!Animal %in% c("F-13.2_B", "F-13.2_R")) # exclude the mouse that died 

######################
## Pathlength       ##
######################
mod = lmer(Distance ~Day * env * treatment +(1|Animal), data = data)

summary(mod)

data_D10 = data %>% filter(Day == 10)
mod_2 = lmer(Distance ~ env * treatment +(1|Animal), data = data_D10)

summary(mod_2)

sum = data %>% group_by(treatment, env, Day) %>% summarise(Distance = mean(Distance))
ggplot(sum)+geom_line(aes(x=Day, y=Distance, color = env, linetype = treatment, group = interaction(env, treatment)))+theme_classic()
  
)

######################
## distance to target#
######################
mod = lmer(Dis_TQ ~Day * env * treatment +(1|Animal), data = data)

summary(mod)
