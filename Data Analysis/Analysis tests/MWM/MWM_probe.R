
#####################
## Load in packages##
#####################
pacman::p_load(tidyverse, ggplot2, here, emmeans, logspline, cmdstanr, posterior, rstanarm, HDInterval, patchwork, effectsize, readxl)
path = here("Data Analysis", "Analysis tests", "MWM", "fits")
source(here("Data Analysis", "Utility.R"))

########################################
## Load in data and mutate conditions ##
########################################
data = read_excel(here("Data", "Statistics-Sabina Probe.xlsx"))
meta_data = read.csv(here("Data", "data.csv"))

data = data[-c(1:3),]
names(data)[17] = "ACQ_dur"
names(data)[23] = "REV_dur"
names(data)[2] = "Animal"
names(data)[3] = "Treatment"


data_P1 = data %>% filter(Probe == 1) %>% filter(Animal %in% meta_data$mice_ID | Animal == "M-15.2_R") %>% rowwise() %>% rowwise() %>% 
  mutate(treatment = factor(ifelse(str_detect(Treatment, "SAL"), 0, 1))) %>% 
  mutate(env = factor(ifelse(str_detect(Treatment, "PE"), 0,1))) %>% mutate(ACQ_dur = as.numeric(ACQ_dur)) 

data_P2 = data %>% filter(Probe == 2) %>% filter(Animal %in% meta_data$mice_ID | Animal == "M-15.2_R") %>% rowwise() %>% rowwise() %>% 
  mutate(treatment = factor(ifelse(str_detect(Treatment, "SAL"), 0, 1))) %>% 
  mutate(env = factor(ifelse(str_detect(Treatment, "PE"), 0,1))) %>% mutate(ACQ_dur = as.numeric(ACQ_dur)) 

data_P3 = data %>% filter(Probe == 3) %>% filter(Animal %in% meta_data$mice_ID | Animal == "M-15.2_R") %>% rowwise() %>% rowwise() %>% 
  mutate(treatment = factor(ifelse(str_detect(Treatment, "SAL"), 0, 1))) %>% 
  mutate(env = factor(ifelse(str_detect(Treatment, "PE"), 0,1))) %>% mutate(ACQ_dur = as.numeric(ACQ_dur)) 

sum(as.numeric(data_P1$treatment)-1) # Check how many animals CVAD (should be 12)
sum(as.numeric(data_P1$env)-1)  # Check how many animals CVAD (should be 18)


data_P1 = data_P1 %>% filter(!Animal %in% c("F-13.2_B", "F-13.2_R"))  # exclude the mouse that died 
data_P3 = data_P3 %>% filter(!Animal %in% c("F-13.2_B", "F-13.2_R"))

######################
##   Probe  1      ###
######################
ANOVA = two_way_aov(data_P1, "env", "treatment", "ACQ_dur")

emm = emmeans(ANOVA[[3]], ~ env| treatment)
contrast(emm, method = "pairwise", adjust = "tukey")

######################
##   Bayesian     ###
######################
mod = cmdstan_model(here("Data Analysis", "Analysis tests", "MWM", "ANOVA_MWM_2way.stan"))
path_P1 = here(path, "fit_MWM_P1.rds")

fit = two_way_bayes(mod, data_P1, "ACQ_dur", path_P1, 7)

######################
##   Probe  2      ###
######################
ANOVA = two_way_aov(data_P2, "env", "treatment", "ACQ_dur")

emm = emmeans(ANOVA[[3]], ~ env| treatment)
contrast(emm, method = "pairwise", adjust = "tukey")

######################
##   Bayesian     ###
######################
mod = cmdstan_model(here("Data Analysis", "Analysis tests", "MWM", "ANOVA_MWM_2way.stan"))
path_P2 = here(path, "fit_MWM_P2.rds")

fit = two_way_bayes(mod, data_P2, "ACQ_dur", path_P2, 7)

######################
##   Probe  3      ###
######################
ANOVA = two_way_aov(data_P3, "env", "treatment", "REV_dur")

emm = emmeans(ANOVA[[3]], ~ env| treatment)
contrast(emm, method = "pairwise", adjust = "tukey")

######################
##   Bayesian     ###
######################
mod = cmdstan_model(here("Data Analysis", "Analysis tests", "MWM", "ANOVA_MWM_2way.stan"))
path_P3 = here(path, "fit_MWM_P3.rds")

fit = two_way_bayes(mod, data_P3, "ACQ_dur", path_P3, 7)
