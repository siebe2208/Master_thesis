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
pacman::p_load(tidyverse, ggplot2, here, emmeans, logspline, cmdstanr, posterior, rstanarm, HDInterval, patchwork, effectsize)
path = here("Data Analysis", "Analysis tests", "SPSN", "fits")
source(here("Data Analysis", "Utility.R"))

########################################
## Load in data and mutate conditions ##
########################################
data = read.csv(here("Data", "069-2022-2 SPSN data.csv"))
meta_data = read.csv(here("Data", "data.csv"))

data = data %>% filter(Animal %in% meta_data$mice_ID | Animal == "M-15.2_R") %>% rowwise() %>% 
  mutate(treatment = factor(ifelse(str_detect(Treatment, "SAL"), 0, 1))) %>% 
  mutate(env = factor(ifelse(str_detect(Treatment, "PE"), 0,1))) %>% mutate(Stage = ifelse(Stage == 'HABITUATION', 1, ifelse(Stage == "S1", 2, 3)))


sum(as.numeric(data$treatment)-1) # Check how many animals CVAD (should be 12)
sum(as.numeric(data$env)-1)  # Check how many animals CVAD (should be 18)

data = data %>% filter(Animal != "F-13.2_B") # exclude the mouse that died 

######################
##   effect stage  ###
######################

stage = aov(Distance..m. ~ treatment * env * Stage + Error(Animal/Stage), data = data)
sum = summary(stage)

eta_squared(stage, partial = TRUE)

######################
##   path length   ###
######################
data_loco = data %>% filter(Stage == 1)

ANOVA = two_way_aov(data_loco, "env", "treatment", "Distance..m.")

emm = emmeans(ANOVA[[3]], ~ treatment | env)
contrast(emm, method = "pairwise", adjust = "tukey")

######################
##   Bayesian     ###
######################
mod = cmdstan_model(here("Data Analysis", "Analysis tests", "SPSN", "ANOVA_SPSN.stan"))
path_SPSN = here(path, "fit_SPSN_PL.rds")

fit = two_way_bayes(mod, data_loco, "Distance..m.", path_SPSN, 7)

######################
##   path length   ###
######################
data_soc = data %>% filter(Stage == 2)

ANOVA = two_way_aov(data_soc, "env", "treatment", "S1...time..s.")

emm = emmeans(ANOVA[[3]], ~ env | treatment)
contrast(emm, method = "pairwise", adjust = "tukey")

######################
##   Bayesian     ###
######################
mod = cmdstan_model(here("Data Analysis", "Analysis tests", "SPSN", "ANOVA_SPSN.stan"))
path_SPSN = here(path, "fit_SPSN_SOC.rds")

fit = two_way_bayes(mod, data_soc, "S1...time..s.", path_SPSN, 7)
