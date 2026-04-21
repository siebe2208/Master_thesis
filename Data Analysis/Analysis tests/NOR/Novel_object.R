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
pacman::p_load(tidyverse, ggplot2, here, emmeans, logspline, cmdstanr, posterior, rstanarm, HDInterval)
path = here("Data Analysis", "Analysis tests", "NOR", "fit")
source(here("Data Analysis", "Utility.R"))

########################################
## Load in data and mutate conditions ##
########################################
data = read.csv(here("Data", "026-2022-2 NOR DATA.csv"))
meta_data = read.csv(here("Data", "data.csv"))

data = data %>% filter(Animal %in% meta_data$mice_ID | Animal == "M-15.2_R") %>% rowwise() %>% 
  mutate(treatment = factor(ifelse(str_detect(Treatment, "SAL"), 0, 1))) %>% 
  mutate(env = factor(ifelse(str_detect(Treatment, "PE"), 0,1))) %>% mutate(total_exploration = Novel.zone...time..s. + Familiar.zone...time..s.) %>% 
  mutate(NOP = Novel.zone...time..s./total_exploration)


sum(as.numeric(data$treatment)-1) # Check how many animals CVAD (should be 12)
sum(as.numeric(data$env)-1)  # Check how many animals CVAD (should be 18)

exclusion_ID = data %>% filter(Stage == "NOR training" & total_exploration < 15) %>% select(Animal)
data = data %>% filter(Animal != "F-13.2_B") %>% filter(!Animal %in% exclusion_ID$Animal) # exclude the mouse that died and when exploration time < 15

data = data %>% filter(Stage == "NOR test")

######################
## path length freq ##
######################
ANOVA = two_way_aov(data, "treatment", "env", "total_exploration")

con = emmeans(ANOVA[[3]], pairwise ~ treatment*env, adjust = "tukey")

############################
## Two‑way ANOVA Bayesian ##
############################
mod = cmdstan_model(here("Data Analysis", "Analysis tests", "NOR", "ANOVA_NOR.stan"))
path_PL = here(path, "fit_NOR_exp.rds")

fit = two_way_bayes(mod, data, "total_exploration", path_PL, 7)

######################
##    NOP freq      ##
######################

ANOVA = two_way_aov(data, "treatment", "env", "NOP")

con_1 = emmeans(ANOVA[[3]], pairwise ~ treatment*env, adjust = "tukey")

#######################################
## Two‑way ANOVA NOP Bayesian        ##
#######################################

mod = cmdstan_model(here("Data Analysis", "Analysis tests", "NOR", "ANOVA_NOR.stan"))
path_PL = here(path, "fit_NOR_NOP.rds")

fit_1 = two_way_bayes(mod, data, "NOP", path_PL, 7)

######################
## central T freq   ##
######################
ANOVA = two_way_aov(data, "treatment", "env", "Centre...time..s.")

con_2 = emmeans(ANOVA[[3]], pairwise ~ treatment*env, adjust = "tukey")

#######################################
## Two‑way ANOVA  Central Bayesian ####
#######################################
mod = cmdstan_model(here("Data Analysis", "Analysis tests", "OF", "ANOVA_OF.stan"))
path_PL = here(path, "fit_OF_CT.rds")

fit_2 = two_way_bayes(mod, data, "Centre...time..s.", path_PL, 7)
