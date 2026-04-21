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
path = here("Data Analysis", "Analysis tests", "RR", "fits")
source(here("Data Analysis", "unclean_data.R"))
source(here("Data Analysis", "Utility.R"))

########################################
## Load in data and mutate conditions ##
########################################
data = data.frame(animal = mice_ID, treatment = treatment, env = environment, T1 = RR1, 
                  T2 = RR2, T3 = RR3, T4 = RR4)

data = data %>% rowwise() %>% 
  mutate(treatment = factor(ifelse(treatment == "Sal", 0, 1))) %>% 
  mutate(env = factor(ifelse(env == "PE", 0,1))) %>% 
  mutate(mean = (T1 +T2 +T3+T4)/4)


sum(as.numeric(data$treatment)-1) # Check how many animals CVAD (should be 12)
sum(as.numeric(data$env)-1)  # Check how many animals CVAD (should be 18)

data = data %>% filter(animal != "F13.2_B") # exclude the mouse that died 

######################
## Latency     freq ##
######################
ANOVA = two_way_aov(data, "treatment", "env", "mean")

con = emmeans(ANOVA[[3]], pairwise ~ treatment*env, adjust = "tukey")

############################
## Two‑way ANOVA Bayesian ##
############################
mod = cmdstan_model(here("Data Analysis", "Analysis tests", "RR", "ANOVA_RR.stan"))
path_TF = here(path, "fit_RR.rds")

fit = two_way_bayes(mod, data, "mean", path_TF, 7)