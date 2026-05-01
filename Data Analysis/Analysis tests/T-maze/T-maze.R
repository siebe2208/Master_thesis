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
names(data)[11] = "ARM_R"
names(data)[17] = "ARM_L"

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


###########################
## Pathlength   Bayesian ##
###########################
data_bayes = data %>% ungroup()%>% arrange(Stage) %>% mutate(Animal = rep(1:34, 2)) %>% arrange(Animal) %>% 
  mutate(treatment = as.numeric(treatment)-1, env = as.numeric(env)-1)

mod = cmdstan_model(here("Data Analysis", "Analysis tests", "T-maze", "ANOVA_T.stan"))
path_PL = here(path, "fit_T_PL.rds")

z = as.numeric(scale(data_bayes$Distance))
stan_data = list(N =68, J = 34, Stage  = data_bayes$Stage, animal = data_bayes$Animal, 
                 treatment = data_bayes$treatment, environment=data_bayes$env, y = z)

if (file.exists(path_PL)) {fit = readRDS(path_PL)} else {
  fit = mod$sample(data = stan_data,iter_warmup = 1000,iter_sampling = 2000,chains = 4,parallel_chains = 4,
                   output_dir = dirname(path_PL))
  saveRDS(fit, path_PL)}

fit$summary()

prior = as.numeric(fit$draws("prior"))

draws = fit$draws(c('beta_env', 'beta_day', 'beta_e_t', 'beta_t')) %>% as_draws_df() %>% 
  select(-.chain, -.iteration, -.draw)
BFs = draws %>% summarise(across(everything(), ~ get_bf(prior, .x)))

######################
##         NAP      ##
######################
blocked_side = data %>% filter(Stage == 1) %>% select(Animal, blocked_stage1 = `Blocked arm`)

data_NAP = data %>% left_join(blocked_side, by = "Animal") %>% filter(Stage == 2) %>% 
  mutate(NAP = ifelse(blocked_stage1 == "Left", ARM_L/(ARM_R+ARM_L), ifelse(blocked_stage1 == "Right", ARM_R/(ARM_R+ARM_L), NA)))

mod = two_way_aov(data_NAP,"env", "treatment", "NAP")

summary(mod)
r2beta(mod, method = "nsj")

######################
##  NAP Bayesian    ##
######################
mod = cmdstan_model(here("Data Analysis", "Analysis tests", "OF", "ANOVA_OF.stan"))
path_NAP = here(path, "fit_T_NAP.rds")

fit = two_way_bayes(mod, data_NAP, "NAP", path_NAP)