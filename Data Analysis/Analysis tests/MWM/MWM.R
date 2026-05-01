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
pacman::p_load(tidyverse, ggplot2, here, emmeans, logspline, cmdstanr, posterior, lme4,lmerTest, HDInterval, patchwork, effectsize, readxl,r2glmm)
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
r2beta(mod, method = "nsj")

data_D10 = data %>% filter(Day == 10)
mod_2 = lmer(Distance ~ env * treatment +(1|Animal), data = data_D10)

summary(mod_2)
r2beta(mod_2, method = "nsj")

sum = data %>% group_by(treatment, env, Day) %>% summarise(Distance = mean(Distance))
ggplot(sum)+geom_line(aes(x=Day, y=Distance, color = env, linetype = treatment, group = interaction(env, treatment)))+theme_classic()
  
###########################
## Pathlength   Bayesian ##
###########################
data_bayes = data %>% ungroup()%>% arrange(Day, Swim) %>% mutate(Animal = rep(1:33, 10*4)) %>% arrange(Animal) %>% 
  mutate(treatment = as.numeric(treatment)-1, env = as.numeric(env)-1)
mod = cmdstan_model(here("Data Analysis", "Analysis tests", "MWM", "ANOVA_MWM.stan"))
path_PL = here(path, "fit_MWM_PLA.rds")

z = as.numeric(scale(data_bayes$Distance))
stan_data = list(N =1320, J = 33, Day  = data_bayes$Day, animal = data_bayes$Animal, 
                 treatment = data_bayes$treatment, environment=data_bayes$env, y = z)

if (file.exists(path_PL)) {fit = readRDS(path_PL)} else {
  fit = mod$sample(data = stan_data,iter_warmup = 1000,iter_sampling = 2000,chains = 4,parallel_chains = 4,
                   output_dir = dirname(path_PL))
  saveRDS(fit, path_PL)}

fit$summary()

prior = as.numeric(fit$draws("prior"))

BF_env = get_bf(prior, as.numeric(fit$draws("beta_env")))
BF_d_e =  get_bf(prior, as.numeric(fit$draws("beta_d_e")))

mod = cmdstan_model(here("Data Analysis", "Analysis tests", "TF", "ANOVA_TF.stan"))
path_D10 = path_PL = here(path, "fit_MWM_D10.rds")

data_D10 = data_D10 %>% ungroup() %>% arrange(Swim) %>% mutate(Animal = rep(1:33, each = 4))

z = as.numeric(scale(data_D10$Distance))
stan_data = list(N =132, J = 33, animal = data_D10$Animal, 
                 treatment = as.numeric(data_D10$treatment)-1, environment=as.numeric(data_D10$env)-1, y = z)

if (file.exists(path_D10)) {fit = readRDS(path_D10)} else {
  fit = mod$sample(data = stan_data,iter_warmup = 1000,iter_sampling = 2000,chains = 4,parallel_chains = 4,
                   output_dir = dirname(path_D10))
  saveRDS(fit, path_D10)}

fit$summary()

draws = fit$draws(c('beta_env', 'beta_t', 'beta_int')) %>% as_draws_df() %>% 
  select(-.chain, -.iteration, -.draw)
BFs = draws %>% summarise(across(everything(), ~ get_bf(prior, .x)))

######################
## distance to target#
######################
mod = lmer(Dis_TQ ~Day * env * treatment +(1|Animal), data = data)

summary(mod)
r2beta(mod, method = "nsj")

z = as.numeric(scale(data_bayes$Dis_TQ))
stan_data = list(N =1320, J = 33, Day  = data_bayes$Day, animal = data_bayes$Animal, 
                 treatment = data_bayes$treatment, environment=data_bayes$env, y = z)

mod = cmdstan_model(here("Data Analysis", "Analysis tests", "MWM", "ANOVA_MWM.stan"))
fit = mod$sample(data = stan_data,iter_warmup = 1000,iter_sampling = 2000,chains = 4,parallel_chains = 4)

draws = fit$draws(c('beta_day', 'beta_env', 'beta_t', 'beta_d_e', 'beta_d_t', 'beta_e_t', 'beta_int')) %>% as_draws_df() %>% 
  select(-.chain, -.iteration, -.draw)


BFs = draws %>% summarise(across(everything(), ~ get_bf(prior, .x)))


######################
##     Plots        ##
######################
sum = data %>% group_by(treatment, env, Day) %>% summarise(Distance = mean(Distance))
sum_2 = data %>% group_by(treatment, env, Day) %>% summarise(mean = mean(Dis_TQ))

Distance = ggplot(sum)+geom_line(aes(x=Day, y=Distance, color = env, linetype = treatment, group = interaction(env, treatment)), linewidth = 1.1)+
  scale_x_continuous(breaks = sort(unique(sum$Day)))+scale_y_continuous(limits=c(0,2000))+theme_classic(base_size = 20)+
  scale_color_manual(values = c("0" = "salmon2", "1" = "cadetblue3"),labels = c("0" = "PE", "1" = "EE"), name = "Environment")+
  labs(x="Day", y = "Total path length (cm)")+
  scale_linetype_discrete(name = "Treatment",labels = c("SAL", "CVAD"))+
  stat_summary(aes(x = Day, y= Distance, color = env, group = interaction(env, treatment)),
               fun = "mean",
               geom = "point",
               size = 3) + theme(legend.position = "none")+
  annotate("rect", xmin = 9.7, xmax = 10.3,   ymin = 100, ymax = 400, 
    fill = NA,            color = "black",      linewidth = 1.3) +
  annotate("text",x = 10, y = 450, label = "p > 0.05",size = 5,fontface = "bold")

DDT = ggplot(sum_2)+geom_line(aes(x=Day, y=mean, color = env, linetype = treatment, group = interaction(env, treatment)), linewidth = 1.1)+
  scale_x_continuous(breaks = sort(unique(sum$Day)))+theme_classic(base_size = 20)+
  scale_color_manual(values = c("0" = "salmon2", "1" = "cadetblue3"),labels = c("0" = "PE", "1" = "EE"), name = "Environment")+
  labs(x="Day", y = "Distance to target (cm)")+
  scale_linetype_discrete(name = "Treatment",labels = c("SAL", "CVAD"))+
  stat_summary(aes(x = Day, y= mean, color = env, group = interaction(env, treatment)),
               fun = "mean",
               geom = "point",
               size = 3) 

Distance|DDT

