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
pacman::p_load(tidyverse, ggplot2, here, emmeans, logspline, cmdstanr, posterior, HDInterval, lme4, lmerTest, r2glmm)
path = here("Data Analysis", "Analysis tests", "PA", "fits")
source(here("Data Analysis", "unclean_data.R"))
source(here("Data Analysis", "Utility.R"))

########################################
## Load in data and mutate conditions ##
########################################
data_M = data.frame(animal = mice_ID, treatment = treatment, env = environment, T1 = PA1, 
                    T2 = PA2) %>% mutate(mean = (T1+T2)/2) %>% 
  mutate(treatment = factor(treatment, levels = c("Sal", "CVAD"), labels = c("SAL", "CVAD")), env = factor(env, levels=c("PE", "EE")))

data = data_M %>% rowwise() %>% 
  mutate(treatment = factor(ifelse(treatment == "SAL", 0, 1))) %>% 
  mutate(env = factor(ifelse(env == "PE", 0,1))) %>% pivot_longer(cols = c(T1,T2), names_to = "trial", values_to = "latency") %>% 
  mutate(trial = as.numeric(str_remove(trial, "T")), treatment = factor(treatment),env = factor(env)) %>% mutate(trial = trial -1)


sum(as.numeric(data$treatment)-1) # Check how many animals CVAD (should be 12)
sum(as.numeric(data$env)-1)  # Check how many animals CVAD (should be 18)

data = data %>% filter(animal != "F13.2_B") # exclude the mouse that died 

######################
## Pathlength       ##
######################
mod = lmer(latency ~ trial * env * treatment + (1|animal), data = data)

summary(mod)
r2beta(mod, method = "nsj")


###########################
## Pathlength   Bayesian ##
###########################
data_bayes = data %>% ungroup() %>% arrange(trial) %>% mutate(animal = rep(1:34, 2)) %>% arrange(animal) %>% 
  mutate(treatment = as.numeric(treatment)-1, env = as.numeric(env)-1)

mod = cmdstan_model(here("Data Analysis", "Analysis tests", "PA", "ANOVA_PA.stan"))
path_PA = here(path, "fit_PA.rds")

z = as.numeric(scale(data_bayes$latency))
stan_data = list(N =68, J = 34, Day  = data_bayes$trial, animal = data_bayes$animal, 
                 treatment = data_bayes$treatment, environment=data_bayes$env, y = z)

if (file.exists(path_PA)) {fit = readRDS(path_PA)} else {
  fit = mod$sample(data = stan_data,iter_warmup = 1000,iter_sampling = 2000,chains = 4,parallel_chains = 4,
                   output_dir = dirname(path_PA))
  saveRDS(fit, path_PA)}

fit$summary()


prior = as.numeric(fit$draws("prior"))

draws = fit$draws(c('beta_day', 'beta_env', 'beta_t', 'beta_d_e', 'beta_d_t', 'beta_e_t', 'beta_int')) %>% as_draws_df() %>% 
  select(-.chain, -.iteration, -.draw)


BFs = draws %>% summarise(across(everything(), ~ get_bf(prior, .x)))

######################
##   Plots         ###
######################
sum = data %>% group_by(treatment, env, trial) %>% summarise(Distance = mean(latency), SEM = sd(latency)/sqrt(n())) %>% 
  mutate(trial = ifelse(trial == 0, "Shock", "No-Shock")) %>% mutate(SEM_l = Distance - SEM, SEM_h = Distance + SEM) %>% mutate(trial = factor(trial, levels = c("Shock", "No-Shock")))

Stage = ggplot(sum)+geom_line(aes(x=trial, y=Distance, color = env, linetype = treatment, group = interaction(env, treatment)), linewidth = 1.5)+
  xlab("")+ylab("Latency (s)")+theme_classic(base_size =20)+scale_color_manual(values = c("0" = "salmon2", "1"= "cadetblue3"),labels = c('0' = "PE", "1" = "EE"), name = "Environment")+
  geom_point(aes(x=trial,y=Distance, color = env), size = 4)+ geom_errorbar(aes(x=trial, ymin = SEM_l,ymax=SEM_h, color = env), width =0.12, linewidth =0.8, alpha = 0.6)+
  scale_linetype_manual(values = c("1" = "dashed", "0"= "solid"),labels=c("1" = "SAL", "0" = "CVAD"), name = "Treatment")+
  guides(linetype = guide_legend(override.aes = list(linetype = c("dotted", "solid"))))+
  scale_x_discrete(expand = expansion(mult = c(0.15, 0.15)))+
  theme(legend.text = element_text(size = 12),legend.title = element_text(size = 13),legend.key.size = unit(1, "lines"), legend.key.width = unit(1, "cm"))   

