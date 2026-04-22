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
data_M = data.frame(animal = mice_ID, treatment = treatment, env = environment, T1 = RR1, 
                  T2 = RR2, T3 = RR3, T4 = RR4)%>% mutate(mean = (T1+T2+T3+T4)/4) %>% 
                  mutate(treatment = factor(treatment, levels = c("Sal", "CVAD"), labels = c("SAL", "CVAD")), env = factor(env, levels=c("PE", "EE")))


data = data_M %>% rowwise() %>% 
  mutate(treatment = factor(ifelse(treatment == "SAL", 0, 1))) %>% 
  mutate(env = factor(ifelse(env == "PE", 0,1))) %>% pivot_longer(cols = c(T1,T2, T3, T4), names_to = "trial", values_to = "latency") %>% 
  mutate(trial = as.numeric(str_remove(trial, "T")), treatment = factor(treatment),env = factor(env))


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

######################
## hierarchical     ##
######################
ANOVA = ANOVA_multi(data,"env", "treatment", "latency", "animal")

con = emmeans(ANOVA[[3]], pairwise ~ treatment*env, adjust = "tukey")

############################
## Bayesian hierarcical   ##
############################
mod = cmdstan_model(here("Data Analysis", "Analysis tests", "RR", "ANOVA_RR.stan"))
path_RR = here(path, "fit_RR_multi.rds")

fit = bayes_multi(mod, data, "latency", path_RR, 34, 4)


############
## plots#### Combine with RR
############
data_M = data_M %>% filter(animal != "F13.2_B")
RR= ggplot(data_M, aes(x = treatment, y = mean, fill = env)) +
  stat_summary(fun = "mean",
               geom = "bar",
               colour = "black",        
               linewidth = 1,          
               width = 0.6,
               position = position_dodge(width = 0.6)) +
  geom_jitter(aes(color = env), size =2.5, alpha = 0.5, show.legend = F, stroke =1.2, shape =21, color = "black",
              position = position_jitterdodge(jitter.width =0.15, dodge.width =0.6))+
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               linewidth = 1,
               width = 0.3,
               position = position_dodge(width = 0.6)) +
  scale_y_continuous(breaks = seq(0,300,50),expand = expansion(mult = c(0, 0.02))) +
  #scale_x_discrete(labels = c("0" = "SAL", "1"= "CVAD"))+
  scale_fill_manual(values = c("PE" = "salmon2", "EE" = "cadetblue3"), name = "")+
  labs(x="Treatment", y = "Latency (s)")+
  theme_classic(base_size =20)