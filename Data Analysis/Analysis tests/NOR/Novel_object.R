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
pacman::p_load(tidyverse, ggplot2, here, emmeans, logspline, cmdstanr, posterior, rstanarm, HDInterval, patchwork)
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

############
## plots#### 
############

EXP = ggplot(data, aes(x = treatment, y = total_exploration, fill = env)) +
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
               linewidth = 0.8,
               width = 0.3,
               position = position_dodge(width = 0.6)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_x_discrete(labels = c("0" = "SAL", "1"= "CVAD"))+
  scale_fill_manual(values = c("0" = "salmon2", "1" = "cadetblue3"),labels = c("0" = "PE", "1" = "EE"), name = "")+
  labs(x="Treatment", y = "Exploration time (s)")+
  theme_classic(base_size = 20)+theme(legend.position = "none")

NOP = ggplot(data, aes(x = treatment, y = NOP, fill = env)) +
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
               linewidth = 0.8,
               width = 0.3,
               position = position_dodge(width = 0.6)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_x_discrete(labels = c("0" = "SAL", "1"= "CVAD"))+
  scale_fill_manual(values = c("0" = "salmon2", "1" = "cadetblue3"),labels = c("0" = "PE", "1" = "EE"), name = "")+
  labs(x="Treatment", y = "NAP (%)")+
  theme_classic(base_size = 20)+theme()

EXP|NOP