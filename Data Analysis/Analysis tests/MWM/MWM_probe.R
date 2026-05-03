
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
  mutate(env = factor(ifelse(str_detect(Treatment, "PE"), 0,1))) %>% mutate(REV_dur = as.numeric(REV_dur)) 

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

############
## plots####
############
P1= ggplot(data_P1, aes(x = treatment, y = ACQ_dur, fill = env)) +
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
  labs(x="Treatment", y = "Duration in TQ (s)")+
  theme_classic(base_size = 20)+theme(legend.position = "none")



P2 = ggplot(data_P2, aes(x = treatment, y = ACQ_dur, fill = env)) +
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
  labs(x="Treatment", y = "Duration in TQ (s)")+
  theme_classic(base_size = 20)

P1|P2



P3 = ggplot(data_P3, aes(x = treatment, y = REV_dur, fill = env)) +
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
  labs(x="Treatment", y = "Duration in RQ (s)")+
  theme_classic(base_size = 20)
