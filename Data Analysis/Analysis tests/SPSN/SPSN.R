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

######################
##   time novel    ###
######################
data_mem = data %>% filter(Stage == 3)

ANOVA = two_way_aov(data_mem, "env", "treatment", "E.S2...time..s.")

emm = emmeans(ANOVA[[3]], ~ treatment | env)
contrast(emm, method = "pairwise", adjust = "tukey")

######################
##   Bayesian     ###
######################
mod = cmdstan_model(here("Data Analysis", "Analysis tests", "SPSN", "ANOVA_SPSN.stan"))
path_SPSN = here(path, "fit_SPSN_MEM.rds")

fit = two_way_bayes(mod, data_mem, "E.S2...time..s.", path_SPSN, 7)

######################
##   Plots         ###
######################
sum = data %>% group_by(treatment, env, Stage) %>% summarise(Distance = mean(Distance..m.), SEM = sd(Distance..m.)/sqrt(n())) %>% 
  mutate(Stage = ifelse(Stage == 1, "Habituation", ifelse(Stage == 2, "S1", "S2"))) %>% mutate(SEM_l = Distance - SEM, SEM_h = Distance + SEM)

Stage = ggplot(sum)+geom_line(aes(x=Stage, y=Distance, color = env, linetype = treatment, group = interaction(env, treatment)), linewidth = 1.5)+
  xlab("")+ylab("Total path length (m)")+theme_classic(base_size =20)+scale_color_manual(values = c("0" = "salmon2", "1"= "cadetblue3"),labels = c('0' = "PE", "1" = "EE"), name = "Environment")+
  geom_point(aes(x=Stage,y=Distance, color = env), size = 4)+ geom_errorbar(aes(x=Stage, ymin = SEM_l,ymax=SEM_h, color = env), width =0.12, linewidth =0.8, alpha = 0.6)+
  scale_linetype_manual(values = c("0" = "dashed", "1"= "solid"),labels=c("0" = "SAL", "1" = "CVAD"), name = "Treatment")+
  guides(linetype = guide_legend(override.aes = list(linetype = c("dotted", "solid"))))+
  theme(legend.text = element_text(size = 12),legend.title = element_text(size = 13),legend.key.size = unit(1, "lines"), legend.key.width = unit(1, "cm"))   



SPSN_PL = ggplot(data_loco, aes(x = treatment, y = Distance..m., fill = env)) +
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
  labs(x="Treatment", y = "Total path length (m)")+
  theme_classic(base_size = 20)+theme(legend.position = "none")

SPSN1 = ggplot(data_soc, aes(x = treatment, y = S1...time..s., fill = env)) +
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
  labs(x="Treatment", y = "Time Stranger 1 (s)")+
  theme_classic(base_size = 20)+theme(legend.position = "none")


SPSN2 = ggplot(data_mem, aes(x = treatment, y = E.S2...time..s., fill = env)) +
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
  labs(x="Treatment", y = "Time Stranger 2 (s)")+
  theme_classic(base_size = 20)
Stage|SPSN_PL
SPSN1|SPSN2
