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
path = here("Data Analysis", "Analysis tests", "BFT", "fits")
source(here("Data Analysis", "unclean_data.R"))
source(here("Data Analysis", "Utility.R"))

########################################
## Load in data and mutate conditions ##
########################################
data = data.frame(animal = mice_ID, treatment = treatment, env = environment, latency = BFT) %>% 
  mutate(treatment = factor(treatment, levels = c("Sal", "CVAD"), labels = c("SAL", "CVAD")), env = factor(env, levels=c("PE", "EE")))


sum(as.numeric(data$treatment)-1) # Check how many animals CVAD (should be 12)
sum(as.numeric(data$env)-1)  # Check how many animals CVAD (should be 18)

data = data %>% filter(animal != "F13.2_B") # exclude the mouse that died 


######################
## path length freq ##
######################
ANOVA = two_way_aov(data, "treatment", "env", "latency")

con = emmeans(ANOVA[[3]], pairwise ~ treatment*env, adjust = "tukey")

############################
## Two‑way ANOVA Bayesian ##
############################
mod = cmdstan_model(here("Data Analysis", "Analysis tests", "BFT", "ANOVA_BFT.stan"))
path_PL = here(path, "fit_BFT.rds")

fit = two_way_bayes(mod, data, "latency", path_PL, 7)
############
## plots#### 
############
BFT= ggplot(data, aes(x = treatment, y = latency, fill = env)) +
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
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  #scale_x_discrete(labels = c("0" = "SAL", "1"= "CVAD"))+
  scale_fill_manual(values = c("PE" = "salmon2", "EE" = "cadetblue3"), name = "")+
  labs(x="Treatment", y = "Latency (s)")+
  theme_classic(base_size =20)
