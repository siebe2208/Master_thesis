pacman::p_load(tidyverse, here, gt,zoo,ggplot2,lme4,lmerTest,r2glmm, logspline, cmdstanr, posterior)
source(here("Data Analysis", "weight_curve_script.R"))
source(here("Data Analysis", "unclean_data.R"))
source(here("Data Analysis", "Utility.R"))
#################################################################
####### Import data #############################################
#################################################################

data_nodig = read.csv(here("Data", "Data.csv"))

data_BW = bind_cols(data_F11.2, data_F11.3, data_F12, data_F13.2,data_F14,data_M13,data_M14,data_M15.2,data_M16) %>% 
  mutate(DAS = days_after...1) %>%select(-matches("days_after")) %>% mutate(across(-DAS, ~ na.spline(.x, x = DAS, na.rm = FALSE))) %>% 
  pivot_longer(cols = -DAS, names_to="animal",values_to = "weight") %>% mutate(treatment = rep(forcats::fct_rev(treatment),26), env = rep(factor(environment),26)) 

#################################################################
####### Meta-Analysis ###########################################
#################################################################

n = ncol(data_BW[ , !(names(data_BW) %in% "days_after") ])      # Sample size
paste0("Sample size: ", n)

## Data males
subset_M = data_BW %>% select(!days_after) %>% select(contains("M")) # Sample size males
n_males = ncol(subset_M)

paste0("Sample size males: ", n_males)

subset_M %>% slice(1) %>% summarise( 
mean_BW_M = mean(c_across(everything())),    # Starting weight M       
  SD_BW_M = sd(c_across(everything()))
)

## Data females
subset_F = data_BW %>% select(!days_after) %>% select(contains("F")) # Sample size females
n_females = ncol(subset_F)

paste0("Sample size males: ", n_females)

subset_F %>% slice(1) %>% summarise(
  mean_BW_F = mean(c_across(everything())),     #Starting weight females
  SD_BW_F = sd(c_across(everything()))
)

#################################################################
####### table counts ###########################################
#################################################################

table_FM = data_nodig %>% 
  count(treatment, env, sex) %>% pivot_wider(names_from = sex,values_from = n,values_fill = 0) %>% 
  mutate(PE_F = ifelse(env == "PE", F, NA),PE_M = ifelse(env == "PE", M, NA),EE_F = ifelse(env == "EE", F, NA),EE_M = ifelse(env == "EE", M, NA)) %>% 
  mutate(treatment = ifelse(treatment == "Sal", "Saline", treatment)) %>% group_by(treatment) %>% summarise(
    PE_F = sum(PE_F, na.rm = TRUE),
    PE_M = sum(PE_M, na.rm = TRUE),
    EE_F = sum(EE_F, na.rm = TRUE),
    EE_M = sum(EE_M, na.rm = TRUE))


table_sample = table_FM %>% gt() %>% tab_spanner(label = "Male", columns = c(PE_M, EE_M)) %>%
  tab_spanner(label = "Female", columns = c(PE_F, EE_F)) %>% 
  cols_label(treatment = "Treatment",PE_F = "PE",PE_M = "PE",EE_F = "EE",EE_M = "EE") %>% 
  cols_align(align = "center",columns = c(PE_F, EE_F, PE_M, EE_M)) %>% 
  tab_style(style = cell_text(size = px(16)),   locations = cells_title(groups = "title")) %>% 
  tab_style(style = cell_text(size = px(13)),locations = cells_source_notes()) %>% 
  tab_style(style = cell_text(style = "italic"),locations = cells_body(columns = treatment,rows = treatment %in% c("CVAD", "Saline"))) %>% 
  tab_options(column_labels.border.top.color = "black",column_labels.border.bottom.color = "black",
              table_body.border.bottom.color = "black", table_body.hlines.color = "transparent",
              table.border.top.color = "transparent",table.border.bottom.color = "transparent",
              heading.border.bottom.color = "black", heading.align = "left", table.width = px(700), table.font.names = "Times new Roman") %>% 
  tab_style(style = cell_text(style = "italic"),locations = cells_column_labels(columns = c(PE_F, EE_F, PE_M, EE_M))) %>% 
  tab_header(title = html("<b>Table 1</b><div style='margin-top:12px;'><i>Distribution of mice across treatment, environment, and sex.</i></div>")) %>% # div to insure according to APA spacing
  tab_source_note(source_note = md("*Note.* Values represent sample sizes (n) of pseudo‑randomly assigned mice across treatment, environment, and sex.
                                    PE = poor environment; EE = enriched environment; CVAD = chemotherapy; 
                                   Saline = control treatment. Due to limited male offspring in these cohorts, data from male and female mice were combined
                                   in subsequent analyses to maximize statistical power."))

gtsave(table_sample, "table_sample.png", zoom = 3) # increase DPI 

#################################################################
####### Analysis weights        #################################
#################################################################
mod = lmer(weight ~DAS * env * treatment +(1|animal), data = data_BW)
summary(mod)

r2beta(mod, method = "nsj")


data_bayes = data_BW %>% ungroup() %>% mutate(animal = rep(1:35, 26)) %>% arrange(animal) %>% 
  mutate(treatment = as.numeric(treatment)-1, env = as.numeric(env)-1)
mod = cmdstan_model(here("Data Analysis", "Analysis tests", "MWM", "ANOVA_MWM.stan"))


z = as.numeric(scale(data_bayes$weight))
stan_data = list(N =nrow(data_bayes), J = 35, Day  = data_bayes$DAS, animal = data_bayes$animal, 
                 treatment = data_bayes$treatment, environment=data_bayes$env, y = z)

fit = mod$sample(data = stan_data,iter_warmup = 1000,iter_sampling = 2000,chains = 4,parallel_chains = 4)

fit$summary()



prior = as.numeric(fit$draws("prior"))

draws = fit$draws(c('beta_day', 'beta_env', 'beta_t', 'beta_d_e', 'beta_d_t', 'beta_e_t', 'beta_int')) %>% as_draws_df() %>% 
  select(-.chain, -.iteration, -.draw)


BFs = draws %>% summarise(across(everything(), ~ get_bf(prior, .x)))


#################################################################
####### plots         ###########################################
#################################################################
ggplot(data_BW) + stat_summary(aes(x=DAS,y=weight,color = env, linetype = treatment,group=interaction(treatment, env)),
                                 fun = mean, geom = "line", linewidth =1.2)+
  stat_summary(aes(x = DAS, y= weight, color = env, group = interaction(env, treatment)),
               fun = "mean",
               geom = "point",
               size = 3)+ labs(x="Days after separation", y = "Weight (g)")+theme_classic(base_size = 20)+
  scale_x_continuous(breaks = c(1,7,14,21,26))+
  scale_color_manual(values = c("PE" = "salmon2", "EE" = "cadetblue3"), name = "Environment")+
  scale_linetype_discrete(name = "Treatment",labels = c("SAL", "CVAD"))

sum = data_BW %>% group_by(treatment, env,DAS) %>% summarise(mean = mean(weight))

scale_x_continuous(breaks = sort(unique(sum$Day)))+theme_classic(base_size = 20)+
  scale_color_manual(values = c("0" = "salmon2", "1" = "cadetblue3"),labels = c("0" = "PE", "1" = "EE"), name = "Environment")+
  labs(x="Day", y = "Distance to target (cm)")+
  scale_linetype_discrete(name = "Treatment",labels = c("SAL", "CVAD"))+
  stat_summary(aes(x = Day, y= mean, color = env, group = interaction(env, treatment)),
               fun = "mean",
               geom = "point",
               size = 3) 

