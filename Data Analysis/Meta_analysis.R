pacman::p_load(tidyverse, here, gt)

#################################################################
####### Import data #############################################
#################################################################

data_nodig = read.csv(here("Data", "Data.csv"))


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

