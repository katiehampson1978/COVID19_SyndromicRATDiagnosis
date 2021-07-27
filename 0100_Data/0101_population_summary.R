# Load libraries
source("0000_HelperCode_Libraries/0001_Libraries.R")
# Import Data
dat <- readRDS("0100_Data/0103_nasal_dat.RDS")

# Calculate population characteristics
pop_summ <- dat %>% 
  group_by(gender) %>% 
  summarise("Gender_Count" = summary(gender),
            "Age_Mean" = mean(age),
            "Age_SD" = sd(age)) %>%
  filter(!Gender_Count == 0) # Remove annoying duplications

saveRDS(pop_summ, "0100_Data/0103_pop_summ.RDS")
