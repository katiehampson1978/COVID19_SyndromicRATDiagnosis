library(dummies)
library(tidyverse)
library(future)
library(cmdstanr)
library(posterior)
library(furrr)
source("RATNeg_Workflow/0300_ModelSelection/0300_Functions.R") 
# Read in data frame
nasal_dat <- readRDS("Workflow/0100_DataPreparation/0100_nasal_dat.RDS") %>%
  filter(nasal_ag == 0)
saliv_dat <- readRDS("Workflow/0100_DataPreparation/0100_saliv_dat.RDS") %>%
  filter(sal_ag == 0)

# Recode 
nasal_dat$age <- (nasal_dat$age - mean(nasal_dat$age))/sd(nasal_dat$age)
saliv_dat$age <- (saliv_dat$age - mean(saliv_dat$age))/sd(saliv_dat$age)
nasal_dat$time_w <- (nasal_dat$week - mean(nasal_dat$week))/sd(nasal_dat$week)
saliv_dat$time_w <- (saliv_dat$week - mean(saliv_dat$week))/sd(saliv_dat$week)
nasal_dat$gender <- dummy(nasal_dat$gender)[,1] # Make women baseline
saliv_dat$gender <- dummy(saliv_dat$gender)[,1] # Make women baseline

# Handy variable generation ------------------------------------------------

# Our data are temporally structured so we will carry out cross-validation
# by removing 1 of 5 weeks of data at a time, fitting our model to the remaining
# 4 weeks of data and predicting for the removed week of data.

num_weeks <- length(unique(nasal_dat$week))
if(num_weeks != length(unique(saliv_dat$week))){
  print("Oh no! Different numbers of weeks for datasets")
}
# Create 5 columns each of which has one cell marked, this will be the testing 
# set for that cross validation 
cross_val_tab <- matrix(0, nrow = num_weeks, ncol = num_weeks)
diag(cross_val_tab) <- 1
colnames(cross_val_tab) <- paste0("cv_", 1:num_weeks)

# Covariate list
covariate_nam <- c("time_w", "age", "gender")
saveRDS(covariate_nam, "RATNeg_Workflow//0100_DataPreparation/0103_covariate_nam.rds")
# Symptoms List
symptom_nam <- saliv_dat %>% 
  select(-c(id, sal_ag, all_of(covariate_nam), result, month, week, time_w, day)) %>% 
  names()
saveRDS(symptom_nam, "RATNeg_Workflow/0100_DataPreparation/0103_symptom_nam.rds")

# Generate initial round of model selection candidate data lists

standat_saliv_all <- stan_dat_prep(dat = saliv_dat, cross_val_table = cross_val_tab, 
                       swab_type = "saliv", covariate_names = covariate_nam, 
                       response_names = c("result", symptom_nam), model_type = "MV")

fit_MVcv_all <- run_MV_cv_model_parr(standat = standat_saliv_all, n_cv_folds = 5,
                                     swab_type = "saliv", fit_type = "MV",
                                     modsel = TRUE)

output_list_complex <- list()
for(i in 1:5){
  output_list_complex[[i]] <- as_draws_df(fit_MVcv_all[[i]]$draws(variable = 
                                                                                    c("beta", "Omega", 
                                                                                      "logprob", "cv_fold",
                                                                                      "log_loss",  
                                                                                      "y_true")))
}

tidy_MVcv_complex <- tidy_output_MV(output_list = output_list_complex, 
                                swab_type = "saliv", 
                                fit_type = "MV_complex",
                                pred_names = covariate_nam,
                                response_names = c("result", symptom_nam), 
                                n_cv_folds = 5)
#---
standat_saliv_simple <- stan_dat_prep(dat = saliv_dat, cross_val_table = cross_val_tab, 
                                   swab_type = "saliv", covariate_names = covariate_nam, 
                                   response_names = c("result", "diarrhoea"), model_type = "MV")

fit_MVcv_simple <- run_MV_cv_model_parr(standat = standat_saliv_simple, n_cv_folds = 5,
                                     swab_type = "saliv", fit_type = "MV",
                                     modsel = TRUE)

output_list_simple <- list()
for(i in 1:5){
  output_list_simple[[i]] <- fit_MVcv_simple[[i]]$draws() %>% as_draws_df(fit$draws(variable = 
                                                                                    c("beta", "Omega", 
                                                                                      "logprob", "cv_fold",
                                                                                      "log_loss",  
                                                                                      "y_true")))
}

tidy_MVcv_simple <- tidy_output_MV(output_list = output_list_simple, 
                                swab_type = "saliv", 
                                fit_type = "MV_simple",
                                pred_names = covariate_nam,
                                response_names = c("result", "diarrhoea"), 
                                n_cv_folds = 5)

simple <- tidy_MVcv_simple$cv_log_loss
complex <- tidy_MVcv_complex$cv_log_loss

comp <- rbind(simple, complex)

ggplot(comp, aes(x = log_loss, fill = FitType)) + 
  geom_histogram(position = "identity", alpha = 0.6) +
  facet_grid(vars(CV))
