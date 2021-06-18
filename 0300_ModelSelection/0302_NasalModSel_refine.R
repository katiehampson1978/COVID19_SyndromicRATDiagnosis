library(dummies)
library(tidyverse)
library(future)
library(cmdstanr)
library(posterior)
library(furrr)
source("RATNeg_Workflow/0300_ModelSelection/0300_Functions.R") 
# Read in data frame
nasal_dat <- readRDS("RATNeg_Workflow/0100_DataPreparation/0100_nasal_dat.RDS") 

# Recode 
nasal_dat$age <- (nasal_dat$age - mean(nasal_dat$age))/sd(nasal_dat$age)
nasal_dat$time_w <- (nasal_dat$week - mean(nasal_dat$week))/sd(nasal_dat$week)
nasal_dat$gender <- dummy(nasal_dat$gender)[,1] # Make women baseline

# Handy variable generation ------------------------------------------------

# Our data are temporally structured so we will carry out cross-validation
# by removing 1 of 5 weeks of data at a time, fitting our model to the remaining
# 4 weeks of data and predicting for the removed week of data.

num_weeks <- length(unique(nasal_dat$week))

# Create 5 columns each of which has one cell marked, this will be the testing 
# set for that cross validation 
cross_val_tab <- matrix(0, nrow = num_weeks, ncol = num_weeks)
diag(cross_val_tab) <- 1
colnames(cross_val_tab) <- paste0("cv_", 1:num_weeks)

# Covariate list
covariate_nam <- c("time_w", "age", "gender")
covariate_combos <- list(covariate_nam[c(1:3)],
                         covariate_nam[c(1:2)],
                         covariate_nam[c(2:3)],
                         covariate_nam[c(1,3)],
                         covariate_nam[c(1)],
                         covariate_nam[c(2)],
                         covariate_nam[c(3)]
)

# Symptoms List
symptom_nam <- nasal_dat %>% 
  select(-c(id, nasal_ag, all_of(covariate_nam), result, month, week, time_w, day)) %>% 
  names()

## nasal dat function
nasal_tidy_run <- function(best_symptoms_so_far, round, covariate_names){
  tidy_run <- tidy_run_cv(dat = nasal_dat, 
              swab_type = "nasal", 
              covariate_names = covariate_names,
              symptom_nam = best_symptoms_so_far, 
              model_type = paste(round, 
                                  paste(best_symptoms_so_far, 
                                        collapse = "."),
                                 paste(covariate_names, 
                                       collapse = "."), sep = "_"),
              cross_val_table = cross_val_tab)
  saveRDS(tidy_run, 
          paste0("RATNeg_Workflow/0300_ModelSelection/Output/nasal_Round", 
                 paste(round, 
                       paste(best_symptoms_so_far, 
                             collapse = "."),
                       paste(covariate_names, 
                             collapse = "."), sep = "_"), ".rds"))
  tidy_run
}

# Top 4 symptoms
reject_symptoms <- c("headache", "loss_of_smell", "runny_nose", "breath_problem",
                     "cough", "loss_of_taste", "sore_throat", "muscle_pain",
                     "red_eye", "vomit")



# A for loop is justified here as the number of cores available means the models
# must run sequentially if CV and chains are to be parallelised
for(i in seq_along(covariate_combos)){
covariate_nam <- covariate_combos[[i]]

round <- 4
best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
assign(paste0("nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round, covariate_names = covariate_nam))

# Round 3 -----------------------------------------------------------------

round <- 3
best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
best_symptoms_so_far <- best_symptoms_so_far[1:3]
assign(paste0("nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round, covariate_names = covariate_nam))


best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
best_symptoms_so_far <- best_symptoms_so_far[2:4]
assign(paste0("nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round, covariate_names = covariate_nam))

best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
best_symptoms_so_far <- best_symptoms_so_far[c(1,3:4)]
assign(paste0("nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round, covariate_names = covariate_nam))


best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
best_symptoms_so_far <- best_symptoms_so_far[c(1:2,4)]
assign(paste0("nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round, covariate_names = covariate_nam))


# Round 2 -----------------------------------------------------------------

round <- 2
best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
best_symptoms_so_far <- best_symptoms_so_far[1:2]
assign(paste0("nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round, covariate_names = covariate_nam))


best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
best_symptoms_so_far <- best_symptoms_so_far[2:3]
assign(paste0("nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round, covariate_names = covariate_nam))


best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
best_symptoms_so_far <- best_symptoms_so_far[3:4]
assign(paste0("nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round, covariate_names = covariate_nam))


best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
best_symptoms_so_far <- best_symptoms_so_far[c(1,3)]
assign(paste0("nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round, covariate_names = covariate_nam))


best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
best_symptoms_so_far <- best_symptoms_so_far[c(1,4)]
assign(paste0("nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round, covariate_names = covariate_nam))


best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
best_symptoms_so_far <- best_symptoms_so_far[c(2,4)]
assign(paste0("nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round, covariate_names = covariate_nam))


# Round 1 -----------------------------------------------------------------

round <- 1

best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
best_symptoms_so_far <- best_symptoms_so_far[1]
assign(paste0("nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round, covariate_names = covariate_nam))


best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
best_symptoms_so_far <- best_symptoms_so_far[2]
assign(paste0("nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round, covariate_names = covariate_nam))


best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
best_symptoms_so_far <- best_symptoms_so_far[3]
assign(paste0("nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round, covariate_names = covariate_nam))


best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
best_symptoms_so_far <- best_symptoms_so_far[4]
assign(paste0("nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round, covariate_names = covariate_nam))


}


# Result Only -------------------------------------------------------------


round <- 0
assign(paste0("nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = NULL,
                      round = round, covariate_names = NULL))

covariate_nam <- c("time_w", "age", "gender")
covariate_combos <- list(covariate_nam[c(1:3)],
                         covariate_nam[c(1:2)],
                         covariate_nam[c(2:3)],
                         covariate_nam[c(1,3)],
                         covariate_nam[c(1)],
                         covariate_nam[c(2)],
                         covariate_nam[c(3)]
)

round <- 0
for(i in seq_along(covariate_combos)){
covariate_nam <-  covariate_combos[[i]]
assign(paste0("nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = NULL,
                      round = round, covariate_names = covariate_nam))


}
assign(paste0("nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = NULL,
                      round = round, covariate_names = NULL))
