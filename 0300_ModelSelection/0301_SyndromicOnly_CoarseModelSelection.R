library(dummies)
library(tidyverse)
library(future)
library(cmdstanr)
library(posterior)
library(furrr)
source("RATNeg_Workflow/0300_ModelSelection/0300_Functions.R") 
# Read in data frame
nasal_dat <- readRDS("Workflow/0100_DataPreparation/0100_nasal_dat.RDS")

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
covariate_nam <- c("age", "gender")
saveRDS(covariate_nam, "RATNeg_Workflow//0100_DataPreparation/0103_covariate_nam.rds")
# Symptoms List
symptom_nam <- nasal_dat %>% 
  select(-c(id, nasal_ag, all_of(covariate_nam), result, month, week, time_w, day)) %>% 
  names()
saveRDS(symptom_nam, "RATNeg_Workflow/0100_DataPreparation/0103_symptom_nam.rds")

## Nasal dat function
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

# Generate initial round of model selection candidate data lists
best_symptoms_so_far <- symptom_nam
round <- 14
assign(paste0("Nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round,
                      covariate_names = NULL))

omegas <- Nasal_Round14$omegas %>% select(starts_with("result"),CV ) %>%
  pivot_longer(-c(CV), names_to = "symptoms")

ggplot(omegas, aes(x = value, fill = symptoms)) +
  geom_histogram() +
  facet_grid(vars(CV), vars(symptoms)) +
  geom_vline(xintercept = 0)

reject_symptoms <- c("headache", "loss_of_smell", "runny_nose")

# Round 11 ----------------------------------------------------------------

best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
round <- 11
round == length(best_symptoms_so_far)
assign(paste0("Nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round))

omegas <- get(paste0("Nasal_Round", round))$omegas %>% select(starts_with("result"),CV ) %>%
  pivot_longer(-c(CV), names_to = "symptoms")

ggplot(omegas, aes(x = value, fill = symptoms)) +
  geom_histogram() +
  facet_grid(vars(CV), vars(symptoms)) +
  geom_vline(xintercept = 0)

reject_symptoms <- c("headache", "loss_of_smell", "runny_nose", "breath_problem",
                     "cough", "loss_of_taste", "sore_throat")

# Round 07 ----------------------------------------------------------------

best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
round <- 7
round == length(best_symptoms_so_far)
assign(paste0("Nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round))

saveRDS(get(paste0("Nasal_Round", round)), 
        paste0("RATNeg_Workflow/0300_ModelSelection/Output/Nasal_Round", 
               round, ".rds"))
omegas <- get(paste0("Nasal_Round", round))$omegas %>% select(starts_with("result"),CV ) %>%
  pivot_longer(-c(CV), names_to = "symptoms")

ggplot(omegas, aes(x = value, fill = symptoms)) +
  geom_histogram() +
  facet_grid(vars(CV), vars(symptoms)) +
  geom_vline(xintercept = 0)

reject_symptoms <- c("headache", "loss_of_smell", "runny_nose", "breath_problem",
                     "cough", "loss_of_taste", "sore_throat", "muscle_pain",
                     "red_eye", "vomit")

# Round 04 ----------------------------------------------------------------

best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
round <- 4
round == length(best_symptoms_so_far)
assign(paste0("Nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round))

saveRDS(get(paste0("Nasal_Round", round)), 
        paste0("RATNeg_Workflow/0300_ModelSelection/Output/Nasal_Round", 
               round, ".rds"))
omegas <- get(paste0("Nasal_Round", round))$omegas %>% select(starts_with("result"),CV ) %>%
  pivot_longer(-c(CV), names_to = "symptoms")

ggplot(omegas, aes(x = value, fill = symptoms)) +
  geom_histogram() +
  facet_grid(vars(CV), vars(symptoms)) +
  geom_vline(xintercept = 0)

reject_symptoms <- c("headache", "loss_of_smell", "runny_nose", "breath_problem",
                     "cough", "loss_of_taste", "sore_throat", "muscle_pain",
                     "red_eye", "vomit", "tired")

# Round 03 ----------------------------------------------------------------

best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
round <- 3
round == length(best_symptoms_so_far)
assign(paste0("Nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round))

saveRDS(get(paste0("Nasal_Round", round)), 
        paste0("RATNeg_Workflow/0300_ModelSelection/Output/Nasal_Round", 
               round, ".rds"))
omegas <- get(paste0("Nasal_Round", round))$omegas %>% select(starts_with("result"),CV ) %>%
  pivot_longer(-c(CV), names_to = "symptoms")

ggplot(omegas, aes(x = value, fill = symptoms)) +
  geom_histogram() +
  facet_grid(vars(CV), vars(symptoms)) +
  geom_vline(xintercept = 0)

reject_symptoms <- c("headache", "loss_of_smell", "runny_nose", "breath_problem",
                     "cough", "loss_of_taste", "sore_throat", "muscle_pain",
                     "red_eye", "vomit", "tired", "wet_cough")

# Round 02 ----------------------------------------------------------------

best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
round <- 2
round == length(best_symptoms_so_far)
assign(paste0("Nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round))

saveRDS(get(paste0("Nasal_Round", round)), 
        paste0("RATNeg_Workflow/0300_ModelSelection/Output/Nasal_Round", 
               round, ".rds"))
omegas <- get(paste0("Nasal_Round", round))$omegas %>% select(starts_with("result"),CV ) %>%
  pivot_longer(-c(CV), names_to = "symptoms")

ggplot(omegas, aes(x = value, fill = symptoms)) +
  geom_histogram() +
  facet_grid(vars(CV), vars(symptoms)) +
  geom_vline(xintercept = 0)

reject_symptoms <- c("headache", "loss_of_smell", "runny_nose", "breath_problem",
                     "cough", "loss_of_taste", "sore_throat", "muscle_pain",
                     "red_eye", "vomit", "tired", "wet_cough", "diarrhoea")

# Round 01 ----------------------------------------------------------------

best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
round <- 1
round == length(best_symptoms_so_far)
assign(paste0("Nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round))

saveRDS(get(paste0("Nasal_Round", round)), 
        paste0("RATNeg_Workflow/0300_ModelSelection/Output/Nasal_Round", 
               round, ".rds"))
omegas <- get(paste0("Nasal_Round", round))$omegas %>% select(starts_with("result"),CV ) %>%
  pivot_longer(-c(CV), names_to = "symptoms")

ggplot(omegas, aes(x = value, fill = symptoms)) +
  geom_histogram() +
  facet_grid(vars(CV), vars(symptoms)) +
  geom_vline(xintercept = 0)

reject_symptoms <- c("headache", "loss_of_smell", "runny_nose", "breath_problem",
                     "cough", "loss_of_taste", "sore_throat", "muscle_pain",
                     "red_eye", "vomit", "tired", "wet_cough", "diarrhoea", "fever")

# Round 0 ----------------------------------------------------------------

best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
round <- 0
round == length(best_symptoms_so_far)
assign(paste0("Nasal_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round))

saveRDS(get(paste0("Nasal_Round", round)), 
        paste0("RATNeg_Workflow/0300_ModelSelection/Output/Nasal_Round", 
               round, ".rds"))
omegas <- get(paste0("Nasal_Round", round))$omegas %>% select(starts_with("result"),CV ) %>%
  pivot_longer(-c(CV), names_to = "symptoms")

ggplot(omegas, aes(x = value, fill = symptoms)) +
  geom_histogram() +
  facet_grid(vars(CV), vars(symptoms)) +
  geom_vline(xintercept = 0)


# Compare -----------------------------------------------------------------

nas0 <- readRDS("RATNeg_Workflow/0300_ModelSelection/Output/nasal_Round0__.rds")
nas1 <- readRDS("RATNeg_Workflow/0300_ModelSelection/Output/nasal_Round1.rds")
nas2 <- readRDS("RATNeg_Workflow/0300_ModelSelection/Output/nasal_Round2__.rds")
nas3 <- readRDS("RATNeg_Workflow/0300_ModelSelection/Output/nasal_Round3.rds")
nas4 <- readRDS("RATNeg_Workflow/0300_ModelSelection/Output/nasal_Round4.rds")
nas7 <- readRDS("RATNeg_Workflow/0300_ModelSelection/Output/nasal_Round7.rds")
nas11 <- readRDS("RATNeg_Workflow/0300_ModelSelection/Output/nasal_Round11.rds")
nas14 <- readRDS("RATNeg_Workflow/0300_ModelSelection/Output/nasal_Round14.rds")


nas <- rbind(nas0$cv_log_loss,
             nas1$cv_log_loss,
             nas2$cv_log_loss,
             nas3$cv_log_loss,
             nas4$cv_log_loss,
             nas7$cv_log_loss,
             nas11$cv_log_loss,
             nas14$cv_log_loss)

ggplot(nas, aes(x = log_loss, fill = factor(FitType))) +
  # geom_histogram(position = "identity", alpha = 0.5) +
  geom_boxplot() +
  facet_grid(vars(FitType)) +
  coord_cartesian(xlim = c(0,1)) +
  ggtitle("Nasal Swabs")


