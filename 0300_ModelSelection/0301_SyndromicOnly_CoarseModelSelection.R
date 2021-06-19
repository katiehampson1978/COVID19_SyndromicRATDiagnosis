#' In this script, we conduct our broadest level model selection to determine 
#' which symptom combinations are most predictive for the syndromic only model.
#' We are exploring a large model space here with the aim of narrowing down to
#' the 5 best predictive models in this class. To do this (without fitting every
#' possible model which would both put us at risk of overfitting, even with the
#' demanding temporal cross validation we are subjecting the model too, and would
#' not be computationally feasible when we need a fast answer), we fit the most
#' complex model we have, select the symptoms with non-zero correlations with 
#' COVID-19 status, and repeat. We use the same covariates (age and gender) for 
#' all model fits here (again to reduce the model space we need to explore). This
#' isn't statistically state of the art but is a quick and robust route to model
#' selection that is easily replicated.

# Source helper code and libraries
source("0000_HelperCode_Libraries/0001_Libraries.R")
source("0000_HelperCode_Libraries/0002_GlobalStanFittingParams.R")
source("0000_HelperCode_Libraries/0003_HelperFunctions.R")
# Read in data frame
nasal_dat <- readRDS("0100_Data/0103_nasal_dat.RDS")

# Recode and standardise covariates
nasal_dat$age <- (nasal_dat$age - mean(nasal_dat$age))/sd(nasal_dat$age)
nasal_dat$gender <- dummy(nasal_dat$gender)[,1] # Make women baseline gender 
# note that there is a downstream package that causes a warning message here, it
# is annoying but doesn't cause any problems.

# Handy variable generation ------------------------------------------------

# Our data are temporally structured so we will carry out temporal cross-validation
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
# Symptoms List
symptom_nam <- nasal_dat %>% 
  select(-c(id, nasal_ag, all_of(covariate_nam), result, month, week, day)) %>% 
  names()

## Nasal dat function
nasal_tidy_run <- function(best_symptoms_so_far, round){
  tidy_run <- tidy_run_cv(dat = nasal_dat, 
                          swab_type = "nasal", 
                          covariate_names = covariate_nam,
                          symptom_nam = best_symptoms_so_far, 
                          model_type = paste(round, 
                                             paste(best_symptoms_so_far, 
                                                   collapse = "."),
                                             paste(covariate_nam, 
                                                   collapse = "."), sep = "_"),
                          cross_val_table = cross_val_tab, 
                          syndromic_only = TRUE)
  saveRDS(tidy_run, 
          paste0("0300_ModelSelection/Output/SyndromicOnly_Coarse_Round", 
                 paste(round, 
                       paste(best_symptoms_so_far, 
                             collapse = "."),
                       paste(covariate_nam, 
                             collapse = "."), sep = "_"), ".rds"))
  tidy_run
}

# Generate initial round of model selection candidate data lists
best_symptoms_so_far <- symptom_nam
# Number of symptoms
round <- 14
# Sense check!
round == length(best_symptoms_so_far)
# Fit model 
assign(paste0("SyndromicOnly_Coarse_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round))

omegas <- get(paste0("SyndromicOnly_Coarse_Round", round))$omegas %>% 
  select(starts_with("result"),CV ) %>%
  pivot_longer(-c(CV), names_to = "symptoms")

ggplot(omegas, aes(x = value, fill = symptoms)) +
  geom_histogram() +
  facet_grid(vars(CV), vars(symptoms)) +
  geom_vline(xintercept = 0)

reject_symptoms <- c("cough", "diarrhoea", "headache", "vomit")

# Round 10 ----------------------------------------------------------------

best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
round <- 10
round == length(best_symptoms_so_far)
assign(paste0("SyndromicOnly_Coarse_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round))

omegas <- get(paste0("SyndromicOnly_Coarse_Round", round))$omegas %>% select(starts_with("result"),CV ) %>%
  pivot_longer(-c(CV), names_to = "symptoms")

ggplot(omegas, aes(x = value)) +
  geom_histogram() +
  facet_grid(vars(CV), vars(symptoms)) +
  geom_vline(xintercept = 0)

reject_symptoms <- c("cough", "diarrhoea", "headache", "vomit", 
                     "loss_of_taste", "muscle_pain", "red_eye", "tired")

# Round 06 ----------------------------------------------------------------

best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
round <- 6
round == length(best_symptoms_so_far)
assign(paste0("SyndromicOnly_Coarse_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round))

omegas <- get(paste0("SyndromicOnly_Coarse_Round", round))$omegas %>% select(starts_with("result"),CV ) %>%
  pivot_longer(-c(CV), names_to = "symptoms")

ggplot(omegas, aes(x = value)) +
  geom_histogram() +
  facet_grid(vars(CV), vars(symptoms)) +
  geom_vline(xintercept = 0)

reject_symptoms <- c("cough", "diarrhoea", "headache", "vomit", 
                     "loss_of_taste", "muscle_pain", "red_eye", "tired",
                     "fever", "sore_throat")

# Round 04 ----------------------------------------------------------------

best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
round <- 4
round == length(best_symptoms_so_far)
assign(paste0("SyndromicOnly_Coarse_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round))


omegas <- get(paste0("SyndromicOnly_Coarse_Round", round))$omegas %>% 
  select(starts_with("result"),CV ) %>%
  pivot_longer(-c(CV), names_to = "symptoms")

ggplot(omegas, aes(x = value)) +
  geom_histogram() +
  facet_grid(vars(CV), vars(symptoms)) +
  geom_vline(xintercept = 0)

reject_symptoms <- c("cough", "diarrhoea", "headache", "vomit", 
                     "loss_of_taste", "muscle_pain", "red_eye", "tired",
                     "fever", "sore_throat", "breath_problem")

# Round 03 ----------------------------------------------------------------

best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
round <- 3
round == length(best_symptoms_so_far)
assign(paste0("SyndromicOnly_Coarse_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round))


omegas <- get(paste0("SyndromicOnly_Coarse_Round", round))$omegas %>% 
  select(starts_with("result"),CV ) %>%
  pivot_longer(-c(CV), names_to = "symptoms")

ggplot(omegas, aes(x = value)) +
  geom_histogram() +
  facet_grid(vars(CV), vars(symptoms)) +
  geom_vline(xintercept = 0)

reject_symptoms <- c("cough", "diarrhoea", "headache", "vomit", 
                     "loss_of_taste", "muscle_pain", "red_eye", "tired",
                     "fever", "sore_throat", "breath_problem", "runny_nose")

# Round 02 ----------------------------------------------------------------

best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
round <- 2
round == length(best_symptoms_so_far)
assign(paste0("SyndromicOnly_Coarse_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round))

omegas <- get(paste0("SyndromicOnly_Coarse_Round", round))$omegas %>% select(starts_with("result"),CV ) %>%
  pivot_longer(-c(CV), names_to = "symptoms")

ggplot(omegas, aes(x = value)) +
  geom_histogram() +
  facet_grid(vars(CV), vars(symptoms)) +
  geom_vline(xintercept = 0)

reject_symptoms <- c("cough", "diarrhoea", "headache", "vomit", 
                     "loss_of_taste", "muscle_pain", "red_eye", "tired",
                     "fever", "sore_throat", "breath_problem", "runny_nose",
                     "wet_cough")

# Round 01 ----------------------------------------------------------------

best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
round <- 1
round == length(best_symptoms_so_far)
assign(paste0("SyndromicOnly_Coarse_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round))


omegas <- get(paste0("SyndromicOnly_Coarse_Round", round))$omegas %>% select(starts_with("result"),CV ) %>%
  pivot_longer(-c(CV), names_to = "symptoms")

ggplot(omegas, aes(x = value, fill = symptoms)) +
  geom_histogram() +
  facet_grid(vars(CV), vars(symptoms)) +
  geom_vline(xintercept = 0)

reject_symptoms <- c("cough", "diarrhoea", "headache", "vomit", 
                     "loss_of_taste", "muscle_pain", "red_eye", "tired",
                     "fever", "sore_throat", "breath_problem", "runny_nose",
                     "wet_cough", "loss_of_smell")

# Round 0 ----------------------------------------------------------------

best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
round <- 0
round == length(best_symptoms_so_far)
assign(paste0("SyndromicOnly_Coarse_Round", round), 
       nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                      round = round))


# Compare -----------------------------------------------------------------

nas0 <- readRDS("0300_ModelSelection/Output/SyndromicOnly_Coarse_Round0__age.gender.rds")
nas1 <- readRDS("0300_ModelSelection/Output/SyndromicOnly_Coarse_Round1_loss_of_smell_age.gender.rds")
nas2 <- readRDS("0300_ModelSelection/Output/SyndromicOnly_Coarse_Round2_wet_cough.loss_of_smell_age.gender.rds")
nas3 <- readRDS("0300_ModelSelection/Output/SyndromicOnly_Coarse_Round3_wet_cough.runny_nose.loss_of_smell_age.gender.rds")
nas4 <- readRDS("0300_ModelSelection/Output/SyndromicOnly_Coarse_Round4_breath_problem.wet_cough.runny_nose.loss_of_smell_age.gender.rds")
nas6 <- readRDS("0300_ModelSelection/Output/SyndromicOnly_Coarse_Round6_fever.breath_problem.wet_cough.sore_throat.runny_nose.loss_of_smell_age.gender.rds")
nas10 <- readRDS("0300_ModelSelection/Output/SyndromicOnly_Coarse_Round10_fever.breath_problem.wet_cough.sore_throat.tired.runny_nose.loss_of_smell.muscle_pain.red_eye.loss_of_taste_age.gender.rds")
nas14 <- readRDS("0300_ModelSelection/Output/SyndromicOnly_Coarse_Round14_fever.headache.cough.breath_problem.wet_cough.sore_throat.tired.runny_nose.loss_of_smell.muscle_pain.diarrhoea.vomit.red_eye.loss_of_taste_age.gender.rds")


synd_coarse <- rbind(nas0$cv_log_loss,
             nas1$cv_log_loss,
             nas2$cv_log_loss,
             nas3$cv_log_loss,
             nas4$cv_log_loss,
             nas6$cv_log_loss,
             nas10$cv_log_loss,
             nas14$cv_log_loss)

synd_coarse_wide <- synd_coarse %>% 
  pivot_wider(names_from = CV, values_from = log_loss, values_fn = list) %>%
  unnest(c(`1`, `2` , `3`, `4`,  `5`))

synd_coarse_wide <- synd_coarse_wide %>% 
  mutate(ModelLogLoss = synd_coarse_wide %>%
           select(c(`1`, `2` , `3`,  `4`, `5`)) %>%
           rowSums())


synd_coarse_wide$SympNum <- parse_number(synd_coarse_wide$FitType)
ggplot(synd_coarse_wide, aes(x = ModelLogLoss)) +
  geom_boxplot() +
  facet_grid(vars(SympNum)) +
  ggtitle("Nasal Swabs") +
  theme(legend.position = "none") 


