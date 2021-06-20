#' Now that we have narrowed down to the best set of symptoms for the syndromic
#' only model. This has reduced the potential model space to manageable level
#' so we can now explore the symptom-covariate combinations to the 5 best models
#' of each complexity level. 

# Source helper code and libraries
source("0000_HelperCode_Libraries/0001_Libraries.R")
source("0000_HelperCode_Libraries/0002_GlobalStanFittingParams.R")
source("0000_HelperCode_Libraries/0003_HelperFunctions.R")
# Read in data frame
nasal_dat <- readRDS("0100_Data/0103_nasal_dat.RDS")

# Recode 
nasal_dat$age <- (nasal_dat$age - mean(nasal_dat$age))/sd(nasal_dat$age)
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
covariate_combos <- list(covariate_nam[c(1:2)],
                         covariate_nam[c(1)],
                         covariate_nam[c(2)]
)

# Symptoms List
symptom_nam <- nasal_dat %>% 
  select(-c(id, nasal_ag, all_of(covariate_nam), result, month, week, day)) %>% 
  names()

## Nasal dat function
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
                          cross_val_table = cross_val_tab, 
                          syndromic_only = FALSE)
  saveRDS(tidy_run, 
          paste0("0300_ModelSelection/Output/SyndromicRAT_Fine_Round", 
                 paste(round, 
                       paste(best_symptoms_so_far, 
                             collapse = "."),
                       paste(covariate_names, 
                             collapse = "."), sep = "_"), ".rds"))
  tidy_run
}

# Top 4 symptoms
reject_symptoms <- c("loss_of_smell", "loss_of_taste", "muscle_pain", 
                     "red_eye", "runny_nose", "cough", "headache",
                     "vomit", "breath_problem", "sore_throat")



# A for loop is justified here as the number of cores available means the models
# must run sequentially if CV and chains are to be parallelised
for(i in seq_along(covariate_combos)){
  covariate_nam <- covariate_combos[[i]]
  
  round <- 4
  best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
  assign(paste0("SyndromicRAT_Coarse_Round", round), 
         nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                        round = round, covariate_names = covariate_nam))
  
  # Round 3 -----------------------------------------------------------------
  
  round <- 3
  best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
  best_symptoms_so_far <- best_symptoms_so_far[1:3]
  assign(paste0("SyndromicRAT_Coarse_Round", round), 
         nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                        round = round, covariate_names = covariate_nam))
  
  
  best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
  best_symptoms_so_far <- best_symptoms_so_far[2:4]
  assign(paste0("SyndromicRAT_Coarse_Round", round), 
         nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                        round = round, covariate_names = covariate_nam))
  
  best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
  best_symptoms_so_far <- best_symptoms_so_far[c(1,3:4)]
  assign(paste0("SyndromicRAT_Coarse_Round", round), 
         nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                        round = round, covariate_names = covariate_nam))
  
  
  best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
  best_symptoms_so_far <- best_symptoms_so_far[c(1:2,4)]
  assign(paste0("SyndromicRAT_Coarse_Round", round), 
         nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                        round = round, covariate_names = covariate_nam))
  
  
  # Round 2 -----------------------------------------------------------------
  
  round <- 2
  best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
  best_symptoms_so_far <- best_symptoms_so_far[1:2]
  assign(paste0("SyndromicRAT_Coarse_Round", round), 
         nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                        round = round, covariate_names = covariate_nam))
  
  
  best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
  best_symptoms_so_far <- best_symptoms_so_far[2:3]
  assign(paste0("SyndromicRAT_Coarse_Round", round), 
         nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                        round = round, covariate_names = covariate_nam))
  
  
  best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
  best_symptoms_so_far <- best_symptoms_so_far[3:4]
  assign(paste0("SyndromicRAT_Coarse_Round", round), 
         nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                        round = round, covariate_names = covariate_nam))
  
  
  best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
  best_symptoms_so_far <- best_symptoms_so_far[c(1,3)]
  assign(paste0("SyndromicRAT_Coarse_Round", round), 
         nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                        round = round, covariate_names = covariate_nam))
  
  
  best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
  best_symptoms_so_far <- best_symptoms_so_far[c(1,4)]
  assign(paste0("SyndromicRAT_Coarse_Round", round), 
         nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                        round = round, covariate_names = covariate_nam))
  
  
  best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
  best_symptoms_so_far <- best_symptoms_so_far[c(2,4)]
  assign(paste0("SyndromicRAT_Coarse_Round", round), 
         nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                        round = round, covariate_names = covariate_nam))
  
  
  # Round 1 -----------------------------------------------------------------
  
  round <- 1
  
  best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
  best_symptoms_so_far <- best_symptoms_so_far[1]
  assign(paste0("SyndromicRAT_Coarse_Round", round), 
         nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                        round = round, covariate_names = covariate_nam))
  
  
  best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
  best_symptoms_so_far <- best_symptoms_so_far[2]
  assign(paste0("SyndromicRAT_Coarse_Round", round), 
         nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                        round = round, covariate_names = covariate_nam))
  
  
  best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
  best_symptoms_so_far <- best_symptoms_so_far[3]
  assign(paste0("SyndromicRAT_Coarse_Round", round), 
         nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                        round = round, covariate_names = covariate_nam))
  
  
  best_symptoms_so_far <- symptom_nam[!symptom_nam %in% reject_symptoms]
  best_symptoms_so_far <- best_symptoms_so_far[4]
  assign(paste0("SyndromicRAT_Coarse_Round", round), 
         nasal_tidy_run(best_symptoms_so_far = best_symptoms_so_far,
                        round = round, covariate_names = covariate_nam))
  

# Round 0 -----------------------------------------------------------------

  assign(paste0("SyndromicRAT_Coarse_Round", 0), 
         nasal_tidy_run(best_symptoms_so_far = NULL,
                        round = 0, covariate_names = covariate_nam))
  
  
}



