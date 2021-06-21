#' Holistic model performance can be measured using cross entropy (aka log loss)
 
# Load libraries and helper functions
source("0000_HelperCode_Libraries/0001_Libraries.R")
source("0000_HelperCode_Libraries/0003_HelperFunctions.R")

# Read in dataframe of best models
best_models_syndonly <- readRDS("0300_ModelSelection/0310_SyndromicOnly_BestModels.rds")
best_models_syndrat <- readRDS("0300_ModelSelection/0310_SyndromicRAT_BestModels.rds")

# rds files for best models
best_model_files_syndonly <- paste0("0300_ModelSelection/Output/SyndromicOnly_Fine_Round", 
                                    best_models_syndonly$ModelName, ".rds")
best_model_files_syndrat  <- paste0("0300_ModelSelection/Output/SyndromicRAT_Fine_Round",
                           best_models_syndrat $ModelName, ".rds")

## Read each file in parallel, extract log loss df, bind into single data frame,
## calculate model level log loss
# Syndromic only models
best_valid_syndonly <- best_model_files_syndonly %>% 
  future_map_dfr(read_logloss) %>% 
  pivot_wider(names_from = CV, values_from = log_loss, values_fn = list) %>%
  unnest(c(`1`, `2` , `3`, `4`,  `5`))  
best_valid_syndonly <- best_valid_syndonly %>%  
  mutate(ModelLogLoss = best_valid_syndonly %>%
           select(c(`1`, `2` , `3`,  `4`, `5`)) %>%
           rowSums())
best_valid_syndonly$ModelClass  <- "SyndOnly"
# Syndromic plus RAT models
best_valid_syndrat <- best_model_files_syndrat %>% 
  future_map_dfr(read_logloss)%>% 
  pivot_wider(names_from = CV, values_from = log_loss, values_fn = list) %>%
  unnest(c(`1`, `2` , `3`, `4`,  `5`))  
best_valid_syndrat <- best_valid_syndrat %>%  
  mutate(ModelLogLoss = best_valid_syndrat %>%
           select(c(`1`, `2` , `3`,  `4`, `5`)) %>%
           rowSums())
best_valid_syndrat$ModelClass  <- "SyndRAT"
# RAT only model
RATonly <- readRDS("0100_Data/0103_nasal_dat.RDS")
best_valid_RATonly <- data.frame("SwabType" = "nasal",
                                 "FitType" = "RATonly",
                                 "log_loss" = NA,
                                 "ModelLogLoss" = LogLoss(y_pred = RATonly$nasal_ag,
                                                             y_true = RATonly$result),
                                 "Iter" = NA,
                                 "Chain" = NA,
                                 "ModelClass" = "RATonly"
                                 )
  
# Bind into one data frame  
best_valid <- rbind(best_valid_syndonly %>% select(FitType, ModelLogLoss, ModelClass),
                    best_valid_syndrat %>% select(FitType, ModelLogLoss, ModelClass),
                    best_valid_RATonly %>% select(FitType, ModelLogLoss, ModelClass)) 
# Tidy names
best_valid$FitType <- paste0(parse_number(best_valid$FitType), 
                              "Symptom", "+Age")
best_valid$FitType[20001] <- "RATonly"
best_valid$FitType <- paste(best_valid$ModelClass, best_valid$FitType, sep = "_")

# Plot
ggplot(best_valid %>% group_by(ModelClass), 
       aes(x = ModelLogLoss, y = FitType, fill = FitType)) +
  geom_boxplot() 

best_valid %>% group_by(FitType) %>% summarise(MedLogLoss = median(ModelLogLoss, na.rm = TRUE),
                                               SDLogLoss = sd(ModelLogLoss, na.rm = TRUE))
  
