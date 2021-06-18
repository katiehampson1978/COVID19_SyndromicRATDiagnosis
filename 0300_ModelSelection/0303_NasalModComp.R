library(dtplyr)
library(tidyverse)
library(furrr)
rds_files <- fs::dir_ls("RATNeg_Workflow/0300_ModelSelection/Output/", regexp = "\\.rds$")

read_fun <- function(file){
  tmp <- readRDS(file)[["cv_log_loss"]]
  tmp$FitType <- as.character(tmp$FitType)
  tmp
}

final_models <- rds_files %>% 
  future_map_dfr(read_fun)

final_models1 <- final_models %>% 
  dplyr::filter((SwabType == "nasal")) %>% 
  dplyr::filter(!str_detect(FitType, "time_w"))

final_models_wide <- final_models1 %>% 
  pivot_wider(names_from = CV, values_from = log_loss, values_fn = list) %>%
  unnest(c(`1`, `2` , `3`, `4`,  `5`))

final_models_wide <- final_models_wide %>% 
  mutate(ModelLogLoss = final_models_wide %>%
           select(c(`1`, `2` , `3`,  `4`, `5`)) %>%
           rowSums())

nasal_log_loss <- final_models_wide %>% filter(SwabType == "nasal")

fittypes <- nasal_log_loss$FitType %>% unique()

round0s <- nasal_log_loss %>% filter(FitType %in%  fittypes[1:4])

ggplot(round0s, aes(x = ModelLogLoss, y = FitType)) +
  geom_boxplot()
bestround1 <- '0__'
round1s <- nasal_log_loss %>% filter(FitType %in%  fittypes[6:17])
bestround1 <- '1_fever_age'
round1s$accept <- ifelse(round1s$FitType == bestround1, 1, 0)
ggplot(round1s, aes(x = ModelLogLoss, y = FitType, fill = factor(accept))) +
  geom_boxplot()


round2s<- nasal_log_loss %>% filter(FitType %in%  fittypes[18:35])
bestround2 <- '2_fever.wet_cough_age'
round2s$accept <- ifelse(round2s$FitType == bestround2, 1, 0)
ggplot(round2s, aes(x = ModelLogLoss, y = FitType, fill = factor(accept))) +
  geom_boxplot()


round3s<- nasal_log_loss %>% filter(FitType %in%  fittypes[36:47])
bestround3 <- '3_fever.wet_cough.tired_age'
round3s$accept <- ifelse(round3s$FitType == bestround3, 1, 0)
ggplot(round3s, aes(x = ModelLogLoss, y = FitType, fill = factor(accept))) +
  geom_boxplot()

round4s<- nasal_log_loss %>% filter(FitType %in%  fittypes[48:50])
bestround4 <- '4_fever.wet_cough.tired.diarrhoea_age'
round4s$accept <- ifelse(round4s$FitType == bestround4, 1, 0)
ggplot(round4s, aes(x = ModelLogLoss, y = FitType, fill = factor(accept))) +
  geom_boxplot()

best_each_round <- c('0__', 
                     '1_fever_age', 
                     '2_fever.wet_cough_age', 
                     '3_fever.wet_cough.tired_age', 
                     '4_fever.wet_cough.tired.diarrhoea_age')

library(MLmetrics)
RAT_logloss <- LogLoss(nas_dat$nasal_ag, nas_dat$result)
 
nasal_log_loss_best <- nasal_log_loss %>% filter(FitType %in% best_each_round) 

modnam <- c("RAT + 4 Symptoms + Age", "RAT + 3 Symptoms + Age",
            "RAT + 2 Symptoms + Age", "RAT + 1 Symptoms + Age",
            "Intercept Only")

ggplot(nasal_log_loss_best, aes(x = ModelLogLoss, y = FitType)) +
  geom_boxplot() +
  ggtitle("Posterior Cross Entropy") +
  ylab("Model Type") +
  xlab("Cross Entropy") +
  scale_y_discrete(labels = modnam) +
  theme_bw()

# ROC Curves --------------------------------------------------------------

source("RATNeg_Workflow/0300_ModelSelection/0300_Functions.R")
library(tidyverse)
library(dtplyr)
library(furrr)
library(ROCR)
library(cowplot)
best_each_round_nas <- c('0__', '1_fever_age', 
                                         '2_fever.wet_cough_age', 
                                         '3_fever.wet_cough.tired_age', 
                                         '4_fever.wet_cough.tired.diarrhoea_age')

# rds_files <- fs::dir_ls("RATNeg_Workflow/0300_ModelSelection/Output", regexp = "\\.rds$")

# read_fun <- function(file){
#   readRDS(file)[["validation"]]
# }
# 
# final_models <- rds_files %>% 
#   future_map_dfr(read_fun)


nas0 <- readRDS("RATNeg_Workflow/0300_ModelSelection/Output/nasal_Round0__.rds")
nas1 <- readRDS("RATNeg_Workflow/0300_ModelSelection/Output/nasal_Round1_fever_age.rds")
nas2 <- readRDS("RATNeg_Workflow/0300_ModelSelection/Output/nasal_Round2_fever.wet_cough_age.rds")
nas3 <- readRDS("RATNeg_Workflow/0300_ModelSelection/Output/nasal_Round3_fever.wet_cough.tired_age.rds")
nas4 <- readRDS("RATNeg_Workflow/0300_ModelSelection/Output/nasal_Round4_fever.wet_cough.tired.diarrhoea_age.rds")

# nas_val <- final_models %>% 
#   dplyr::filter((SwabType == "nasal")) %>% 
#   filter(FitType %in%  best_each_round_nas)
nas_dat <- readRDS("RATNeg_Workflow/0100_DataPreparation/0100_nasal_dat.RDS")
nas_null <- data.frame("SwabType" = "nasal", "FitType" = "RATonly", "CV" = NA,
                       "Prob" = nas_dat$nasal_ag, "Truth" = nas_dat$result,
                       "Iter" = NA, "Chain" = NA)


nas_val <- rbind(nas0$validation,
                 nas1$validation,
                 nas2$validation,
                 nas3$validation,
                 nas4$validation,
                 nas_null)

nas_auc <- auc_calc(validation_df = nas_val, 
                    prob_range = c(0.00,1), 
                    increment_size = 0.05)

ggplot(nas_auc, aes(x = auc, fill = FitType)) +
  geom_histogram(position = "identity", alpha =  0.3) +
  ggtitle("Posterior AUC for Nasal Swabs")


nas_ROC <- ROC_diagnose(validation_df = nas_val, 
                        prob_range = c(0.00,1), 
                        increment_size = 0.05)

saveRDS(nas_ROC3, "RATNeg_Workflow/0300_ModelSelection/plotly/nas_ROC3.rds")

nas_ROC2 <- nas_ROC  %>% group_by(SwabType, FitType, threshold) %>%
  summarise(MedFalseNegRate = median(FalseNegRate, na.rm = TRUE),
            SDFalseNegRate = sd(FalseNegRate, na.rm = TRUE),
            MedFalsePosRate = median(FalsePosRate, na.rm = TRUE),
            SDFalsePosRate = sd(FalsePosRate, na.rm = TRUE),
            MedTrueNegRate = median(TrueNegRate, na.rm = TRUE),
            SDTrueNegRate = sd(TrueNegRate, na.rm = TRUE),
            MedTruePosRate = median(TruePosRate, na.rm = TRUE),
            SDTruePosRate = sd(TruePosRate, na.rm = TRUE))



ggplot(nas_ROC2, aes(x = MedFalsePosRate, y = MedTruePosRate, 
                     colour = FitType)) +
  geom_point() +
  geom_line(data = nas_ROC2 %>% filter(!FitType %in% c("0__", "RATonly")), 
            aes(x = MedFalsePosRate, y = MedTruePosRate, colour = FitType)) +
  geom_errorbar(aes(ymin = MedTruePosRate - SDTruePosRate, 
                    ymax = MedTruePosRate + SDTruePosRate)) +
  
  geom_errorbarh(aes(xmin = MedFalsePosRate - SDFalsePosRate, 
                     xmax = MedFalsePosRate + SDFalsePosRate)) +
  geom_line() +
  geom_abline(slope = 1, alpha = 0.1) +
  ggtitle("Nasal Swabs")

ggplot(nas_ROC2, aes(x = MedFalseNegRate, y = MedTrueNegRate, colour = FitType)) +
  geom_point() +
  geom_errorbar(aes(ymin = MedTrueNegRate - SDTrueNegRate, 
                    ymax = MedTrueNegRate + SDTrueNegRate)) +
  geom_errorbarh(aes(xmin = MedFalseNegRate - SDFalseNegRate, 
                     xmax = MedFalseNegRate + SDFalseNegRate)) +
  geom_line(data = nas_ROC2 %>% filter(!FitType %in% c("0__", "RATonly")), 
            aes(x = MedFalseNegRate, y = MedTrueNegRate, colour = FitType)) +
  geom_abline(slope = 1, alpha = 0.1) +
  ggtitle("Nasal Swabs") 


nas_ROC3 <- nas_ROC2
nas_ROC3$threshold <- parse_number(nas_ROC2$threshold)

acceptable_false_neg <- 0.125
minimise_falseneg <- nas_ROC3 %>% group_by(FitType) %>%
  slice(which.min(abs(MedFalseNegRate - acceptable_false_neg)))
model_names <- c("Intercept Only", "1 Symptom", "2 Symptom", "3 Symptom",
                 "4 Symptom", "RAT Only")
a <- ggplot(nas_ROC3, aes(x = threshold/100, y = MedTruePosRate, colour = FitType)) +
  geom_point() + 
  geom_vline(data = minimise_falseneg, alpha = 0.5,
             aes(xintercept = threshold/100,  colour = FitType)) +
  ggtitle("True Positive Rate") +
  ylab("") +
  xlab("") +
  coord_cartesian(xlim = c(0, 0.3))
legend <- get_legend(a)
a <- a + theme(legend.position="none")
a2 <- ggplot(minimise_falseneg, aes(x =FitType, 
                                    y = MedTruePosRate, 
                                    colour = FitType)) +
  geom_point() +
  geom_errorbar(aes(ymin = MedTruePosRate - SDTruePosRate, 
                    ymax = MedTruePosRate + SDTruePosRate)) +
  ylab("") +
  xlab("") +
  theme(legend.position="none") + 
  coord_flip() +
  scale_x_discrete(labels= model_names, position = "top")

b <- ggplot(nas_ROC3, aes(x = threshold/100, y = MedFalsePosRate, colour = FitType)) +
  geom_point() + 
  geom_vline(data = minimise_falseneg, alpha = 0.5,
             aes(xintercept = threshold/100,  colour = FitType)) +
  ggtitle("False Positive Rate") + 
  ylab("") +
  xlab("") +
  theme(legend.position="none") +
  coord_cartesian(xlim = c(0, 0.3))
b2 <- ggplot(minimise_falseneg, aes(x =FitType, 
                                    y = MedFalsePosRate, 
                                    colour = FitType)) +
  geom_point() +
  geom_errorbar(aes(ymin = MedFalsePosRate - SDFalsePosRate, 
                    ymax = MedFalsePosRate + SDFalsePosRate)) +
  ylab("") +
  xlab("") +
  theme(legend.position="none") + 
  coord_flip() +
  scale_x_discrete(labels= model_names, position = "top")

c <- ggplot(nas_ROC3, aes(x = threshold/100, y = MedTrueNegRate, colour = FitType)) +
  geom_point() + 
  geom_vline(data = minimise_falseneg, alpha = 0.5, 
             aes(xintercept = threshold/100, colour = FitType)) +
  ggtitle("True Negative Rate") + 
  ylab("") +
  xlab("") +
  theme(legend.position="none") +
  coord_cartesian(xlim = c(0, 0.3))
c2 <- ggplot(minimise_falseneg, aes(x =FitType, 
                                    y = MedTrueNegRate,
                                    colour = FitType)) +
  geom_point() +
  geom_errorbar(aes(ymin = MedTrueNegRate - SDTrueNegRate, 
                    ymax = MedTrueNegRate + SDTrueNegRate)) +
  ylab("") +
  xlab("") +
  theme(legend.position="none") + 
  coord_flip() +
  scale_x_discrete(labels= model_names, position = "top")

d <- ggplot(nas_ROC3, aes(x = threshold/100, y = MedFalseNegRate, colour = FitType)) +
  geom_point()+ 
  geom_vline(data = minimise_falseneg, alpha = 0.5, 
             aes(xintercept = threshold/100, colour = FitType)) +
  geom_hline(yintercept = acceptable_false_neg, alpha = 0.5) +
  ggtitle("False Negative Rate") + 
  xlab("Threshold") +
  ylab("") +
  theme(legend.position="none") +
  coord_cartesian(xlim = c(0, 0.3))

d2 <- ggplot(minimise_falseneg, aes(x =FitType, 
                                    y = MedFalseNegRate, 
                                    colour = FitType)) +
  ylab("") +
  xlab("Model") +
  geom_point() +
  geom_errorbar(aes(ymin = MedFalseNegRate - SDFalseNegRate, 
                    ymax = MedFalseNegRate + SDFalseNegRate)) +
  theme(legend.position="none") + 
  coord_flip() +
  scale_x_discrete(labels= model_names, position = "top")




plot_grid(a, a2, b, b2, c,c2, d, d2, 
          ncol = 2, 
          # labels = c("A", "", "B", "", "C", "",  "D", ""),
          # label_size = 10,
          # rel_widths = c(2, 1)
          align = "h"
          )

# ggplot(nas_ROC, aes(x = FalsePosRate, y = TruePosRate, 
#                      colour = FitType)) +
#   geom_point(alpha = 0.01) +
#   # geom_line(nas_ROC2, aes(x = MedFalsePosRate, y = MedTruePosRate,
#                           # colour = FitType)) +
#   geom_abline(slope = 1, alpha = 0.5) +
#   facet_grid(vars(FitType)) +
#   ggtitle("Nasal Swabs")
# 
# ggplot(nas_ROC, aes(x = FalseNegRate, y = TrueNegRate, colour = FitType)) +
#   geom_point(alpha = 0.1) +
#   geom_abline(slope = 1, alpha = 0.1) +
#   facet_grid(vars(FitType)) +
#   ggtitle("Nasal Swabs") + 
#   theme(legend.position="none")


