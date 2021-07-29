#' The RAT model (RAT +ve = COVID +ve and RAT -ve = COVID -ve) has very straight 
#' forward calculatation of ROCs as we do not need to work with an entire posterior
#' distribution.

# Source helper code and libraries
source("0000_HelperCode_Libraries/0001_Libraries.R")
source("0000_HelperCode_Libraries/0003_HelperFunctions.R")
RATonly <- readRDS("0100_Data/0103_nasal_dat.RDS") %>% select(nasal_ag, result, week)
RATonly <- RATonly %>% rename(CV = week)
# Create data frame to match format of other ROC dateframes.
classifications <- RATonly$nasal_ag - RATonly$result
classifications[classifications == 0 & RATonly$result == 1] <- 10 #"TruePos" 
classifications[classifications == 0 & RATonly$result == 0] <- -10 #"TrueNeg" 
RATonly <- cbind(RATonly, classifications)

RAT_only_ROC <- RATonly %>%
  group_by(CV) %>% 
  summarise(FalsePos = sum(classifications == 1), 
            TruePos  = sum(classifications == 10), 
            FalseNeg = sum(classifications == -1), 
            TrueNeg  = sum(classifications == -10)) %>% 
  ungroup() %>%
  mutate(FitType = "RATonly", SwabType = "nasal", threshold = NA)

RAT_only_ROC$FalseNegRate = RAT_only_ROC$FalseNeg/(RAT_only_ROC$FalseNeg + RAT_only_ROC$TruePos)
RAT_only_ROC$FalsePosRate = RAT_only_ROC$FalsePos/(RAT_only_ROC$FalsePos + RAT_only_ROC$TrueNeg)
RAT_only_ROC$TrueNegRate  = RAT_only_ROC$TrueNeg /(RAT_only_ROC$FalsePos + RAT_only_ROC$TrueNeg)
RAT_only_ROC$TruePosRate  = RAT_only_ROC$TruePos /(RAT_only_ROC$FalseNeg + RAT_only_ROC$TruePos)


RAT_only_ROC <- RAT_only_ROC  %>% 
  group_by(SwabType, FitType, threshold) %>%
  summarise(MedFalseNegRate = mean(FalseNegRate, na.rm = TRUE),
            SDFalseNegRate = NA,
            my_80CI(varnam = "FalseNegRate", x = FalseNegRate),
            MedFalsePosRate = mean(FalsePosRate, na.rm = TRUE),
            SDFalsePosRate = NA,
            my_80CI(varnam = "FalsePosRate",x = FalsePosRate),
            MedTrueNegRate = mean(TrueNegRate, na.rm = TRUE),
            SDTrueNegRate = NA,
            my_80CI(varnam = "TrueNegRate", x = TrueNegRate),
            MedTruePosRate = mean(TruePosRate, na.rm = TRUE),
            SDTruePosRate = NA,
            my_80CI(varnam = "TruePosRate", x = TruePosRate)) %>%
  pivot_wider(names_from = probs, values_from = ends_with("_CI"))



# Save output
saveRDS(RAT_only_ROC, "0400_ModelAssessment/0410_RATOnly_ROCrate.rds")
