source("0000_HelperCode_Libraries/0001_Libraries.R")
source("0000_HelperCode_Libraries/0003_HelperFunctions.R")

RAT_only_ROC <- data.frame("SwabType" = "nasal",
                           "FitType" = "RATonly",
                           "FalsePos" = sum(ifelse((RATresult - truth) == 1, 1, 0)),
                           "TruePos" = sum(ifelse(RATresult == 1 & truth == 1, 1, 0)),
                           "FalseNeg" = sum(ifelse((RATresult - truth) == -1, 1, 0)),
                           "TrueNeg" = sum(ifelse(RATresult == 0 & truth == 0, 1, 0)))

RAT_only_ROC$MedFalseNegRate = RAT_only_ROC$FalseNeg/(RAT_only_ROC$FalseNeg +
                                                     RAT_only_ROC$TruePos)
RAT_only_ROC$MedFalsePosRate = RAT_only_ROC$FalsePos/(RAT_only_ROC$FalsePos +
                                                     RAT_only_ROC$TrueNeg)
RAT_only_ROC$MedTrueNegRate =  RAT_only_ROC$TrueNeg/(RAT_only_ROC$FalsePos +
                                                    RAT_only_ROC$TrueNeg)
RAT_only_ROC$MedTruePosRate =  RAT_only_ROC$TruePos/(RAT_only_ROC$FalseNeg +
                                                    RAT_only_ROC$TruePos)

saveRDS(RAT_only_ROC, "0400_ModelAssessment/0410_RATOnly_ROCrate.rds")
