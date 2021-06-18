nas_ROC <- readRDS("RATNeg_Workflow/0300_ModelSelection/nas_ROC_0.02.rds") %>%
  filter(FitType != "0__" & FitType != "RATonly")

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
  geom_line(#data = nas_ROC2 %>% filter(!FitType %in% c("0__", "RATonly")), 
    #aes(x = MedFalsePosRate, y = MedTruePosRate, colour = FitType)
  ) +
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
  geom_line(#data = nas_ROC2 %>% filter(!FitType %in% c("0__", "RATonly")), 
    #aes(x = MedFalseNegRate, y = MedTrueNegRate, colour = FitType)
  ) +
  geom_abline(slope = 1, alpha = 0.1) +
  ggtitle("Nasal Swabs") 


nas_mods <- nas_ROC2 
nas_mods$threshold <- parse_number(nas_ROC2$threshold)
nas_mods <-  nas_mods %>%
  filter(threshold > 0 & threshold < 100) 
acceptable_false_neg <- 0.2 
minimise_falseneg <- nas_mods %>% group_by(FitType) %>%
  slice(which.min((MedFalseNegRate - acceptable_false_neg) > 0))
acceptable_false_pos <- 0.2 
minimise_falsepos <- nas_mods %>% group_by(FitType) %>%
  slice(which.min((MedFalsePosRate - acceptable_false_pos) > 0))
max_perf <- nas_mods %>% group_by(FitType) %>%
  slice(which.max(sum(MedTruePosRate + MedTrueNegRate)))

scenario_outcomes <- rbind(max_perf, minimise_falseneg, minimise_falsepos)
scenario_outcomes$scenario <- rep(c("Maximise Correct Rates", "Minimise False Negatives", "Minimise False Positives"), each = 4)
scenario_outcomes_rounded <- scenario_outcomes%>% 
  mutate_if(is.numeric, ~round(.x, 2))

View(scenario_outcomes_rounded)
write_clip(scenario_outcomes_rounded)
?write.table
