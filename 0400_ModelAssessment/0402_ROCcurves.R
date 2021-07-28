#' Working with our ROCs we can now generate ROC curves and identify the probability
#' thresholds at which key criteria are met (i.e. our scenarios)
# Load libraries and helper functions
source("0000_HelperCode_Libraries/0001_Libraries.R")
source("0000_HelperCode_Libraries/0003_HelperFunctions.R")
# LoadROC outputs
SyndOnlyROC <- readRDS("0400_ModelAssessment/0410_SyndromicOnly_ROCrate.rds")
SyndRATROC <- readRDS("0400_ModelAssessment/0410_SyndromicRAT_ROCrate.rds")
RATOnlyROC <- readRDS("0400_ModelAssessment/0410_RATOnly_ROCrate.rds")

# Add model class column to help keep track
SyndOnlyROC$ModelClass <- "SyndOnly"
SyndRATROC$ModelClass <- "SyndRAT"
RATOnlyROC$ModelClass <- "RATonly"

# Combine into one data frame
ROCs <- rbind(SyndOnlyROC,
              SyndRATROC,
              RATOnlyROC)

# Modify fit type name to include class, symptom number and age 
ROCs$FitType[1:500] <- paste0(ROCs$ModelClass[1:500], "_", 
                              parse_number(ROCs$FitType[1:500]), 
                              "Symptom")
# Change RATonly row
ROCs$FitType[501] <- "RATonly"

# Create mirrored data frame to make plotting easier
ROC_plots <- ROCs
ROC_plots$FitType <-  str_replace(ROC_plots$FitType, ".*_", "")

saveRDS(ROC_plots, "0400_ModelAssessment/0420_ROC_plot_dat.rds")
# True Positive vs False Positive plot
ggplot(ROC_plots , 
       aes(x = MedFalsePosRate, y = MedTruePosRate, 
                     colour = FitType)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = MedTruePosRate - SDTruePosRate, 
                    ymax = MedTruePosRate + SDTruePosRate)) +
  
  geom_errorbarh(aes(xmin = MedFalsePosRate - SDFalsePosRate, 
                     xmax = MedFalsePosRate + SDFalsePosRate)) +
  geom_line() +
  geom_abline(slope = 1, alpha = 0.1) +
  ylab("True Positive Rate") +
  xlab("False Positive Rate") +
  facet_wrap(~ModelClass) +
  ggtitle("ROC Curve For Each Model Class") +
  scale_colour_discrete(name = "Model Class", 
                        labels = c("0 Symptoms", 
                                   "1 Symptom ", 
                                   "2 Symptoms", 
                                   "3 Symptoms", 
                                   "4 Symptoms", 
                                   "RAT Only"))

# True Negative vs False Negative plot
ggplot(ROC_plots, 
       aes(x = MedFalseNegRate, y = MedTrueNegRate, 
           colour = FitType)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = MedTrueNegRate - SDTrueNegRate, 
                    ymax = MedTrueNegRate + SDTrueNegRate)) +
  
  geom_errorbarh(aes(xmin = MedFalseNegRate - SDFalseNegRate, 
                     xmax = MedFalseNegRate + SDFalseNegRate)) +
  geom_line() +
  geom_abline(slope = 1, alpha = 0.1) +
  facet_wrap(~ModelClass)

# Change threshold variable to numeric
ROCs$threshold <- parse_number(ROCs$threshold)

# Costly false negatives scenario
acceptable_false_neg <- 0.2 # Acceptable Rate
# Identify threshold for each fit type that most closely exceed criterion 
minimise_falseneg <- ROCs %>% group_by(FitType) %>%
  mutate(tmp = MedFalseNegRate - acceptable_false_neg) %>%
  filter(tmp < 0) %>%
  slice(which.max(tmp)) %>%
  select(-tmp) 

# Costly false positives scenario
acceptable_false_pos <- 0.2 # Acceptable Rate
# Identify threshold for each fit type that most closely exceed criterion
minimise_falsepos <- ROCs %>% group_by(FitType) %>%
  mutate(tmp = MedFalsePosRate - acceptable_false_pos) %>%
  filter(tmp < 0) %>%
  slice(which.max(tmp)) %>%
  select(-tmp) %>%
  rbind(RATOnlyROC) # Need to add in RATOnly as false pos will go negative here

# Agnostic scenario
# Identify threshold for each fit type that most closely exceed criterion
max_perf <- ROCs %>% group_by(FitType) %>%
  mutate(tmp = (MedTruePosRate + MedTrueNegRate)) %>%
  filter(tmp > 0) %>%
  slice(which.max(tmp)) %>%
  select(-tmp) 

# Add scenario labels
max_perf$scenario <- "Agnostic"
minimise_falseneg$scenario <- "Costly False Negatives"
minimise_falsepos$scenario <- "Costly False Positives"

# Calculate median error for each scenario
max_perf$MedError <- (max_perf$MedFalseNegRate + max_perf$MedFalsePosRate)
minimise_falseneg$MedError <- minimise_falseneg$MedFalsePosRate
minimise_falsepos$MedError <- minimise_falsepos$MedFalseNegRate 

# Calculate SD of error for each scenario
max_perf$SDError <- (max_perf$SDFalseNegRate + max_perf$SDFalsePosRate)
minimise_falseneg$SDError <- minimise_falseneg$SDFalsePosRate
minimise_falsepos$SDError <- minimise_falsepos$SDFalseNegRate 

# Remove RATonly error as it is technically infinite
# minimise_falseneg[minimise_falseneg$FitType == "RATonly", "MedError"] <- NA

# Combine best thresholds for each scenario into one dataframe
scenario_outcomes <- rbind(max_perf, minimise_falseneg, minimise_falsepos)

# Round values
scenario_outcomes_rounded <- scenario_outcomes %>% 
  mutate_if(is.numeric, ~round(.x, 2))

# Tidy names
scenario_outcomes_rounded$FitType <- str_replace(scenario_outcomes_rounded$FitType, ".*_", "")
View(scenario_outcomes_rounded)
saveRDS(scenario_outcomes_rounded, "0400_ModelAssessment/0420_scenario_outcomes.rds")
# Plot
ggplot(scenario_outcomes_rounded, 
       aes(x = FitType, y = MedError, colour = ModelClass)) +
  geom_point() +
  geom_errorbar(aes(ymin = MedError - SDError, 
                    ymax = MedError + SDError)) +

  ylab("Error") + 
  xlab("Model") +
  facet_wrap(~scenario)  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_colour_discrete(name = "Model Class", 
                      labels = c("RAT only", "Syndromic Only", "Syndromic + RAT"))

