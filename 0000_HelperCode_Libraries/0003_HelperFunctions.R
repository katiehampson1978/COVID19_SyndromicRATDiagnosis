#' A collection of helper functions to declutter code elsewhere. These functions
#' were not written as part of a package so their use outwith the confines of
#' this analysis is not recommended!


# Raw Data Preparation ----------------------------------------------------

#' Prepare data for cross validation and model fitting in CmdStanR This code takes
#' the raw data, a table defining cross validation sets, and the names of response
#' and explanatory variables, and returns a list of lists ready for input to 
#' Stan fitting code.
stan_dat_prep <- 
  function(dat, # Raw data
           swab_type, # Name of swab type
           cross_val_table, # Assignment table of cross validation sets
           covariate_names, # Names of explanatory variables
           response_names # Names of response variables
           ){
    # Convert swab name to match raw data
      swab_name <- "nasal_ag"
    # Nest data to give one row per week and bind with cross validation definition
    nested_dat <- dat %>%
      nest(data = c(-week)) %>%
      cbind(cross_val_table)
    # Initialise storage
    output_list <- list()
    # Create test and training objects
    for(i in seq_along(1:ncol(cross_val_table))){
      ## Select current cross validation set
      var <- paste0("cv_", i)
      # Identify object name for test data
      test_name <- paste0(swab_type,  paste0("_cv_test_", i))
      # Identify object name for training data
      train_name <- paste0(swab_type,  paste0("_cv_train_", i))
      # Create list to pass to Stan code
      dat_name <- paste0(swab_type,  paste0("_stan_dat_", i))
      # Generate test data object based on name above
      assign(
        paste(test_name), 
        nested_dat %>% 
          filter(get(var) == 1) %>%
          unnest(cols = data)
      )
      # Generate training data object based on name above
      assign(
        paste(train_name), 
        nested_dat %>% 
          filter(get(var) != 1) %>%
          unnest(cols = data)
      )
      # Generate Stan data list based on name above
      assign(paste(dat_name),
             # If there are no covariates, the intercept requires an extra
             # argument to keep the matrix format
               if (is.null(covariate_names)) {
                 list(
                   # Size of training set
                   N_train = nrow(get(paste(train_name))),
                   # Size of testing set
                   N_test = nrow(get(paste(test_name))),
                   # Response data for training set
                   y_train = get(paste(train_name)) %>%
                     dplyr::select(starts_with(response_names)),
                   # Response data for testing set
                   y_test = get(paste(test_name)) %>%
                     dplyr::select(starts_with(response_names)),
                   # COVID result for test set
                   y_test_ll = unlist(get(paste(test_name)) %>%
                                        dplyr::select(starts_with(
                                          response_names[1]
                                        )),
                                      use.names = FALSE),
                   # Intercept only for training set
                   x_train = matrix(1, nrow = nrow(get(paste(train_name))), 
                                    ncol = 1),
                   # Intercept only for testing set
                   x_test = matrix(1, nrow = nrow(get(paste(test_name))), 
                                   ncol = 1),
                   # Rapid antigen test result for training set
                   RAT_train = get(paste(train_name)) %>%
                     dplyr::select(starts_with(swab_name)) %>%
                     unlist(use.names = FALSE),
                   # Rapid antigen test result for testing set
                   RAT_test = get(paste(test_name)) %>%
                     dplyr::select(starts_with(swab_name)) %>%
                     unlist(use.names = FALSE),
                   # Number of covariates (including intercept)
                   K = length(covariate_names) + 1,
                   # Number of response dimensions (including COVID result)
                   D = length(response_names),
                   # Cross validation fold ID
                   fold = i
                 )
               } else
                 list(
                   # Size of training set
                   N_train = nrow(get(paste(train_name))),
                   # Size of testing set
                   N_test = nrow(get(paste(test_name))),
                   # Response data for training set
                   y_train = get(paste(train_name)) %>%
                     dplyr::select(starts_with(response_names)),
                   # Response data for testing set
                   y_test = get(paste(test_name)) %>%
                     dplyr::select(starts_with(response_names)),
                   # COVID result for test set
                   y_test_ll = unlist(get(paste(test_name)) %>%
                                        dplyr::select(starts_with(
                                          response_names[1]
                                        )),
                                      use.names = FALSE),
                   # Intercept only for training set
                   x_train = cbind(1, get(paste(train_name)) %>%
                                     dplyr::select(starts_with(
                                       covariate_names
                                     ))),
                   # Intercept only for testing set
                   x_test = cbind(1, get(paste(test_name)) %>%
                                    dplyr::select(starts_with(
                                      covariate_names
                                    ))),
                   # Rapid antigen test result for training set
                   RAT_train = get(paste(train_name)) %>%
                     dplyr::select(starts_with(swab_name)) %>%
                     unlist(use.names = FALSE),
                   # Rapid antigen test result for testing set
                   RAT_test = get(paste(test_name)) %>%
                     dplyr::select(starts_with(swab_name)) %>%
                     unlist(use.names = FALSE),
                   # Number of covariates (including intercept)
                   K = length(covariate_names) + 1,
                   # Number of response dimensions (including COVID result)
                   D = length(response_names),
                   # Cross validation fold ID
                   fold = i
                 )
      )
      
      output_list[[i]] <- get(paste(dat_name))
    }
    output_list
  }


# Model Fitting -----------------------------------------------------------

#' Run multivariate probit code in Stan for a given output of *stan_dat_prep*. 
#' Function uses the furrr package (wrapper for purrr and the future package) to
#' run different cross validation sets in parallel.
run_MV_cv_model_parr <- function(standat, # stan ready data (output from stan_dat_prep)
                                 swab_type, # name of rapid antigen test swab
                                 fit_type, # name of model fit
                                 n_cv_folds, # number of cross validation folds
                                 syndromic_only = TRUE # Is model syndromic only or does it include RAT
                                 ){
  # Load the global fitting params
  source("0000_HelperCode_Libraries/0002_GlobalStanFittingParams.R")
  # If syndromic only model we want the standard MV Probit, if RAT positive we
  # want the two step model formulation
  if(syndromic_only){
    model_file <- "0200_ModelDefinitions/0201_ModelDef_MVProbit.stan"
  }else{
    model_file <- "0200_ModelDefinitions/0202_ModelDef_MVProbit_Split.stan"
  }
  # Compile model
  model <- cmdstan_model(model_file) 
  # Default to future multisession plan (system agnostic)
  plan(multisession)
  # In parallel, fit model in cmdstan to each cross validation set of data
  fit <- future_map(1:n_cv_folds, ~model$sample(standat[.x][[1]],
                          chains = N_Chains,
                          parallel_chains = N_Chains,
                          refresh = 100,
                          # iter_warmup=Iter,
                          # iter_sampling = Iter,
                          ## Note: using params from the "fitting params" file here
                          # restrict the initial range to limit numerical overflow 
                          # (inits being rejected)
                          init = Init_R,
                          # Only effects performance if reached
                          max_treedepth = TreeDepth,
                          # Think of like step length in Gibbs, this will slow things down
                          adapt_delta = AdaptDelta,
                          save_warmup=TRUE),
                    .options = furrr_options(seed = NULL))
  # Return list of model fits
 fit
}

# Model Tidying ----------------------------------------------

#' Take ouput of *run_MV_cv_model_parr* and return list of tidy data frames for 
#' easier downstream analyses
tidy_output_MV <- function(output_list, # list output from *run_MV_cv_model_parr*
                           swab_type, # name of swab type
                           fit_type  = "MV", # name of fit type
                           pred_names, # names of covariates
                           response_names, # names of response dimensions
                           n_cv_folds # number of cross validation folds
                           ){
  ## Initialise data frames
  # This will store the log loss (aka cross entropy)
  cv_log_loss <- data.frame("CV" = c(), "log_loss" = c(), 
                            "Iter" = c(), "Chain" = c())
  # This will store the beta coefficients
  # Generate covariate-response dimension combinations
  resp_cov_combos <- expand.grid(c(response_names),
                                 c("Intercept", pred_names), 
                                 stringsAsFactors = FALSE, 
                                 KEEP.OUT.ATTRS = TRUE) %>% 
    unite(col = "temp") %>%
    unlist(use.names = FALSE)
  # Generate dataframe of correct dimensions to overwrite
  betas <- matrix(NA, 
                  ncol = 3 + #+ Cross Validation +iter +chain
                    length(resp_cov_combos)) %>% 
    as.data.frame()
  # Fix variable names
  names(betas) <- c("CV", resp_cov_combos, "Iter", "Chain")
  # this will store the correlation coefficients
  omegas <- matrix(NA, 
                   ncol = 3 + #+ Cross Validation
                     length(response_names)^2) %>% 
    as.data.frame()
  # Fix variable names
  names(omegas) <- c("CV", expand.grid(response_names, response_names) %>% 
                           unite("Omega_nam", Var1:Var2, sep="_") %>% 
                       unlist(use.names = FALSE), "Iter", "Chain")
  # This will store the predicted probability and true values
  validation <- data.frame("CV" = c(), "Prob" = c(),
                           "Truth" = c(),
                           "Pat" = c(),
                           "Iter" = c(),
                           "Chain" = c()
  )
  
  # Iterate through the model output for each fold
  for(i in seq_along(1:n_cv_folds)){
    # Create object name
    obj_name <- paste0(swab_type, i)
    # Generate a unique data frame from the list object
    assign(paste(obj_name), enframe(output_list) %>% 
             filter(name == i) %>%
             unnest(c(value)))
    # Extract log loss and bind to other log losses
    CV <- get(obj_name) %>% select(cv_fold) %>% unlist(use.names = FALSE)
    temp_log_loss  <- data.frame("CV" = CV, 
                                 "log_loss" = get(obj_name) %>% 
                                   select(log_loss) %>% 
                                   unlist(use.names = FALSE),
                                "Iter"= get(obj_name) %>% 
                                   select(.iteration) %>% 
                                   unlist(use.names = FALSE),
                                "Chain"= get(obj_name) %>% 
                                  select(.chain) %>% 
                                  unlist(use.names = FALSE))
    cv_log_loss <- rbind(cv_log_loss, temp_log_loss)
    ## Extract beta estimates and bind to other betas
    # Extract results for this round
    temp_beta <- data.frame("CV" = CV, get(obj_name) %>% 
                              select(starts_with("beta")),
                                     "Iter"= get(obj_name) %>% 
                                       select(c(".iteration")),
                                     "Chain"= get(obj_name) %>% 
                                       select(c(".chain")))
    # Fix variable names
    names(temp_beta) <- c("CV", resp_cov_combos, "Iter", "Chain")
    # Bind to previous results
    betas <- rbind(betas, temp_beta)
    
    # Extract correlation estimates
    if(length(response_names) == 1){
      corr_names <- "result_result"
    }else{
      corr_names <- expand.grid(response_names, response_names) %>% 
        unite("Omega_nam", Var1:Var2, sep="_") %>% 
        unlist(use.names = FALSE)
    }
    # Extract results for this round
    temp_omegas <- data.frame("CV" = CV, get(obj_name) %>% 
                                select(starts_with("Omega")),
                              "Iter"= get(obj_name) %>% 
                                select(c(".iteration")),
                              "Chain"= get(obj_name) %>% 
                                select(c(".chain")))
    # Fix variable names
    names(temp_omegas) <- c("CV", corr_names,
                            "Iter",
                            "Chain"
                            )
    # Bind to previous results
    omegas <- rbind(omegas, temp_omegas)
    ## Extract probabilities and true values
    probs <- get(obj_name) %>% 
      select(cv_fold, starts_with(c("logprob", ".iteration", ".chain"))) %>% 
      pivot_longer(cols = starts_with("logprob"), values_to = "logprob")
    # Convert from log probability to probability
    probs$Prob <- exp(probs$logprob)
    fold_size <- get(obj_name) %>% select(starts_with("y_true")) %>% ncol()
    truth <- get(obj_name) %>%
      select(starts_with("y_true")) %>%
      pivot_longer(cols = starts_with("y_true"), values_to = "Truth")
    validation_temp <- data.frame("CV" = probs$cv_fold, 
                                  "Prob" = probs$Prob,
                                  "Truth" = truth$Truth,
                                  "Iter" = probs$.iteration,
                                  "Chain" = probs$.chain)
    validation <- rbind(validation, validation_temp)
  }
  # Drop NA row from betas and omegas
  betas <- betas[-1,] 
  omegas <- omegas[-1,] 
  # Add swab type and fit type metadata
  omegas <- cbind("SwabType" = swab_type, "FitType" = fit_type, omegas)
  betas <- cbind("SwabType" = swab_type, "FitType" = fit_type, betas)
  validation <- cbind("SwabType" = swab_type, "FitType" = fit_type, validation)
  cv_log_loss <- cbind("SwabType" = swab_type, "FitType" = fit_type, cv_log_loss)
  # combine into list for export
  results_list <- list("betas" = betas, 
                       "validation" = validation, 
                       "cv_log_loss" = cv_log_loss,
                       "omegas" = omegas)
  results_list
}


# One stop model shop -----------------------------------------------------

#' This function combines the previous functions into a single, overarching
#' statement that takes raw data, prepares it for Stan, fits the model and
#' returns a list of tidy data frames,

## Run cv model
tidy_run_cv <- function(dat, swab_type,  covariate_names, 
                        symptom_nam, model_type, cross_val_table, 
                        syndromic_only = TRUE){
  # Wrap data for stan
  standat <- stan_dat_prep(dat = dat, 
                           cross_val_table = cross_val_table, 
                           swab_type = swab_type, 
                           covariate_names = covariate_names, 
                           response_names = c("result", symptom_nam)
                           )
  # Run model for each CV fold
  fit <- run_MV_cv_model_parr(standat = standat, swab_type = swab_type,
                              fit_type = model_type,
                              n_cv_folds = 5,
                              syndromic_only = syndromic_only)
  # Extract draws data frame
  temp_fun <- function(x, fit){
    fit_temp <- fit[[x]]
    draws <- fit_temp$draws(c("beta", "Omega",
                              "logprob", "cv_fold",
                              "log_loss",
                              "y_true"))
    draws_df <- as_draws_df(draws)
    draws_df
  }
  # Extract draws for each CV fold
  output_list <- map(1:5, ~temp_fun(.x, fit = fit))
  # Convert draws to data frame
  output_df <- tidy_output_MV(output_list = output_list,
                              swab_type = swab_type,
                              fit_type = model_type,
                              pred_names = covariate_names,
                              response_names = c("result", symptom_nam),
                              n_cv_folds = 5)
  output_df
}


# Measuring Classifier Performance ----------------------------------------

#' Threshold takes a series of probabilities (such as those provided by the
#' validation dataframe output of *tidy_run_cv*) and converts them to classifications
#' based on a given range and resolution of threshold values
threshold <- function(validation_df, prob_range, increment_size){
  # Generate vector of probability values
  thresholds <- seq(from = prob_range[1], to = prob_range[2], by = increment_size)
  # For each probability determine whether classifier will return 1 or 0
  classifier_df <- map_dfc(thresholds, ~if_else(validation_df$Prob >= .x, 1, 0))
  # Make names nicer
  names(classifier_df) <- paste0("threshold", thresholds*100)
  # Combine with original dataframe
  output_df <- cbind(validation_df, classifier_df)
  output_df
}

#' This function generates classifications based on a given range and resolution 
#' of threshold values and calculates the classification error rates for these
#' thresholds (i.e. false and true positives and negatives) 
ROC_diagnose <- function(validation_df, prob_range, increment_size){
  classification_df <- threshold(validation_df = validation_df, 
                                 prob_range = c(prob_range[1], prob_range[2]),
                                 increment_size = increment_size)
  # Subtract true value from each classified value
  classification_df <- classification_df %>% 
    mutate(across(starts_with("threshold"), ~ .x - Truth))
  # Extract classified values and convert to false & true positives and false & true negatives
  classifications <- classification_df %>% select(starts_with("threshold"))
  # Encode as +- 10 as it takes up less memory than character vectors -
  # could do this more neatly for sure!
  classifications[classifications == 0 & classification_df$Truth == 1] <- 10 #"TruePos" 
  classifications[classifications == 0 & classification_df$Truth == 0] <- -10 #"TrueNeg" 
  # Combinewith meta data
  diagnosed_df <- cbind(classification_df %>% 
                          select(SwabType, FitType, CV, Iter, Chain), 
                        classifications)
  # Have one column for threshold and one column for value at threshold rather
  # than one column per threshold.
  diagnosed_grouped_df <- diagnosed_df %>% 
    pivot_longer(cols = starts_with("threshold"), 
                 names_to = "threshold", 
                 values_to = "classification") 
  # Get classification counts
  diag_df <- lazy_dt(diagnosed_grouped_df) %>%
    group_by(SwabType, FitType, Iter, Chain, threshold) %>%
    summarise(FalsePos = sum(classification == 1), #"FalsePos"),
              TruePos = sum(classification == 10), #"TruePos"),
              FalseNeg = sum(classification == -1), #"FalseNeg"),
              TrueNeg = sum(classification == -10)) %>% #"TrueNeg")) %>%
    ungroup() %>%
    as.data.frame()
  # Get classification rates
  diag_df$FalseNegRate = diag_df$FalseNeg/(diag_df$FalseNeg + diag_df$TruePos)
  diag_df$FalsePosRate = diag_df$FalsePos/(diag_df$FalsePos + diag_df$TrueNeg)
  diag_df$TrueNegRate = diag_df$TrueNeg/(diag_df$FalsePos + diag_df$TrueNeg)
  diag_df$TruePosRate = diag_df$TruePos/(diag_df$FalseNeg + diag_df$TruePos)

  diag_df
}




# Function to read in file and extract validation dataframe
read_validation <- function(file){
  tmp <- readRDS(file)[["validation"]]
  tmp
}
# Function to read in file and extract log loss dataframe
read_logloss <- function(file){
  tmp <- readRDS(file)[["cv_log_loss"]]
  tmp
}
# Extract 80% credible interval
my_80CI <- function(varnam, x) {
  probs <- c(0.1,0.9)
  tibble("{varnam}_CI" := quantile(x, probs, na.rm = TRUE, names = FALSE), probs = probs)
}


