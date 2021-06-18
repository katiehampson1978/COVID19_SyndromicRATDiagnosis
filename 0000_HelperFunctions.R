# Functions to take data and cross validation structure and return a list of stan-
# ready input lists.
stan_dat_prep <- 
  function(dat, cross_val_table, swab_type, covariate_names, response_names,
           model_type = "MV_ModSel"){
    if(swab_type == "nasal"){
      swab_name <- "nasal_ag"
    } else {
      swab_name = "sal_ag"
    }
    nested_dat <- dat %>%
      nest(data = c(-week)) %>%
      cbind(cross_val_table)
    output_list <- list()
    # Create test and training objects
    for(i in seq_along(1:ncol(cross_val_table))){
      var <- paste0("cv_", i)
      test_name <- paste0(swab_type,  paste0("_cv_test_", i))
      train_name <- paste0(swab_type,  paste0("_cv_train_", i))
      dat_name <- paste0(swab_type,  paste0("_stan_dat_", i))
      assign(
        paste(test_name), 
        nested_dat %>% 
          filter(get(var) == 1) %>%
          unnest(cols = data)
      )
      assign(
        paste(train_name), 
        nested_dat %>% 
          filter(get(var) != 1) %>%
          unnest(cols = data)
      )
      
      assign(paste(dat_name),
               if (is.null(covariate_names)) {
                 list(
                   N_train = nrow(get(paste(train_name))),
                   N_test = nrow(get(paste(test_name))),
                   y_train = get(paste(train_name)) %>%
                     dplyr::select(starts_with(response_names)),
                   y_test = get(paste(test_name)) %>%
                     dplyr::select(starts_with(response_names)),
                   y_test_ll = unlist(get(paste(test_name)) %>%
                                        dplyr::select(starts_with(
                                          response_names[1]
                                        )),
                                      use.names = FALSE),
                   x_train = matrix(1, nrow = nrow(get(paste(train_name))), ncol = 1),
                   x_test = matrix(1, nrow = nrow(get(paste(test_name))), ncol = 1),
                   RAT_train = get(paste(train_name)) %>%
                     dplyr::select(starts_with(swab_name)) %>%
                     unlist(use.names = FALSE),
                   RAT_test = get(paste(test_name)) %>%
                     dplyr::select(starts_with(swab_name)) %>%
                     unlist(use.names = FALSE),
                   K = length(covariate_names) + 1,
                   D = length(response_names),
                   fold = i
                 )
               } else
                 list(
                   N_train = nrow(get(paste(train_name))),
                   N_test = nrow(get(paste(test_name))),
                   y_train = get(paste(train_name)) %>%
                     dplyr::select(starts_with(response_names)),
                   y_test = get(paste(test_name)) %>%
                     dplyr::select(starts_with(response_names)),
                   y_test_ll = unlist(get(paste(test_name)) %>%
                                        dplyr::select(starts_with(
                                          response_names[1]
                                        )),
                                      use.names = FALSE),
                   x_train = cbind(1, get(paste(train_name)) %>%
                                     dplyr::select(starts_with(
                                       covariate_names
                                     ))),
                   x_test = cbind(1, get(paste(test_name)) %>%
                                    dplyr::select(starts_with(
                                      covariate_names
                                    ))),
                   RAT_train = get(paste(train_name)) %>%
                     dplyr::select(starts_with(swab_name)) %>%
                     unlist(use.names = FALSE),
                   RAT_test = get(paste(test_name)) %>%
                     dplyr::select(starts_with(swab_name)) %>%
                     unlist(use.names = FALSE),
                   K = length(covariate_names) + 1,
                   D = length(response_names),
                   fold = i
                 )
      )
      
      output_list[[i]] <- get(paste(dat_name))
    }
    output_list
  }

source("RATNeg_Workflow/0200_Models/0202_FittingParams.R")
# Code to run a univariate model 
run_MV_cv_model <- function(standat, n_cv_folds, swab_type, fit_type, modsel = FALSE){
  model_file <- "RATNeg_Workflow/0200_Models/0203_ModelDef_MVProbit.stan"
  model <- cmdstan_model(model_file) 
  output_list <- list()
  for(i in seq_along(1:n_cv_folds)){
    fit <- model$sample(standat[i][[1]],
                        chains = if(modsel){2}else {4},
                        parallel_chains = 4,
                        refresh = 100,
                        ## Note: using params from the "fitting params" file here
                        # restrict the initial range to limit numerical overflow 
                        # (inits being rejected)
                        init = Init_R,
                        # Only effects performance if reached
                        max_treedepth = TreeDepth,
                        # Think of like step length in Gibbs, this will slow things down
                        adapt_delta = AdaptDelta,
                        save_warmup=TRUE)
    temp <- paste0(paste(swab_type, fit_type, sep = "_"), i, sep = "_")
    # Extract samples for key params 
    assign(paste(temp),
           as_draws_df(fit$draws(variable = 
                                   c("beta", "Omega", 
                                     "logprob", "cv_fold",
                                     "log_loss",  
                                     "y_true"))))
    output_list[[i]] <- get(paste(temp))
  }
  output_list  
}

run_MV_cv_model_parr <- function(standat, swab_type, fit_type, modsel = FALSE, n_cv_folds){
  model_file <- "RATNeg_Workflow/0200_Models/0203_ModelDef_MVProbit.stan"
  model <- cmdstan_model(model_file) 
  plan(multisession)
  fit <- future_map(1:n_cv_folds, ~model$sample(standat[.x][[1]],
                          chains = if(modsel){2}else {4},
                          parallel_chains = 4,
                          refresh = 100,
                          # iter_warmup=10,
                          # iter_sampling = 10,
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
 fit
}

# Multivariate Model Tidying ----------------------------------------------



tidy_output_MV <- function(output_list, swab_type, fit_type  = "MV", pred_names, 
                           response_names, n_cv_folds){
  ## Initialise data frames
  # This will store the log loss (aka cross entropy)
  cv_log_loss <- data.frame("CV" = c(), "log_loss" = c(), "Iter" = c(), "Chain" = c())
  # This will store the beta coefficients
  # Generate covariate-response dimension combinations
  resp_cov_combos <- expand.grid(c(response_names),
                                 c("Intercept", pred_names), 
                                 stringsAsFactors = FALSE, 
                                 KEEP.OUT.ATTRS = TRUE) %>% 
    unite(col = "temp") %>%
    unlist(use.names = FALSE)
  
  betas <- matrix(NA, 
                  ncol = 3 + #+ Cross Validation +iter +chain
                    length(resp_cov_combos)) %>% 
    as.data.frame()
  
  names(betas) <- c("CV", resp_cov_combos, "Iter", "Chain")
  # this will store the correlation coefficients
  omegas <- matrix(NA, 
                   ncol = 3 + #+ Cross Validation
                     length(response_names)^2) %>% 
    as.data.frame()
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
    # Extract beta estimates and bind to other betas
    temp_beta <- data.frame("CV" = CV, get(obj_name) %>% 
                              select(starts_with("beta")),
                                     "Iter"= get(obj_name) %>% 
                                       select(c(".iteration")),# %>% 
                                       #unlist(use.names = FALSE),
                                     "Chain"= get(obj_name) %>% 
                                       select(c(".chain")))
                            #%>% 
                                       #unlist(use.names = FALSE)))
    names(temp_beta) <- c("CV", resp_cov_combos, "Iter", "Chain")
    betas <- rbind(betas, temp_beta)
    #
    if(length(response_names) == 1){
      corr_names <- "result_result"
    }else{
      corr_names <- expand.grid(response_names, response_names) %>% 
        unite("Omega_nam", Var1:Var2, sep="_") %>% 
        unlist(use.names = FALSE)
    }
    
    temp_omegas <- data.frame("CV" = CV, get(obj_name) %>% select(starts_with("Omega")),
                              "Iter"= get(obj_name) %>% 
                                select(c(".iteration")),# %>% 
                              #unlist(use.names = FALSE),
                              "Chain"= get(obj_name) %>% 
                                select(c(".chain")))
    names(temp_omegas) <- c("CV", corr_names,
                            "Iter",
                            "Chain"
                            )
    omegas <- rbind(omegas, temp_omegas)
    # Extract probabilities and true values
    probs <- get(obj_name) %>% 
      select(cv_fold, starts_with(c("logprob", ".iteration", ".chain"))) %>% 
      pivot_longer(cols = starts_with("logprob"), values_to = "logprob")
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
  # Drop NA row from betas
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

## Run cv modelå
tidy_run_cv <- function(dat, swab_type,  covariate_names, 
                        symptom_nam, model_type, cross_val_table){
  standat <- stan_dat_prep(dat = dat, cross_val_table = cross_val_table, 
                           swab_type = swab_type, 
                           covariate_names = covariate_names, 
                           response_names = c("result", symptom_nam), 
                           model_type = model_type)
  fit <- run_MV_cv_model_parr(standat = standat, swab_type = swab_type, 
                              fit_type = model_type, modsel = 1, 
                              n_cv_folds = 5)

  temp_fun <- function(x, fit){
    fit_temp <- fit[[x]]
    draws <- fit_temp$draws(c("beta", "Omega",
                              "logprob", "cv_fold",
                              "log_loss",
                              "y_true"))
    draws_df <- as_draws_df(draws)
    draws_df
  }
  
  output_list <- map(1:5, ~temp_fun(.x, fit = fit))

  output_df <- tidy_output_MV(output_list = output_list,
                              swab_type = swab_type,
                              fit_type = model_type,
                              pred_names = covariate_names,
                              response_names = c("result", symptom_nam),
                              n_cv_folds = 5)
  output_df
}

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
# Need to convert these values to false & true positives and false & true negatives
ROC_diagnose <- function(validation_df, prob_range, increment_size){
  classification_df <- threshold(validation_df = validation_df, 
                                 prob_range = c(prob_range[1], prob_range[2]),
                                 increment_size = increment_size)
  # Subtract true value from each classified value
  classification_df <- classification_df %>% 
    mutate(across(starts_with("threshold"), ~ .x - Truth))
  # Extract classified values and convert to false & true positives and false & true negatives
  classifications <- classification_df %>% select(starts_with("threshold"))
  # classifications[classifications == 1] <- "FalsePos"
  # classifications[classifications == -1] <- "FalseNeg"
  classifications[classifications == 0 & classification_df$Truth == 1] <- 10 #"TruePos" 
  classifications[classifications == 0 & classification_df$Truth == 0] <- -10 #"TrueNeg" 
  
  diagnosed_df <- cbind(classification_df %>% 
                          select(SwabType, FitType, CV, Iter, Chain), 
                        classifications)
  
  
  diagnosed_grouped_df <- diagnosed_df %>% 
    pivot_longer(cols = starts_with("threshold"), 
                 names_to = "threshold", 
                 values_to = "classification") 

  diag_df <- lazy_dt(diagnosed_grouped_df) %>%
    group_by(SwabType, FitType, Iter, Chain, threshold) %>%
    summarise(FalsePos = sum(classification == 1), #"FalsePos"),
              TruePos = sum(classification == 10), #"TruePos"),
              FalseNeg = sum(classification == -1), #"FalseNeg"),
              TrueNeg = sum(classification == -10)) %>% #"TrueNeg")) %>%
    ungroup() %>%
    as.data.frame()
 
  diag_df$FalseNegRate = diag_df$FalseNeg/(diag_df$FalseNeg + diag_df$TruePos)
  diag_df$FalsePosRate = diag_df$FalsePos/(diag_df$FalsePos + diag_df$TrueNeg)
  diag_df$TrueNegRate = diag_df$TrueNeg/(diag_df$FalsePos + diag_df$TrueNeg)
  diag_df$TruePosRate = diag_df$TruePos/(diag_df$FalseNeg + diag_df$TruePos)

  # diag_df <- raw_diag_df %>%
  #   group_by(SwabType, FitType, CV, Iter, Chain, threshold) %>%
  #   mutate(FalseNegRate = FalseNeg/(FalseNeg + TruePos),
  #          FalsePosRate = FalsePos/(FalsePos + TrueNeg),
  #          TrueNegRate = TrueNeg/(FalsePos + TrueNeg),
  #          TruePosRate = TruePos/(FalseNeg + TruePos)) %>%
  #   as.data.frame()
  diag_df
}

auc_calc <- function(validation_df, prob_range, increment_size){
  
  classification_df <- threshold(validation_df = validation_df, 
                                             prob_range = c(prob_range[1], prob_range[2]),
                                             increment_size = increment_size)
  auc_tmp <- classification_df %>% 
    pivot_longer(cols = starts_with("threshold"), 
                 names_to = "threshold", 
                 values_to = "classification") %>%
    na.omit() %>%
    nest(data = -c(SwabType, FitType, Chain, Iter))  %>% 
    mutate(auc = map_dbl(data, ~ performance(
      prediction(
        predictions = .x$classification, labels = .x$Truth), "auc")@y.values %>% 
        unlist(use.names = FALSE)))

auc_tmp
}

