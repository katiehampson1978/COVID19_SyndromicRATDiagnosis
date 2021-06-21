// This code is adapted from Ben Goodrich's https://github.com/stan-dev/example-models/blob/master/misc/multivariate-probit/probit-multi-good.stan
// Any mistakes are my own.

// Create multivariate probit function - this improves code readability and allows
// for easy conversion to reduce_sum within chain parallelisation if needed (not
// implemented here)
functions {
   real multiprobit_lpmf(
   int[] y, // Integer array of binary response dimensions
   vector x, // Covariate vector
   matrix beta, // Matrix of coefficients
   real[] u, // Nuiance parameter to help with numerical instability
   int D, // Number of response dimensions
   matrix L_Omega // Cholesky factor for correlation matrix
   ){
      vector[D] mu; // Sum of linear predictor
      vector[D] z; // Continuous latent state
      vector[D] logLik; // Log likelihood
      real prev; // Boundary of previous dimension
      mu = beta * x; // Sum linear predictor
      prev = 0; // Starting boundary condition
     for (d in 1:D) { // Phi and inv_Phi may overflow and / or be numerically inaccurate
        real bound; // threshold at which utility = 0
        bound = Phi( -(mu[d] + prev) / L_Omega[d,d] );
        if (y[d] == 1) {
          real t;
          t = bound + (1 - bound) * u[d];
          z[d] = inv_Phi(t);       // implies utility is positive
          logLik[d] = log1m(bound);  // Jacobian adjustment
        }
        else {
          real t;
          t = bound * u[d];
          z[d] = inv_Phi(t);     // implies utility is negative
          logLik[d]= log(bound);  // Jacobian adjustment
        }
        if (d < D) prev = L_Omega[d+1,1:d] * head(z, d);
        // Jacobian adjustments imply z is truncated standard normal
        // thus utility --- mu + L_Omega * z --- is truncated multivariate normal
      }
      return sum(logLik);
   }
}
data {
  int<lower=1> K; // Number of covariates
  int<lower=1> D; // Number of response dimensions
  int<lower=0> N_train; // Size of training set
  int<lower=0> N_test; // Size of testing set
  int<lower=1> fold; // Fold ID - helps ease code admin
  int<lower=0,upper=1> y_train[N_train,D]; // Response data for training set
  int<lower=0,upper=1> y_test[N_test,D]; // Response data for testing set
  int<lower=0,upper=1> y_test_ll[N_test]; // Duplicate response data for testing set
  vector[K] x_train[N_train]; // Covariate data for training set
  vector[K] x_test[N_test]; // Covariate data for testing set 
  int<lower=0,upper=1> RAT_train[N_train]; // RAT result for training set
  int<lower=0,upper=1> RAT_test[N_test]; // RAT result for testing set
}
transformed data{
    // Create test set where COVID status is zero and one for comparison later
    // Move out of Stan
    int<lower=0,upper=1>  y_zero[N_test,D]; // Response data if COVID status is 0
    int<lower=0,upper=1>  y_one[N_test,D]; // Response data if COVID status is 1
    y_zero = y_test;
    y_one = y_test;
    for(i in 1:N_test){
      y_zero[i,1] = 0;
      y_one[i,1] = 1;
    }
}
parameters {
  matrix[D,K] beta; // Covariate coefficients
  cholesky_factor_corr[D] L_Omega; // Cholesky factor for correlation matrix
  real<lower=0,upper=1> u_train[N_train,D]; // nuisance that absorbs inequality constraints
  real<lower=0,upper=1> u_test[N_test,D]; // nuisance that absorbs inequality constraints 
  real alpha; // Error rate for RAT positive patients 
}
model {
  L_Omega ~ lkj_corr_cholesky(1); // Prior for correlation matrix prior
  to_vector(beta) ~ normal(0, 1); // Prior for covariate coefficients
    alpha ~ normal(0,1); // Prior for RAT positive error rates 
  // Fit model to training data
  for(m in 1:N_train){
  if(RAT_train[m]){ // If patient is RAT positive, fit intercept only logistic regression
    y_train[m,1] ~ bernoulli_logit(alpha);
  }else{ // If patient is RAT positive, fit multivariate probit 
        target += multiprobit_lpmf(y_train[m] | x_train[m], beta, u_train[m], D, L_Omega);
}
}
}
generated quantities {
  real logprob[N_test]; // Log probability for each test profile
  corr_matrix[D] Omega; // Correlation matrix
  real log_loss_tmp[N_test]; // Per individual log loss (cross entropy) 
  real log_loss; //Model log loss (cross entropy) 
  real<lower=0,upper=1> y_true[N_test]; // True result for testing set
  int<lower=1> cv_fold; // Cross validation fold ID
  Omega = multiply_lower_tri_self_transpose(L_Omega); // Go from cholesky factor to full matrix
    // Generate predictions for test data
  for(j in 1:N_test){
  if(RAT_test[j]){  // If patient is RAT positive, predict from logistic regression
      logprob[j] = log(inv_logit(alpha));
    }else{ // If patient is RAT positive, predict from multivariate probit
      vector[2] lp; // likelihood storage
      // likelihood of COVID = 0
      lp[1] = multiprobit_lpmf(y_zero[j]  | x_test[j], beta, u_test[j], D, L_Omega);
      // likelihood of COVID = 1
      lp[2] = multiprobit_lpmf(y_one[j]  | x_test[j], beta, u_test[j], D, L_Omega);
      // Log probability that COVID = 1
      logprob[j] = (lp[2] - log_sum_exp(lp));
      // Log loss for point j
      log_loss_tmp[j] = binary_log_loss(y_test_ll[j], exp(logprob[j]));
    }
    log_loss_tmp[j] = binary_log_loss(y_test_ll[j], exp(logprob[j]));
  }
    // Model log loss
    log_loss = sum(log_loss_tmp)/N_test;
    // Some inefficient ways of keeping track of the true values and CV ID 
    y_true = y_test_ll;
    cv_fold = fold;
}
