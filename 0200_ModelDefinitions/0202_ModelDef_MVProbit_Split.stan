// This code is adapted from Ben Goodrich's https://github.com/stan-dev/example-models/blob/master/misc/multivariate-probit/probit-multi-good.stan
// Any mistakes are my own.
functions {
   real multiprobit_lpmf(
   int[] y,
   vector x,
   matrix beta,
   real[] u,
   int D,
   matrix L_Omega
   ){
      vector[D] mu;
      vector[D] z;
      vector[D] logLik;
      real prev;
      mu = beta * x;
      prev = 0;
     for (d in 1:D) { // Phi and inv_Phi may overflow and / or be numerically inaccurate
        real bound; // threshold at which utility = 0
        bound = Phi( -(mu[d] + prev) / L_Omega[d,d]  );
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
  int<lower=1> K;
  int<lower=1> D;
  int<lower=0> N_train;
  int<lower=0> N_test;
  int<lower=1> fold;
  int<lower=0,upper=1> y_train[N_train,D];
  int<lower=0,upper=1> y_test[N_test,D];
  int<lower=0,upper=1> y_test_ll[N_test];
  vector[K] x_train[N_train];
  vector[K] x_test[N_test];
  int<lower=0,upper=1> RAT_train[N_train];
  int<lower=0,upper=1> RAT_test[N_test];
}
transformed data{
    // Create test set where COVID status is zero and one for comparison later
    // Move out of Stan
    int<lower=0,upper=1>  y_zero[N_test,D];
    int<lower=0,upper=1>  y_one[N_test,D];
    y_zero = y_test;
    y_one = y_test;
    for(i in 1:N_test){
      y_zero[i,1] = 0;
      y_one[i,1] = 1;
    }
}
parameters {
  matrix[D,K] beta;
  cholesky_factor_corr[D] L_Omega;
  real alpha;
  real<lower=0,upper=1> u_train[N_train,D]; // nuisance that absorbs inequality constraints
  real<lower=0,upper=1> u_test[N_test,D]; 
}
model {
  L_Omega ~ lkj_corr_cholesky(1);//16-K);
  to_vector(beta) ~ normal(0, 1);
  alpha ~ normal(0,1);
  for(m in 1:N_train){
  if(RAT_train[m]){
    y_train[m,1] ~ bernoulli_logit(alpha);
  }else{
        target += multiprobit_lpmf(y_train[m] | x_train[m], beta, u_train[m], D, L_Omega);
}
}
}
generated quantities {
  real logprob[N_test];
  corr_matrix[D] Omega;
  real log_loss_tmp[N_test];
  real log_loss;
  real beta_bar;
  real<lower=0,upper=1> y_true[N_test];
  int<lower=1> cv_fold;
  Omega = multiply_lower_tri_self_transpose(L_Omega);
  for(j in 1:N_test){
if(RAT_test[j]){
      logprob[j] = log(inv_logit(alpha));
    }else{
      vector[2] lp;
      lp[1] = multiprobit_lpmf(y_zero[j]  | x_test[j], beta, u_test[j], D, L_Omega);
      lp[2] = multiprobit_lpmf(y_one[j]  | x_test[j], beta, u_test[j], D, L_Omega);
      logprob[j] = (lp[2] - log_sum_exp(lp));
    }
    log_loss_tmp[j] = binary_log_loss(y_test_ll[j], exp(logprob[j]));
  }
    log_loss = sum(log_loss_tmp)/N_test;
    y_true = y_test_ll;
    cv_fold = fold;
}
