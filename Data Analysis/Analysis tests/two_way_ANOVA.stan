data {
  int<lower=0> N;
  vector[N] y;
  vector[N] treatment;
  vector[N] environment;
}

parameters {
  real mu;
  real beta_int;
  real beta_env;
  real beta_t;
  real sigma;
}


model {
  mu ~ normal(40, 6); // adjust based on outcome variable 
  sigma ~ exponential(1);
  
  beta_int ~ normal(0,7);
  beta_env ~ normal(0,7);
  beta_t ~ normal(0,7);
  
  for (i in 1:N){
    real out_mu = mu + beta_t*treatment[i] + beta_env*environment[i] + beta_int*(environment[i]*treatment[i]); 
    y[i] ~ normal(out_mu, sigma);
  }

}

