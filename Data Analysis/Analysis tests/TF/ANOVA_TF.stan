data {
  int<lower=0> N;
  int<lower =0> J; 
  
  array[N] int animal;
  vector[N] y;
  vector[N] treatment;
  vector[N] environment;
  
}

parameters {
  real mu_alpha;
  real<lower=0> sigma_alpha;
  real beta_int;
  real beta_env;
  real beta_t;
  real<lower=0> sigma;
  vector[J] mu_raw;
}

transformed parameters{ 
  vector[J] mu;
  mu = mu_alpha + sigma_alpha*mu_raw;
  }

model {
  sigma ~ exponential(1);
  
  mu_alpha ~ normal(0,1);
  sigma_alpha ~ exponential(1);
  
  beta_int ~ normal(0,1);
  beta_env ~ normal(0,1);
  beta_t ~ normal(0,1);
  
  mu_raw ~ normal(0, 1);
  
  for (i in 1:N){
    
    real out_mu = mu[animal[i]] + beta_t*treatment[i] + beta_env*environment[i] + beta_int*(environment[i]*treatment[i]); 
    y[i] ~ normal(out_mu, sigma);
  }

}
generated quantities{
  real ME_t = beta_t + 0.5*beta_int;
  real ME_env = beta_env + 0.5*beta_int;
  
  real prior_ME = normal_rng(0,1.25);
  real prior_int = normal_rng(0,1);
  
  real PE_diff = beta_t; //use prior N(0,1)
  real EE_diff = beta_t + beta_int;
  real SAL_diff = beta_env; //use prior N(0,1)
  real CVAD_diff = beta_env + beta_int;
  real int_diff = beta_env+beta_t+beta_int;
  
  real prior_EE_CVAD_diff = normal_rng(0,2);
  real prior_int_diff= normal_rng(0,3);
}

