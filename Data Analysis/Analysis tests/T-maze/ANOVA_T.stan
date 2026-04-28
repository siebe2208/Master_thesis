data {
  int<lower=0> N;
  int<lower =0> J; 
  
  array[N] int animal;
  vector[N] y;
  vector[N] treatment;
  vector[N] environment;
  vector[N] Stage;
  
}

parameters {
  real mu_alpha;
  real<lower=0> sigma_alpha;
  real beta_env;
  real beta_day;
  real beta_e_t;
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

  beta_env ~ normal(0,1);
  beta_t ~ normal(0,1);
  beta_day ~ normal(0,1);
  
  beta_e_t ~ normal(0,1);
  
  mu_raw ~ normal(0, 1);
  
  for (i in 1:N){
    
    real out_mu = mu[animal[i]] + beta_t*treatment[i] + beta_env*environment[i] + beta_day * Stage[i] + beta_e_t*(environment[i]*treatment[i]); 
    y[i] ~ normal(out_mu, sigma);
  }

}
generated quantities{
  real prior = normal_rng(0,1);
}

