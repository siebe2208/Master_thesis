## Frequentist analysis

two_way_aov = function(data, fac1, fac2, out){
  
  f = as.formula(paste(out, "~", fac1, "*", fac2))
  ANOVA = aov(f, data = data)
  coef = coef(ANOVA)
  
  sum = summary(ANOVA)[[1]]
  
  
  res = sum["Residuals", ] %>% select(`Sum Sq`) %>% as.numeric()
  sum = sum %>% mutate(ES = `Sum Sq`/(`Sum Sq` + res)) # Add effect size 
  
  result = list(summary = sum, coef = coef, fit = ANOVA)
  return(result)
}


## bayesian two_way

two_way_bayes = function(model, data, out, path, SD){
  
  normalized = as.numeric(scale(data[[out]]))
  
  data_stan = list(N = nrow(data), y = normalized, treatment = as.numeric(data$treatment)-1, environment = as.numeric(data$env)-1)
  
  if (file.exists(path)) {fit = readRDS(path)} else {
    fit = mod$sample(data = data_stan,iter_warmup = 1000,iter_sampling = 2000,chains = 4,parallel_chains = 4,
                     output_dir = dirname(path))
    saveRDS(fit, path)}
  
  bayes = fit$summary()
  
  post_int = as_draws_df(fit$draws("beta_int"))$beta_int
  post_t = as_draws_df(fit$draws("ME_t"))$ME_t
  post_env = as_draws_df(fit$draws("ME_env"))$ME_env
  post_SAL_diff = as_draws_df(fit$draws("SAL_diff"))$SAL_diff
  post_CVAD_diff = as_draws_df(fit$draws("CVAD_diff"))$CVAD_diff
  post_PE_diff = as_draws_df(fit$draws("PE_diff"))$PE_diff
  post_EE_diff = as_draws_df(fit$draws("EE_diff"))$EE_diff
  post_int_diff = as_draws_df(fit$draws("int_diff"))$int_diff
  
  prior_ME = as_draws_df(fit$draws("prior_ME"))$prior_ME
  prior_int = as_draws_df(fit$draws("prior_int"))$prior_int
  prior_CVAD_EE = as_draws_df(fit$draws("prior_EE_CVAD_diff"))$prior_EE_CVAD_diff
  prior_int_diff = as_draws_df(fit$draws("prior_int_diff"))$prior_int_diff
  
  
  BF_int = get_bf(prior_int, post_int)
  BF_t = get_bf(prior_ME,post_t)
  BF_env = get_bf(prior_ME,post_env)
  BF_SAL_diff = get_bf(prior_int, post_SAL_diff)
  BF_CVAD_diff = get_bf(prior_CVAD_EE, post_CVAD_diff)
  BF_PE_diff = get_bf(prior_int, post_PE_diff)
  BF_EE_diff = get_bf(prior_CVAD_EE, post_EE_diff)
  BF_int_diff = get_bf(prior_int_diff, post_int_diff)
  
  
  BF_tbl = tibble(
    effect = c("interaction", "treatment", "environment", "SAL diff", "CVAD diff", "PE diff", "EE diff", "int diff"),
    BF     = c(BF_int, BF_t, BF_env, BF_SAL_diff, BF_CVAD_diff, BF_PE_diff, BF_EE_diff, BF_int_diff))
  
  result = list(summary = bayes, BF=BF_tbl)
  
  return(result)
}

##Bayes factor fun
get_bf = function(prior,post){
  post_fit  <- logspline(post)
  prior_fit <- logspline(prior)
  
  post_at_0  <- dlogspline(0, post_fit)
  prior_at_0 <- dlogspline(0, prior_fit)
  
  BF_01 <- prior_at_0 / post_at_0
  
  return(BF_01)
}

# Multilevel ANOVA
ANOVA_multi = function(data, fac1, fac2, out, level){
  f = as.formula(paste(out, "~", fac1, "*", fac2, "+(1|", level, ")"))
  
  fit = lmer(f, data)
  
  sum = summary(fit)
  
  ES = r2beta(fit, method = "nsj")
  
  return(list(summary = sum, ES =ES, fit = fit))
}
