data {
  int<lower=0> Nobs; 
  real TempStd[Nobs]; 
  int<lower=0, upper = 1> Y[Nobs]; 
} 
parameters {
  real alpha; 
  real beta;
  real delta;
} 

transformed parameters {
  real p[Nobs];
  
  for (i in 1:Nobs) {
    
    p[i] = inv_logit(
      
      alpha*(TempStd[i])^2 + beta*TempStd[i] + delta
      
      );
      //beta - alpha*(TempStd[i]^2)
      //exp(-((TempStd[i] - alpha)^2) / (2 * beta^2)));
    
  }
}

model {
  // priors
  alpha ~ normal(0, 1); 
  beta ~ normal(0, 1);
  delta ~ normal(0, 1);
  
  // likelihood
  
  
  Y ~ bernoulli( p );
  
}

