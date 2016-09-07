
# likelihood
logLikelihood <- function(h, omegaM, omegaL, z, mu, cov){
  
  -0.5 * ((mu - unlist(lapply(z, function(x) mu_integ(x, omegaM = omegaM, omegaL = omegaL, h=h))))
          %*% solve(cov)
          %*% (mu - unlist(lapply(z, function(x) mu_integ(x, omegaM = omegaM, omegaL = omegaL, h=h)))))
  
}

# sample new values and assess whether they should be accepted
sampler <- function(h, omegaM, omegaL, z, mu, cov, width=0.01){
  
  repeat{
    h_trial = rnorm(1, mean = h, sd = width)
    omegaM_trial = rnorm(1, mean = omegaM, sd = width)
    omegaL_trial = rnorm(1, mean = omegaL, sd = width)
    if(h_trial > 0 & omegaM_trial > 0){
      break
    }
  }
  
  accept = min(1,exp(logLikelihood(h_trial, omegaM_trial, omegaL_trial, z=z, mu=mu, cov=cov) 
                     - logLikelihood(h, omegaM, omegaL, z=z, mu=mu, cov=cov)))
  
  if(runif(1)<accept){
    return(list(h_trial, omegaM_trial, omegaL_trial))
  }else{
    return(list(h, omegaM, omegaL))
  }
}