

Dl_general <- function(z, omegaM, omegaL, h=1){
  
  omega = omegaM + omegaL
  
  rint <- integrate(r_integrand, omegaM = omegaM, omegaL = omegaL, lower=0, upper = z)$value
  
  if(omega > 1){
    r <- sqrt(abs(1-(omegaM+omegaL))) * rint
    return( (1+z) * 2.9979e8 / (h*100000) / sqrt(abs(1-omega)) * sin(r) )
  }else if(omega < 1){
    r <- sqrt(abs(1-(omegaM+omegaL))) * rint
    return( (1+z) * 2.9979e8 / (h*100000) / sqrt(abs(1-omega)) * sinh(r) )
  }else{
    return( ((1+z) * 2.9979e8 * rint) / (h*100000) )
  }
}

r_integrand <- function(z, omegaM, omegaL){
  1 / sqrt(omegaM * (1+z)^3 + omegaL + ((1-(omegaM+omegaL)) * (1+z)^2) )
}

mu_integ <- function(z, omegaM, omegaL, h=0.7){
  25 - 5*log10(h) + 5*log10(Dl_general(z, omegaM, omegaL, h=1))
}

