
Dl <- function(z, omegaM){
  2.9979e8/100000 * (1+z) * (nu1(1,omegaM) - nu1(1/(1+z),omegaM))
}


nu2 <- function(a, omegaM){
  O3 <- (1 - omegaM)/omegaM
  
  return(2 * sqrt(O3^3 + 1) * ( (1/a^4) - (0.154)*(O3/a^3) + 0.4304*((03^2)/(a^2)) 
                        + 0.19097*(O3^3)/a + 0.066941*(O3^4) )^(-1/8))
}

nu1 <- function(a, omegaM){
  O3 <- (1 - omegaM)/omegaM
  
  return(2 * sqrt(O3 + 1) * ( (1/(a^4)) - (0.154)*(O3^(1/3)/a^3) + 0.4304*((03^(2/3))/(a^2)) 
                        + 0.19097*(O3)/a + 0.066941*(O3^(4/3)) )^(-1/8))
}

mu <- function(z, omegaM, h=0.7){
  25 - 5*log10(h) + 5*log10(Dl(z, omegaM))
}

# 1) 
# plot analytical prediction
z <- seq(0,2,0.01)

plot(z, mu(z, omegaM=0.1), type='l')
lines(z, mu(z, omegaM=0.2), col=2)
lines(z, mu(z, omegaM=0.3), col=3)
lines(z, mu(z, omegaM=0.4), col=4)


# 2) 
# read SNIa data
dat <- read.table(file = 'sussex/ICICworkshop/Pre-workshop Exercise/jla_mub_0.txt')
names(dat) <- c("z","mu")

points(dat[['z']],dat[['mu']], pch=5)

# 3)
# generate random SN

z <- runif(20,min=0,max=2)
appm <- mu(z, h=0.7, omegaM=0.3)

appm = appm + rnorm(20,mean=0,sd=0.1)

# 4)
points(z, appm, pch = 4)



# 5)
# integral form for flat universes

Dl_integ <- function(z, omegaM, h=0.7){
  return( 3000 / h * (1+z) * integrate(integrand, omegaM = omegaM, lower = 0, upper = z)$value )
}

integrand <- function(z, omegaM){
  1 / sqrt(omegaM * ((1+z)^3) + 1 - omegaM)
}

mu_integ <- function(z, omegaM, h=0.7){
  25 - 5*log10(h) + 5*log10(Dl_integ(z, omegaM))
}

z <- seq(0,2,0.01)
lapply(z, function(x) mu_integ(x, 0.3))


# integral form for general (non-flat) universes

Dl_general <- function(z, omegaM, omegaL, h=0.7){
  
  omega = omegaM + omegaL

  rint <- integrate(r_integrand, omegaM = omegaM, omegaL = omegaL, lower=0, upper = z)$value
  
  if(omega > 1){
    r <- sqrt(abs(1-(omegaM+omegaL))) * rint
    return( (1+z) * 2.9979e8 / h*100000 / sqrt(abs(1-omega)) * sin(r) )
    
  }else if(omega < 1){
    r <- sqrt(abs(1-(omegaM+omegaL))) * rint
    return( (1+z) * 2.9979e8 / h*100000 / sqrt(abs(1-omega)) * sinh(r) )
    
  }else{
    return( (1+z) * 2.9979e8 * rint / h*100000 )
    
  }
  
  #return( (1+z) * 2.9979e8 / h*100000 / sqrt(abs(1-omega)) * Sk )
}

rz <- function(z, omegaM, omegaL){
  integrate(integrand, omegaM=omegaM, omegaL=omegaL, lower=0, upper=z)$value
}

r_integrand <- function(z, omegaM, omegaL){
  1 / sqrt(omegaM * (1+z)^3 + omegaL + (1-(omegaM+omegaL))*(1+z)^2)
}

Dl_general(0.5, 0.3, 0.7)

