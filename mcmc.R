library(ggplot2)
library(MASS)
library(RColorBrewer)


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

z <- seq(0,2,0.01)
plot(z,lapply(z, function(x) mu_integ(x, 0.3, 0.7)), ylab = 'mu')


# read SNIa data
dat <- read.table(file = 'sussex/ICICworkshop/Pre-workshop Exercise/jla_mub_0.txt')
names(dat) <- c("z","mu")

# read covariance table and assign to matrix
cov <- read.table('sussex/ICICworkshop/Pre-workshop Exercise/jla_mub_covmatrix.txt')
cov <- matrix(as.matrix(cov),nrow=31,ncol=31)


# likelihood
logLikelihood <- function(h, omegaM, omegaL, z, mu, cov){
  
  -0.5 * ((mu - unlist(lapply(z, function(x) mu_integ(x, omegaM = omegaM, omegaL = omegaL, h=h))))
              %*% solve(cov)
              %*% (mu - unlist(lapply(z, function(x) mu_integ(x, omegaM = omegaM, omegaL = omegaL, h=h)))))
  
}


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


# initial values
h = runif(1)
omegaM = runif(1)
omegaL = runif(1)

# number of iterations
N <- 1000

# initialise chain vectors
h_vector = rep(NA,N)
omegaM_vector = rep(NA,N)
omegaL_vector = rep(NA,N)

for(i in seq(N)){
  # sample new values
  params = sampler(h, omegaM, omegaL, z=dat[,1], mu=dat[,2], cov=cov, width=0.008)
  
  # save to chain vector, update current values
  h = h_vector[i] = params[[1]]
  omegaM = omegaM_vector[i] = params[[2]]
  omegaL = omegaL_vector[i] = params[[3]]
}

## acceptance rate
print(length(unique(h_vector))/length(h_vector))

## base plotting
plot(h_vector, omegaM_vector,type = 'b')

# pretty colous for 2d density
k <- 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))

z <- kde2d(h_vector, omegaM_vector, n = 50)
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)

## ggplot2 plotting
df <- data.frame(h = h_vector, omegaM = omegaM_vector, omegaL = omegaL_vector)

m <- ggplot(df, aes(x = omegaL, y = omegaM)) + geom_point() #+ xlim(0.2,0.5) + ylim(0.6,0.8)
m + geom_density2d(aes(colour = ..level..))
#m + stat_density2d(aes(fill = ..level..), geom="polygon")


## Parameter values
mean(h_vector)
var(h_vector)

mean(omegaM_vector)
var(omegaM_vector)

mean(omegaL_vector)
var(omegaL_vector)


# covariances
var(omegaM_vector, h_vector)
var(omegaL_vector, h_vector)
var(omegaM_vector, omegaL_vector)








