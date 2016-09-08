library(ggplot2)
library(MASS)
library(RColorBrewer)
library(car)
library(reshape2)

library(foreach)
library(doParallel)

source('sussex/ICICworkshop/SN_functions.R')
source('sussex/ICICworkshop/MCMC_functions.R')


# read SNIa data
dat <- read.table(file = 'sussex/ICICworkshop/Pre-workshop Exercise/jla_mub_0.txt')
names(dat) <- c("z","mu")

# read covariance table and assign to matrix
cov <- read.table('sussex/ICICworkshop/Pre-workshop Exercise/jla_mub_covmatrix.txt')
cov <- matrix(as.matrix(cov),nrow=31,ncol=31)


cl<-makeCluster(2)
registerDoParallel(cl)


# number of iterations
N <- 6000

# number of chains
n = 4

chains <- foreach(j = seq(n), .combine = rbind)  %dopar% {
  # initial values
  h = runif(1)
  omegaM = runif(1)
  omegaL = runif(1)
  
  # initialise chain vectors
  h_vector = rep(NA,N)
  omegaM_vector = rep(NA,N)
  omegaL_vector = rep(NA,N)
  
  for(i in seq(N)){
    # sample new values
    params = sampler(h, omegaM, omegaL, z=dat[,1], mu=dat[,2], cov=cov, width=0.003)
    
    # save to chain vector, update current values
    h = h_vector[i] = params[[1]]
    omegaM = omegaM_vector[i] = params[[2]]
    omegaL = omegaL_vector[i] = params[[3]]
  }
  
  list(h_vector, omegaM_vector, omegaL_vector)
}


stopImplicitCluster()


## Gelman-Rubin diagnostic

source('sussex/ICICworkshop/GelmanRubin_functions.R')

Vhat_vector <- apply(chains, 2, function(x) Vhat(x))
W_vector <- apply(chains, 2, function(x) W(x))

apply(chains, 2, function(x) sqrt(Vhat(x) / W(x)) )



## Dataframes

melt_dataframes <- function(chains, burnin=F){
  
  h_df = data.frame(chains[,1])
  omegaM_df = data.frame(chains[,2])
  omegaL_df = data.frame(chains[,3])
  
  if(burnin){
    h_df = h_df[burnin:nrow(h_df),]
    omegaM_df = h_df[burnin:nrow(omegaM_df),]
    omegaL_df = h_df[burnin:nrow(omegaL_df),]
  }
  
  # add id column
  ids <- seq.int(nrow(h_df))
  h_df$id = omegaL_df$id = omegaM_df$id = ids
  
  # add name column
  h_df$param = 'h'
  omegaM_df$param = 'omegaM'
  omegaL_df$param = 'omegaL'
  
  h_df <- melt(h_df, id.vars=c('id','param'))
  omegaM_df <- melt(omegaM_df, id.vars=c('id','param'))
  omegaL_df <- melt(omegaL_df, id.vars=c('id','param'))
  
  rbind(h_df, omegaM_df, omegaL_df)
}

df <- melt_dataframes(chains)

## Trace plots
m <- ggplot(df) + geom_line(aes(x=id, y=value, colour=variable))
m + facet_wrap(~ param)

## BURN-IN
burnin <- 1000
df <- melt_dataframes(chains, burnin=burnin)

scatterplotMatrix(df, diagonal="density", smooth=F, reg.line = F, ellipse=T)





