library(ggplot2)
library(MASS)
library(RColorBrewer)
library(car)
library(reshape2)

source('sussex/ICICworkshop/SN_functions.R')
source('sussex/ICICworkshop/MCMC_functions.R')

z <- seq(0,2,0.01)
plot(z,lapply(z, function(x) mu_integ(x, 0.3, 0.7)), ylab = 'mu')


# read SNIa data
dat <- read.table(file = 'sussex/ICICworkshop/Pre-workshop Exercise/jla_mub_0.txt')
names(dat) <- c("z","mu")

# read covariance table and assign to matrix
cov <- read.table('sussex/ICICworkshop/Pre-workshop Exercise/jla_mub_covmatrix.txt')
cov <- matrix(as.matrix(cov),nrow=31,ncol=31)


# initial values
h = runif(1)
omegaM = runif(1)
omegaL = runif(1)

# number of iterations
N <- 10000

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

## save as data frame
df <- data.frame(h = h_vector, omegaM = omegaM_vector, omegaL = omegaL_vector)

## Trace plots
df$id <- seq(nrow(df))

m <- ggplot(melt(df, id.vars = 'id')) + geom_line(aes(x=id, y=value, colour=variable))
m + facet_wrap(~ variable)

## acceptance rate
print(length(unique(h_vector))/length(h_vector))

## BURN-IN
burnin <- 5000
df = df[burnin:nrow(df),]

## scatter plots
scatterplotMatrix(df[,1:3], diagonal="density", smooth=F, reg.line = F, ellipse=T)

## base plotting
plot(h_vector, omegaM_vector,type = 'b')

# pretty colous for 2d density
k <- 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))

z <- kde2d(h_vector, omegaM_vector, n = 50)
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)

## ggplot2 plotting
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


## covariances
var(omegaM_vector, h_vector)
var(omegaL_vector, h_vector)
var(omegaM_vector, omegaL_vector)

