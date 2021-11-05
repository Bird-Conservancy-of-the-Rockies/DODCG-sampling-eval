##############################################
### Simulate sample sizes needed for IMBCR
### Account for density, heterogeneity
### Monitor precision, bias
##############################################

rm(list = ls())
library(jagsUI)

setwd('C:/Users/Quresh.Latif/files/projects/power_analyses/DOD')
mod.file <- "DODCG-sampling-eval/dens_mod_qzhao.jags"

nsim <- 100
#==============
# Basic values
#==============
D.true <- 1
sample.true <- 100 #Number of grids to sample
pts.per.grid <- 15 #Number of points in a grid
cutoff <- 100
sigma <- 55

#____ Detection parameters _____#
area.circle <- pi * (cutoff / 1000) ^ 2 #area of point count circle in km2
nG <- 10 #Number of distance classes
breaks <- seq(0, cutoff, length.out=nG + 1) #Break points for distance classes
area.band <- pi*breaks[-1]^2 - pi*breaks[-(nG+1)]^2
area.prop <- area.band / sum(area.band)
   
int.true <- (sigma^2*(1 - exp(-breaks[-1]^2/(2*sigma^2))) - sigma^2*(1 - exp(-breaks[-(nG+1)]^2/(2*sigma^2)))) #integral of half-normal function
p.true <- 2*pi*int.true/area.band #Complete half-normal calculations
pi.true <- p.true*area.prop #Correct for availability
overall.pi.true <- sum(pi.true) #Overall detection probability

#======================
# Define model in Jags
#======================
sink(file=mod.file)
cat("
  model{
      
  # Priors
  #lambda ~ dgamma(.01, .01)
  for(k in 1:sample.true) {
    lambda[k] ~ dgamma(.01, .01)
  }
  sigma ~ dgamma(.01, .01)

  for(k in 1:nG) {
    int[k] <- sigma^2 * (1-exp(-pow(breaks[k+1],2)/(2*sigma^2))) - sigma^2 *
      (1-exp(-pow(breaks[k],2)/(2*sigma^2)))
    pi[k] <- 2*3.141593*int[k]/area.band[k]*area.prop[k]
  } # k
  pi.sum <- sum(pi[1:nG])
  for(k in 1:nG) {
    pic[k] <- pi[k] / pi.sum
  } # k

  # Process model
  for(i in 1:n.pt) {
    N.pt[i] ~ dpois(lambda[grid.pt.ind[i]])
  } # i

  # Observation model
  for(i in 1:n.pt) {
    y.sum[i] ~ dbinom(pi.sum, N.pt[i])
    y.bin[i,1:nG] ~ dmultinom(pic[1:nG], y.sum[i])
  } # i

  } # model
", fill=TRUE)
sink()

#===============
# Simulate data
#===============
for (s in 1:nsim) {
  n.pt <- sample.true * pts.per.grid #Number of points
  grid.pt.ind <- rep(1:sample.true, each = pts.per.grid) #Index mapping point to grid
  lambda.true <- exp(log(D.true * area.circle))
  N.true <- rpois(n.pt, lambda.true) #Sample point-level abundance
  y.sum <- rbinom(n.pt, N.true, overall.pi.true) #Sample number of observations at each point
  y.bin <- matrix(NA, n.pt, nG) #Create matrix to hold distance class data
  for (j in 1:n.pt) {
    y.bin[j,] <- rmultinom(1, y.sum[j], pi.true)
  } # close j loop


  #==========
  # run Jags
  #==========
  data <- list(nG=nG, n.pt=n.pt, breaks=breaks, area.band=area.band,
               area.prop=area.prop, y.bin=y.bin, grid.pt.ind=grid.pt.ind,
               sample.true=sample.true)

  inits <- function() list(lambda=rep(lambda.true,sample.true),
                           sigma=sigma, y.sum=apply(y.bin, 1, sum),
                           N.pt=apply(y.bin, 1, sum)*2+2)

  parms <- c('lambda', 'sigma')

  fit <- jags(data, inits, parms, mod.file, 
  #            n.chains=1, n.adapt=100, n.burnin=100, n.iter=200, n.thin=1, 
            n.chains=3, n.burnin=1000, n.iter=2000, n.thin=1, 
            parallel=TRUE)

  save(fit, file=paste(c('fit_', s, '.RData'), collapse=''))

} # s

print(fit)
par(mfrow=c(1,3))
traceplot(fit)

#===============
# Check results
#===============
library(vioplot)

lambda.post <- sigma.post <- matrix(, 3000, nsim)
for (s in 1:nsim) {
  load(file=paste(c('fit_', s, '.RData'), collapse=''))
  lambda.post[,s] <- apply(fit$sims.list$lambda / area.circle, 1, mean)
  sigma.post[,s] <- fit$sims.list$sigma
} # s

par(mfrow=c(1,2))
plot(1, ylim=c(0, 2), type='n', main='lambda')
vioplot(apply(lambda.post, 2, median), add=T)
abline(h=D.true, col=2)

plot(1, ylim=c(0, 140), type='n', main='sigma')
vioplot(as.vector(sigma.post), add=T)
abline(h=sigma, col=2)

hist(apply(lambda.post, 2, median) - D.true)
mean(apply(lambda.post, 2, median) - D.true)
