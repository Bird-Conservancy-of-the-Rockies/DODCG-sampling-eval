##############################################
### Simulate sample sizes needed for IMBCR
### Account for density, heterogeneity
### Monitor precision, bias
##############################################

rm(list = ls())
#setwd("/home/RMBO.LOCAL/rstudio05")
setwd("C:/Users/Quresh.Latif/files/projects/power_analyses/DOD")
scripts.loc <- "DODCG-sampling-eval/"
out.loc <- "sim files/"

library(jagsUI)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(stringr)
lognorm.params <- function(mean, var) {
  m <- log(mean^2/sqrt(var+(mean^2)))
  sig <- sqrt(log(var/(mean^2) +1))
  return(list(m=m, sig=sig))
}

D.true <- c(1) #Birds/km2
Het.true <- c(0) #CV around D.true to simulate heterogeneity in density
sample.true <- c(100) #Number of grids to sample
pts.per.grid <- 15 #Number of points in a grid
cutoffs <- c(100)
theta.true <- c(55)
n.det.xtra <- 0 # number of additional detections for estimating p

n.iter <- 20 #number of iterations

### Create scenarios ###
output.table <- expand.grid(Density=D.true,
                          Heterogeneity=Het.true,
                          N.samp=sample.true,
                          theta = theta.true,
                          Iteration=1:n.iter,
                          D.est=NA,
                          D.SD=NA,
                          D.CV=NA,
                          D.bias=NA,
                          D.pctbias=NA,
                          D.coverage=NA)


for(d in 1:length(D.true)) { #Density loop
  D.d <- D.true[d]
  for(h in 1:length(Het.true)) { #Heterogeneity loop
    Het.h <- Het.true[h]
    for(n in 1:length(sample.true)) { #Sample size loop
      samp.n <- sample.true[n]
      for(p in 1:length(theta.true)) {
        theta <- theta.true[p]
        cutoff <- cutoffs[p]
        for(i in 1:n.iter) { #iteration loop
          pt1=proc.time()
          out.row <- which(output.table$Density==D.d &
                             output.table$Heterogeneity==Het.h &
                             output.table$N.samp==samp.n &
                             output.table$theta==theta &
                             output.table$Iteration==i)
          
          #____ Detection parameters _____#
          area.circle <- pi*(cutoff/1000)^2 #area of point count circle in km2
          nG <- 10 #Number of distance classes
          breaks <- seq(0, cutoff, length.out=nG + 1) #Break points for distance classes
          area.band <- pi*breaks[-1]^2 - pi*breaks[-(nG+1)]^2
          area.prop <- area.band/sum(area.band)
          
          int.true <- (theta^2*(1 - exp(-breaks[-1]^2/(2*theta^2))) - theta^2*(1 - exp(-breaks[-(nG+1)]^2/(2*theta^2)))) #integral of half-normal function
          p.true <- 2*pi*int.true/area.band #Complete half-normal calculations
          pi.true <- p.true*area.prop #Correct for availability
          overall.p.true <- sum(pi.true) #Overall detection probability
          #________________________________#
          
          n.pt <- samp.n*pts.per.grid #Number of points
          grid.pt.ind <- rep(1:samp.n,each = pts.per.grid) #Index mapping point to grid
          
          mean.i <- lognorm.params(D.d, (D.d*Het.h)^2)$m
          sd.i <- lognorm.params(D.d, (D.d*Het.h)^2)$sig
          lambda.true <- exp(rnorm(samp.n, mean.i, sd.i))*area.circle #Grid-level mean abundance, given heterogeneity level
          N.true <- rpois(n.pt,lambda.true[grid.pt.ind]) #Sample point-level abundance

          y.sum <- rbinom(n.pt, N.true, overall.p.true) #Sample number of observations at each point
          
          y.bin <- matrix(NA, n.pt, nG) #Create matrix to hold distance class data
          for(j in 1:n.pt) {
            y.bin[j,] <- rmultinom(1, y.sum[j], pi.true)
          } #close j loop
          
          det.xtra <- rmultinom(1, n.det.xtra, pi.true) #sample extra detections
          det.all <- det.xtra + colSums(y.bin) #add point-level detections
          
          for(j in 1:nG) {
            if(j == 1) {dclass <- rep(j, det.all[j])}
            if(j > 1) {dclass <- c(dclass, rep(j, det.all[j]))}
          } #close j loop
          
          ### Initial values ###
          #        Z.init <- y.sum
          #        Z.init[which(Z.init>0)] <- 1
          
          data <- list(n.grid=samp.n, n.pt=n.pt, grid.pt.ind=grid.pt.ind, area.circle=area.circle, nG=nG, breaks=breaks,
                       area.band=area.band, area.prop=area.prop, y.sum=y.sum, y.bin=y.bin, nind=length(dclass), dclass=dclass)
          params <- c('lambda', 'sigma.mn', 'pcap', 'p', 'D.all')
          inits <- list(list('sigma.mn' = 400, 'log.lambda' = lambda.true, 'N.pt' = apply(y.bin, 1, sum)*2+2), # 'N.pt' = y.sum
                        list('sigma.mn' = 400, 'log.lambda' = lambda.true, 'N.pt' = apply(y.bin, 1, sum)*2+2), # 'N.pt' = y.sum
                        list('sigma.mn' = 400, 'log.lambda' = lambda.true, 'N.pt' = apply(y.bin, 1, sum)*2+2)) # 'N.pt' = y.sum
          model.out <- jagsUI(data = data, parameters.to.save = params, inits = inits,
                              model.file = stringr::str_c(scripts.loc, 'sample_size_sim_HNdot_identity.jags'),
                              n.chains = 3, n.iter = 3000, n.burnin = 1000, n.thin = 1, parallel = TRUE, verbose = FALSE)
          # Populate density output row
          output.table$D.est[out.row] <- median(model.out$sims.list$D.all)
          output.table$D.SD[out.row] <- sd(model.out$sims.list$D.all)
          output.table$D.CV[out.row] <- output.table$D.SD[out.row]/(output.table$D.est[out.row]+0.001)
          output.table$D.bias[out.row] <- output.table$D.est[out.row] - output.table$Density[out.row]
          output.table$D.pctbias[out.row] <- output.table$D.bias[out.row]/output.table$Density[out.row] # This is actually a proportion
          output.table$D2.est[out.row] <- (model.out$sims.list$lambda / area.circle) %>% apply(2, median) %>% mean
          output.table$D2.bias[out.row] <- output.table$D2.est[out.row] - output.table$Density[out.row]
          output.table$D2.pctbias[out.row] <- output.table$D2.bias[out.row]/output.table$Density[out.row] # This is actually a proportion
          if(output.table$D.bias[out.row]==0) {output.table$D.pctbias[out.row]<-0}
          output.table$D.coverage[out.row] <- ifelse(quantile(model.out$sims.list$D.all,prob=0.025,type=8) <=
                                                       output.table$Density[out.row] &&
                                                       quantile(model.out$sims.list$D.all,prob=0.975,type=8) >=
                                                       output.table$Density[out.row], 1, 0)
          
          pt2=proc.time()
          pt2-pt1
          
        } #close i loop
      } # close p loop
    } #close n loop
  } #close h loop
} #close d loop

output.table$D.pctbias <- output.table$D.pctbias*100 # Converts proportion to percent bias.
output.table$D2.pctbias <- output.table$D2.pctbias*100 # Converts proportion to percent bias.
#write.csv(output.table, str_c(out.loc, 'IMBCR sample size sim.csv'), row.names = FALSE)


#output.table <- read.csv(str_c(out.loc, 'IMBCR sample size sim.csv'), header = T, stringsAsFactors = F)# If loading from file
sum.table <- output.table %>%
  dplyr::group_by(Density, Heterogeneity, N.samp, theta) %>%
  summarise(
    D.pctbias.median = median(D.pctbias),
    D.bias = median(D.bias),
    D2.pctbias.median = median(D2.pctbias),
    D2.bias = median(D2.bias),
    D.cov = mean(D.coverage),
    D.CV.median = median(D.CV),
    D.CV.LCI = quantile(D.CV, prob = 0.025, type = 8),
    D.CV.UCI = quantile(D.CV, prob = 0.975, type = 8),
    D.RMSE = sqrt(mean(D.bias^2)),
    D.RMSE.pct = sqrt(mean(D.pctbias^2))
  )
