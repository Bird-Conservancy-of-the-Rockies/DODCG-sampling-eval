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

D.true <- c(1, 5, 20) #Birds/km2
Het.true <- c(0.1, 1) #CV around D.true to simulate heterogeneity in density
sample.true <- c(15, 20, 25, 30) #Number of grids to sample
pts.per.grid <- 15 #Number of points in a grid
cutoffs <- c(250, 100)
sigma.true <- c(150, 55)

n.iter <- 500 #number of iterations

### Create scenarios ###
output.table <- expand.grid(Density=D.true,
                          Heterogeneity=Het.true,
                          N.samp=sample.true,
                          sigma = sigma.true,
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
      for(p in 1:length(sigma.true)) {
        sigma <- sigma.true[p]
        cutoff <- cutoffs[p]
        for(i in 1:n.iter) { #iteration loop
          pt1=proc.time()
          out.row <- which(output.table$Density==D.d &
                             output.table$Heterogeneity==Het.h &
                             output.table$N.samp==samp.n &
                             output.table$sigma==sigma &
                             output.table$Iteration==i)
          
          #____ Detection parameters _____#
          area.circle <- pi*(cutoff/1000)^2 #area of point count circle in km2
          nG <- 10 #Number of distance classes
          breaks <- seq(0, cutoff, length.out=nG + 1) #Break points for distance classes
          area.band <- pi*breaks[-1]^2 - pi*breaks[-(nG+1)]^2
          area.prop <- area.band/sum(area.band)
          
          int.true <- (sigma^2*(1 - exp(-breaks[-1]^2/(2*sigma^2))) - sigma^2*(1 - exp(-breaks[-(nG+1)]^2/(2*sigma^2)))) #integral of half-normal function
          p.true <- 2*pi*int.true/area.band #Complete half-normal calculations
          pi.true <- p.true*area.prop #Correct for availability
          overall.p.true <- sum(pi.true) #Overall detection probability
          n.det.xtra <- 200 # number of additional detections for estimating p
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
                       area.band=area.band, area.prop=area.prop, y.sum=y.sum, nind=length(dclass), dclass=dclass)
          params <- c('lambda', 'sigma.mn', 'pcap', 'p', 'D.all')
          inits <- list(list('sigma.mn' = 400, 'log.lambda' = lambda.true, 'N.pt' = N.true), # 'N.pt' = y.sum
                        list('sigma.mn' = 400, 'log.lambda' = lambda.true, 'N.pt' = N.true), # 'N.pt' = y.sum
                        list('sigma.mn' = 400, 'log.lambda' = lambda.true, 'N.pt' = N.true)) # 'N.pt' = y.sum
          model.out <- jagsUI(data = data, parameters.to.save = params, inits = inits,
                              model.file = stringr::str_c(scripts.loc, 'sample_size_sim_HNdot_identity.jags'),
                              n.chains = 3, n.iter = 3000, n.burnin = 1000, n.thin = 1, parallel = TRUE, verbose = FALSE)
          # Populate density output row
          D.est <- (model.out$sims.list$lambda / area.circle) %>% apply(1, mean)
          output.table$D.est[out.row] <- median(D.est)
          output.table$D.SD[out.row] <- sd(D.est)
          output.table$D.CV[out.row] <- output.table$D.SD[out.row]/(output.table$D.est[out.row]) # +0.001
          output.table$D.bias[out.row] <- output.table$D.est[out.row] - output.table$Density[out.row]
          output.table$D.pctbias[out.row] <- output.table$D.bias[out.row]/output.table$Density[out.row] # This is actually a proportion
          if(output.table$D.bias[out.row]==0) {output.table$D.pctbias[out.row]<-0}
          output.table$D.coverage[out.row] <- ifelse(quantile(D.est,prob=0.025,type=8) <=
                                                       output.table$Density[out.row] &&
                                                       quantile(D.est,prob=0.975,type=8) >=
                                                       output.table$Density[out.row], 1, 0)
          
          pt2=proc.time()
          pt2-pt1
          
        } #close i loop
        cat('D = ',d,' of ',length(D.true),', Het =',h,' of ',length(Het.true),', n = ',n,' of ',length(sample.true),', p = ',p,' of ',length(sigma.true),'\n',sep='')
        write.csv(output.table, str_c(out.loc, 'IMBCR sample size sim.csv'), row.names = FALSE)
      } # close p loop
    } #close n loop
  } #close h loop
} #close d loop

output.table$D.pctbias <- output.table$D.pctbias*100 # Converts proportion to percent bias.
write.csv(output.table, str_c(out.loc, 'IMBCR sample size sim.csv'), row.names = FALSE)

output.table <- read.csv(str_c(out.loc, 'IMBCR sample size sim.csv'), header = T, stringsAsFactors = F)# If loading from file
sum.table <- output.table %>%
  dplyr::group_by(Density, Heterogeneity, N.samp, sigma) %>%
  summarise(
    D.pctbias.median = median(D.pctbias),
    D.bias.median = median(D.bias),
    D.cov = mean(D.coverage),
    D.CV.median = median(D.CV),
    D.CV.LCI = quantile(D.CV, prob = 0.025, type = 8),
    D.CV.UCI = quantile(D.CV, prob = 0.975, type = 8),
    D.RMSE = sqrt(mean(D.bias^2)),
    D.RMSE.pct = sqrt(mean(D.pctbias^2))
  )
write.csv(sum.table, str_c(out.loc, 'IMBCR sample size summary.csv'), row.names = FALSE)
sum.table <- read.csv(str_c(out.loc, 'IMBCR sample size summary.csv'), header = T, stringsAsFactors = F)

### Plot % bias
ymin <- min(c(-10, sum.table$D.pctbias.median))
ymax <- max(c(10, sum.table$D.pctbias.median))

sub.tab <- sum.table[which(sum.table$Heterogeneity==Het.true[1] & sum.table$sigma==sigma.true[1]),]
sub.tab$Density <- as.factor(sub.tab$Density)
sub.tab$N.samp[which(sub.tab$Density==D.true[1])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[1])] - 0.5
sub.tab$N.samp[which(sub.tab$Density==D.true[3])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[3])] + 0.5
g.pbias.Het1.sigma150 <- ggplot(sub.tab, aes(x=N.samp, y=D.pctbias.median)) + 
  geom_point(aes(col=Density)) + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  labs(title=paste('Heterogeneity=',Het.true[1],",\nDetection range = ", sigma.true[1], " m",sep=''), x="Number of grid cells", y="Density % bias") + 
  ylim(ymin, ymax) +
  #  theme(legend.position='None') +
  theme(axis.title=element_text(size=32),  # main title
        axis.title.x=element_text(vjust=00,
                                  size=14),  # X axis title
        axis.title.y=element_text(size=14),  # Y axis title
        axis.text.x=element_text(size=14, 
                                 angle = 30,
                                 vjust=.5),  # X axis text
        axis.text.y=element_text(size=14))

sub.tab <- sum.table[which(sum.table$Heterogeneity==Het.true[2] & sum.table$sigma==sigma.true[1]),]
sub.tab$Density <- as.factor(sub.tab$Density)
sub.tab$N.samp[which(sub.tab$Density==D.true[1])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[1])] - 0.5
sub.tab$N.samp[which(sub.tab$Density==D.true[3])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[3])] + 0.5
g.pbias.Het10.sigma150 <- ggplot(sub.tab, aes(x=N.samp, y=D.pctbias.median)) + 
  geom_point(aes(col=Density)) + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  labs(title=paste('Heterogeneity=',Het.true[2],",\nDetection range = ", sigma.true[1], " m",sep=''), x="Number of grid cells", y="Density % bias") + 
  ylim(ymin, ymax) +
  #  theme(legend.position='None') +
  theme(axis.title=element_text(size=32),  # main title
        axis.title.x=element_text(vjust=00,
                                  size=14),  # X axis title
        axis.title.y=element_text(size=14),  # Y axis title
        axis.text.x=element_text(size=14, 
                                 angle = 30,
                                 vjust=.5),  # X axis text
        axis.text.y=element_text(size=14))

sub.tab <- sum.table[which(sum.table$Heterogeneity==Het.true[1] & sum.table$sigma==sigma.true[2]),]
sub.tab$Density <- as.factor(sub.tab$Density)
sub.tab$N.samp[which(sub.tab$Density==D.true[1])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[1])] - 0.5
sub.tab$N.samp[which(sub.tab$Density==D.true[3])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[3])] + 0.5
g.pbias.Het1.sigma55 <- ggplot(sub.tab, aes(x=N.samp, y=D.pctbias.median)) + 
  geom_point(aes(col=Density)) + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  labs(title=paste('Heterogeneity=',Het.true[1],",\nDetection range = ", sigma.true[2], " m",sep=''), x="Number of grid cells", y="Density % bias") + 
  ylim(ymin, ymax) +
  #  theme(legend.position='None') +
  theme(axis.title=element_text(size=32),  # main title
        axis.title.x=element_text(vjust=00,
                                  size=14),  # X axis title
        axis.title.y=element_text(size=14),  # Y axis title
        axis.text.x=element_text(size=14, 
                                 angle = 30,
                                 vjust=.5),  # X axis text
        axis.text.y=element_text(size=14))

sub.tab <- sum.table[which(sum.table$Heterogeneity==Het.true[2] & sum.table$sigma==sigma.true[2]),]
sub.tab$Density <- as.factor(sub.tab$Density)
sub.tab$N.samp[which(sub.tab$Density==D.true[1])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[1])] - 0.5
sub.tab$N.samp[which(sub.tab$Density==D.true[3])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[3])] + 0.5
g.pbias.Het10.sigma55 <- ggplot(sub.tab, aes(x=N.samp, y=D.pctbias.median)) + 
  geom_point(aes(col=Density)) + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  labs(title=paste('Heterogeneity=',Het.true[2],",\nDetection range = ", sigma.true[2], " m",sep=''), x="Number of grid cells", y="Density % bias") + 
  ylim(ymin, ymax) +
  #  theme(legend.position='None') +
  theme(axis.title=element_text(size=32),  # main title
        axis.title.x=element_text(vjust=00,
                                  size=14),  # X axis title
        axis.title.y=element_text(size=14),  # Y axis title
        axis.text.x=element_text(size=14, 
                                 angle = 30,
                                 vjust=.5),  # X axis text
        axis.text.y=element_text(size=14))

### Plot coverage
sub.tab <- sum.table[which(sum.table$Heterogeneity==Het.true[1] & sum.table$sigma==sigma.true[1]),]
sub.tab$Density <- as.factor(sub.tab$Density)
sub.tab$N.samp[which(sub.tab$Density==D.true[1])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[1])] - 0.5
sub.tab$N.samp[which(sub.tab$Density==D.true[3])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[3])] + 0.5
g.cov.Het1.sigma150 <- ggplot(sub.tab, aes(x=N.samp, y=D.cov)) + 
  geom_point(aes(col=Density)) + 
  geom_hline(yintercept = 0.95) + 
  ylim(0,1) + 
  theme_bw() + 
  labs(title=paste('Heterogeneity=',Het.true[1],",\nDetection range = ", sigma.true[1], " m",sep=''), x="Number of grid cells", y="Density coverage") + 
  #  theme(legend.position='None') +
  theme(axis.title=element_text(size=32),  # main title
        axis.title.x=element_text(vjust=00,
                                  size=14),  # X axis title
        axis.title.y=element_text(size=14),  # Y axis title
        axis.text.x=element_text(size=14, 
                                 angle = 30,
                                 vjust=.5),  # X axis text
        axis.text.y=element_text(size=14))


sub.tab <- sum.table[which(sum.table$Heterogeneity==Het.true[2] & sum.table$sigma==sigma.true[1]),]
sub.tab$Density <- as.factor(sub.tab$Density)
sub.tab$N.samp[which(sub.tab$Density==D.true[1])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[1])] - 0.5
sub.tab$N.samp[which(sub.tab$Density==D.true[3])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[3])] + 0.5
g.cov.Het10.sigma150 <- ggplot(sub.tab, aes(x=N.samp, y=D.cov)) + 
  geom_point(aes(col=Density)) + 
  geom_hline(yintercept = 0.95) + 
  ylim(0,1) + 
  theme_bw() + 
  labs(title=paste('Heterogeneity=',Het.true[2],",\nDetection range = ", sigma.true[1], " m",sep=''), x="Number of grid cells", y="Density coverage") + 
  #  theme(legend.position='None') +
  theme(axis.title=element_text(size=32),  # main title
        axis.title.x=element_text(vjust=00,
                                  size=14),  # X axis title
        axis.title.y=element_text(size=14),  # Y axis title
        axis.text.x=element_text(size=14, 
                                 angle = 30,
                                 vjust=.5),  # X axis text
        axis.text.y=element_text(size=14))

sub.tab <- sum.table[which(sum.table$Heterogeneity==Het.true[1] & sum.table$sigma==sigma.true[2]),]
sub.tab$Density <- as.factor(sub.tab$Density)
sub.tab$N.samp[which(sub.tab$Density==D.true[1])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[1])] - 0.5
sub.tab$N.samp[which(sub.tab$Density==D.true[3])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[3])] + 0.5
g.cov.Het1.sigma55 <- ggplot(sub.tab, aes(x=N.samp, y=D.cov)) + 
  geom_point(aes(col=Density)) + 
  geom_hline(yintercept = 0.95) + 
  ylim(0,1) + 
  theme_bw() + 
  labs(title=paste('Heterogeneity=',Het.true[1],",\nDetection range = ", sigma.true[2], " m",sep=''), x="Number of grid cells", y="Density coverage") + 
  #  theme(legend.position='None') +
  theme(axis.title=element_text(size=32),  # main title
        axis.title.x=element_text(vjust=00,
                                  size=14),  # X axis title
        axis.title.y=element_text(size=14),  # Y axis title
        axis.text.x=element_text(size=14, 
                                 angle = 30,
                                 vjust=.5),  # X axis text
        axis.text.y=element_text(size=14))


sub.tab <- sum.table[which(sum.table$Heterogeneity==Het.true[2] & sum.table$sigma==sigma.true[2]),]
sub.tab$Density <- as.factor(sub.tab$Density)
sub.tab$N.samp[which(sub.tab$Density==D.true[1])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[1])] - 0.5
sub.tab$N.samp[which(sub.tab$Density==D.true[3])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[3])] + 0.5
g.cov.Het10.sigma55 <- ggplot(sub.tab, aes(x=N.samp, y=D.cov)) + 
  geom_point(aes(col=Density)) + 
  geom_hline(yintercept = 0.95) + 
  ylim(0,1) + 
  theme_bw() + 
  labs(title=paste('Heterogeneity=',Het.true[2],",\nDetection range = ", sigma.true[2], " m",sep=''), x="Number of grid cells", y="Density coverage") + 
  #  theme(legend.position='None') +
  theme(axis.title=element_text(size=32),  # main title
        axis.title.x=element_text(vjust=00,
                                  size=14),  # X axis title
        axis.title.y=element_text(size=14),  # Y axis title
        axis.text.x=element_text(size=14, 
                                 angle = 30,
                                 vjust=.5),  # X axis text
        axis.text.y=element_text(size=14))


### Plot CV
ymin <- min(sum.table$D.CV.LCI)
ymax <- max(sum.table$D.CV.UCI)

sub.tab <- sum.table[which(sum.table$Heterogeneity==Het.true[1] & sum.table$sigma==sigma.true[1]),]
sub.tab$Density <- as.factor(sub.tab$Density)
sub.tab$N.samp[which(sub.tab$Density==D.true[1])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[1])] - 0.5
sub.tab$N.samp[which(sub.tab$Density==D.true[3])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[3])] + 0.5
g.CV.Het1.sigma150 <- ggplot(sub.tab, aes(x=N.samp, y=D.CV.median)) + 
  geom_point(aes(col=Density)) + 
  geom_errorbar(aes(x=N.samp,ymin=D.CV.LCI,ymax=D.CV.UCI,width=0)) +
  coord_cartesian(ylim=c(0,1)) + 
  theme_bw() + 
  labs(title=paste('Heterogeneity=',Het.true[1],",\nDetection range = ", sigma.true[1], " m",sep=''), x="Number of grid cells", y="Density CV") + 
  ylim(ymin, ymax) +
  #  theme(legend.position='None') +
  theme(axis.title=element_text(size=32),  # main title
        axis.title.x=element_text(vjust=00,
                                  size=14),  # X axis title
        axis.title.y=element_text(size=14),  # Y axis title
        axis.text.x=element_text(size=14, 
                                 angle = 30,
                                 vjust=.5),  # X axis text
        axis.text.y=element_text(size=14))

sub.tab <- sum.table[which(sum.table$Heterogeneity==Het.true[2] & sum.table$sigma==sigma.true[1]),]
sub.tab$Density <- as.factor(sub.tab$Density)
sub.tab$N.samp[which(sub.tab$Density==D.true[1])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[1])] - 0.5
sub.tab$N.samp[which(sub.tab$Density==D.true[3])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[3])] + 0.5
g.CV.Het10.sigma150 <- ggplot(sub.tab, aes(x=N.samp, y=D.CV.median)) + 
  geom_point(aes(col=Density)) + 
  geom_errorbar(aes(x=N.samp,ymin=D.CV.LCI,ymax=D.CV.UCI,width=0)) +
  coord_cartesian(ylim=c(0,1)) + 
  theme_bw() + 
  labs(title=paste('Heterogeneity=',Het.true[2],",\nDetection range = ", sigma.true[1], " m",sep=''), x="Number of grid cells", y="Density CV") + 
  ylim(ymin, ymax) +
  #  theme(legend.position='None') +
  theme(axis.title=element_text(size=32),  # main title
        axis.title.x=element_text(vjust=00,
                                  size=14),  # X axis title
        axis.title.y=element_text(size=14),  # Y axis title
        axis.text.x=element_text(size=14, 
                                 angle = 30,
                                 vjust=.5),  # X axis text
        axis.text.y=element_text(size=14))

sub.tab <- sum.table[which(sum.table$Heterogeneity==Het.true[1] & sum.table$sigma==sigma.true[2]),]
sub.tab$Density <- as.factor(sub.tab$Density)
sub.tab$N.samp[which(sub.tab$Density==D.true[1])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[1])] - 0.5
sub.tab$N.samp[which(sub.tab$Density==D.true[3])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[3])] + 0.5
g.CV.Het1.sigma55 <- ggplot(sub.tab, aes(x=N.samp, y=D.CV.median)) + 
  geom_point(aes(col=Density)) + 
  geom_errorbar(aes(x=N.samp,ymin=D.CV.LCI,ymax=D.CV.UCI,width=0)) +
  coord_cartesian(ylim=c(0,1)) + 
  theme_bw() + 
  labs(title=paste('Heterogeneity=',Het.true[1],",\nDetection range = ", sigma.true[2], " m",sep=''), x="Number of grid cells", y="Density CV") + 
  ylim(ymin, ymax) +
  #  theme(legend.position='None') +
  theme(axis.title=element_text(size=32),  # main title
        axis.title.x=element_text(vjust=00,
                                  size=14),  # X axis title
        axis.title.y=element_text(size=14),  # Y axis title
        axis.text.x=element_text(size=14, 
                                 angle = 30,
                                 vjust=.5),  # X axis text
        axis.text.y=element_text(size=14))

sub.tab <- sum.table[which(sum.table$Heterogeneity==Het.true[2] & sum.table$sigma==sigma.true[2]),]
sub.tab$Density <- as.factor(sub.tab$Density)
sub.tab$N.samp[which(sub.tab$Density==D.true[1])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[1])] - 0.5
sub.tab$N.samp[which(sub.tab$Density==D.true[3])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[3])] + 0.5
g.CV.Het10.sigma55 <- ggplot(sub.tab, aes(x=N.samp, y=D.CV.median)) + 
  geom_point(aes(col=Density)) + 
  geom_errorbar(aes(x=N.samp,ymin=D.CV.LCI,ymax=D.CV.UCI,width=0)) +
  coord_cartesian(ylim=c(0,1)) + 
  theme_bw() + 
  labs(title=paste('Heterogeneity=',Het.true[2],",\nDetection range = ", sigma.true[2], " m",sep=''), x="Number of grid cells", y="Density CV") + 
  ylim(ymin, ymax) +
  #  theme(legend.position='None') +
  theme(axis.title=element_text(size=32),  # main title
        axis.title.x=element_text(vjust=00,
                                  size=14),  # X axis title
        axis.title.y=element_text(size=14),  # Y axis title
        axis.text.x=element_text(size=14, 
                                 angle = 30,
                                 vjust=.5),  # X axis text
        axis.text.y=element_text(size=14))


# ### Plot RMSE
# ymin <- min(sum.table$D.RMSE)
# ymax <- max(sum.table$D.RMSE)
# 
# sub.tab <- sum.table[which(sum.table$Heterogeneity==Het.true[1] & sum.table$sigma==sigma.true[1]),]
# sub.tab$Density <- as.factor(sub.tab$Density)
# sub.tab$N.samp[which(sub.tab$Density==D.true[1])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[1])] - 0.5
# sub.tab$N.samp[which(sub.tab$Density==D.true[3])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[3])] + 0.5
# g.RMSE.Het1.sigma150 <- ggplot(sub.tab, aes(x=N.samp, y=D.RMSE)) + 
#   geom_point(aes(col=Density)) + 
#   theme_bw() + 
#   labs(title=paste('Heterogeneity=',Het.true[1],",\nDetection range = ", sigma.true[1], " m",sep=''), x="Number of grid cells", y="RMSE") + 
#   ylim(ymin, ymax) +
#   #  theme(legend.position='None') +
#   theme(axis.title=element_text(size=32),  # main title
#         axis.title.x=element_text(vjust=00,
#                                   size=14),  # X axis title
#         axis.title.y=element_text(size=14),  # Y axis title
#         axis.text.x=element_text(size=14, 
#                                  angle = 30,
#                                  vjust=.5),  # X axis text
#         axis.text.y=element_text(size=14))
# 
# 
# sub.tab <- sum.table[which(sum.table$Heterogeneity==Het.true[2] & sum.table$sigma==sigma.true[1]),]
# sub.tab$Density <- as.factor(sub.tab$Density)
# sub.tab$N.samp[which(sub.tab$Density==D.true[1])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[1])] - 0.5
# sub.tab$N.samp[which(sub.tab$Density==D.true[3])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[3])] + 0.5
# g.RMSE.Het10.sigma150 <- ggplot(sub.tab, aes(x=N.samp, y=D.RMSE)) + 
#   geom_point(aes(col=Density)) + 
#   theme_bw() + 
#   labs(title=paste('Heterogeneity=',Het.true[2],",\nDetection range = ", sigma.true[1], " m",sep=''), x="Number of grid cells", y="RMSE") + 
#   ylim(ymin, ymax) +
#   #  theme(legend.position='None') +
#   theme(axis.title=element_text(size=32),  # main title
#         axis.title.x=element_text(vjust=00,
#                                   size=14),  # X axis title
#         axis.title.y=element_text(size=14),  # Y axis title
#         axis.text.x=element_text(size=14, 
#                                  angle = 30,
#                                  vjust=.5),  # X axis text
#         axis.text.y=element_text(size=14))
# 
# 
# sub.tab <- sum.table[which(sum.table$Heterogeneity==Het.true[1] & sum.table$sigma==sigma.true[2]),]
# sub.tab$Density <- as.factor(sub.tab$Density)
# sub.tab$N.samp[which(sub.tab$Density==D.true[1])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[1])] - 0.5
# sub.tab$N.samp[which(sub.tab$Density==D.true[3])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[3])] + 0.5
# g.RMSE.Het1.sigma55 <- ggplot(sub.tab, aes(x=N.samp, y=D.RMSE)) + 
#   geom_point(aes(col=Density)) + 
#   theme_bw() + 
#   labs(title=paste('Heterogeneity=',Het.true[1],",\nDetection range = ", sigma.true[2], " m",sep=''), x="Number of grid cells", y="RMSE") + 
#   ylim(ymin, ymax) +
#   #  theme(legend.position='None') +
#   theme(axis.title=element_text(size=32),  # main title
#         axis.title.x=element_text(vjust=00,
#                                   size=14),  # X axis title
#         axis.title.y=element_text(size=14),  # Y axis title
#         axis.text.x=element_text(size=14, 
#                                  angle = 30,
#                                  vjust=.5),  # X axis text
#         axis.text.y=element_text(size=14))
# 
# 
# sub.tab <- sum.table[which(sum.table$Heterogeneity==Het.true[2] & sum.table$sigma==sigma.true[2]),]
# sub.tab$Density <- as.factor(sub.tab$Density)
# sub.tab$N.samp[which(sub.tab$Density==D.true[1])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[1])] - 0.5
# sub.tab$N.samp[which(sub.tab$Density==D.true[3])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[3])] + 0.5
# g.RMSE.Het10.sigma55 <- ggplot(sub.tab, aes(x=N.samp, y=D.RMSE)) + 
#   geom_point(aes(col=Density)) + 
#   theme_bw() + 
#   labs(title=paste('Heterogeneity=',Het.true[2],",\nDetection range = ", sigma.true[2], " m",sep=''), x="Number of grid cells", y="RMSE") + 
#   ylim(ymin, ymax) +
#   #  theme(legend.position='None') +
#   theme(axis.title=element_text(size=32),  # main title
#         axis.title.x=element_text(vjust=00,
#                                   size=14),  # X axis title
#         axis.title.y=element_text(size=14),  # Y axis title
#         axis.text.x=element_text(size=14, 
#                                  angle = 30,
#                                  vjust=.5),  # X axis text
#         axis.text.y=element_text(size=14))


### Plot RMSE.pct
ymin <- min(sum.table$D.RMSE.pct)
ymax <- max(sum.table$D.RMSE.pct)

sub.tab <- sum.table[which(sum.table$Heterogeneity==Het.true[1] & sum.table$sigma==sigma.true[1]),]
sub.tab$Density <- as.factor(sub.tab$Density)
sub.tab$N.samp[which(sub.tab$Density==D.true[1])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[1])] - 0.5
sub.tab$N.samp[which(sub.tab$Density==D.true[3])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[3])] + 0.5
g.pctRMSE.Het1.sigma150 <- ggplot(sub.tab, aes(x=N.samp, y=D.RMSE.pct)) + 
  geom_point(aes(col=Density)) + 
  theme_bw() + 
  labs(title=paste('Heterogeneity=',Het.true[1],",\nDetection range = ", sigma.true[1], " m",sep=''), x="Number of grid cells", y="RMSE (%)") + 
  ylim(ymin, ymax) +
  #  theme(legend.position='None') +
  theme(axis.title=element_text(size=32),  # main title
        axis.title.x=element_text(vjust=00,
                                  size=14),  # X axis title
        axis.title.y=element_text(size=14),  # Y axis title
        axis.text.x=element_text(size=14, 
                                 angle = 30,
                                 vjust=.5),  # X axis text
        axis.text.y=element_text(size=14))


sub.tab <- sum.table[which(sum.table$Heterogeneity==Het.true[2] & sum.table$sigma==sigma.true[1]),]
sub.tab$Density <- as.factor(sub.tab$Density)
sub.tab$N.samp[which(sub.tab$Density==D.true[1])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[1])] - 0.5
sub.tab$N.samp[which(sub.tab$Density==D.true[3])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[3])] + 0.5
g.pctRMSE.Het10.sigma150 <- ggplot(sub.tab, aes(x=N.samp, y=D.RMSE.pct)) + 
  geom_point(aes(col=Density)) + 
  theme_bw() + 
  labs(title=paste('Heterogeneity=',Het.true[2],",\nDetection range = ", sigma.true[1], " m",sep=''), x="Number of grid cells", y="RMSE (%)") + 
  ylim(ymin, ymax) +
  #  theme(legend.position='None') +
  theme(axis.title=element_text(size=32),  # main title
        axis.title.x=element_text(vjust=00,
                                  size=14),  # X axis title
        axis.title.y=element_text(size=14),  # Y axis title
        axis.text.x=element_text(size=14, 
                                 angle = 30,
                                 vjust=.5),  # X axis text
        axis.text.y=element_text(size=14))

sub.tab <- sum.table[which(sum.table$Heterogeneity==Het.true[1] & sum.table$sigma==sigma.true[2]),]
sub.tab$Density <- as.factor(sub.tab$Density)
sub.tab$N.samp[which(sub.tab$Density==D.true[1])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[1])] - 0.5
sub.tab$N.samp[which(sub.tab$Density==D.true[3])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[3])] + 0.5
g.pctRMSE.Het1.sigma55 <- ggplot(sub.tab, aes(x=N.samp, y=D.RMSE.pct)) + 
  geom_point(aes(col=Density)) + 
  theme_bw() + 
  labs(title=paste('Heterogeneity=',Het.true[1],",\nDetection range = ", sigma.true[2], " m",sep=''), x="Number of grid cells", y="RMSE (%)") + 
  ylim(ymin, ymax) +
  #  theme(legend.position='None') +
  theme(axis.title=element_text(size=32),  # main title
        axis.title.x=element_text(vjust=00,
                                  size=14),  # X axis title
        axis.title.y=element_text(size=14),  # Y axis title
        axis.text.x=element_text(size=14, 
                                 angle = 30,
                                 vjust=.5),  # X axis text
        axis.text.y=element_text(size=14))


sub.tab <- sum.table[which(sum.table$Heterogeneity==Het.true[2] & sum.table$sigma==sigma.true[2]),]
sub.tab$Density <- as.factor(sub.tab$Density)
sub.tab$N.samp[which(sub.tab$Density==D.true[1])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[1])] - 0.5
sub.tab$N.samp[which(sub.tab$Density==D.true[3])] <- sub.tab$N.samp[which(sub.tab$Density==D.true[3])] + 0.5
g.pctRMSE.Het10.sigma55 <- ggplot(sub.tab, aes(x=N.samp, y=D.RMSE.pct)) + 
  geom_point(aes(col=Density)) + 
  theme_bw() + 
  labs(title=paste('Heterogeneity=',Het.true[2],",\nDetection range = ", sigma.true[2], " m",sep=''), x="Number of grid cells", y="RMSE (%)") + 
  ylim(ymin, ymax) +
  #  theme(legend.position='None') +
  theme(axis.title=element_text(size=32),  # main title
        axis.title.x=element_text(vjust=00,
                                  size=14),  # X axis title
        axis.title.y=element_text(size=14),  # Y axis title
        axis.text.x=element_text(size=14, 
                                 angle = 30,
                                 vjust=.5),  # X axis text
        axis.text.y=element_text(size=14))


pdf(str_c(out.loc, 'Density sim summary plots.pdf'))
grid.arrange(g.pbias.Het1.sigma150, g.pbias.Het10.sigma150, g.pbias.Het1.sigma55, g.pbias.Het10.sigma55, nrow=2, ncol=2)
grid.arrange(g.cov.Het1.sigma150, g.cov.Het10.sigma150, g.cov.Het1.sigma55, g.cov.Het10.sigma55, nrow=2, ncol=2)
grid.arrange(g.CV.Het1.sigma150, g.CV.Het10.sigma150, g.CV.Het1.sigma55, g.CV.Het10.sigma55, nrow=2, ncol=2)
#grid.arrange(g.RMSE.Het1.sigma150, g.RMSE.Het10.sigma150, g.RMSE.Het1.sigma55, g.RMSE.Het10.sigma55, nrow=2, ncol=2)
grid.arrange(g.pctRMSE.Het1.sigma150, g.pctRMSE.Het10.sigma150, g.pctRMSE.Het1.sigma55, g.pctRMSE.Het10.sigma55, nrow=2, ncol=2)
dev.off()
