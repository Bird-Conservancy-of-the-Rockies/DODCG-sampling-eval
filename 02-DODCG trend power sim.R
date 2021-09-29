###########################################################
### Simulate sample sizes needed to estimate IMBCR trends
### Account for density, spatial and annual variation
### Monitor precision, bias, RMSE
###########################################################

rm(list = ls())
#setwd("C:\\Users\\Adam.Green\\Documents\\IMBCR\\Power analysis\\IMBCR simulations")
setwd("C:/Users/Quresh.Latif/files/projects/power_analyses/DOD")
scripts.loc <- "DODCG-sampling-eval/"
out.loc <- "sim files/"

library(dplyr)
library(stringr)
library(jagsUI)
library(ggplot2)
library(gridExtra)
lognorm.params <- function(mean, var) {
  m <- log(mean^2/sqrt(var+(mean^2)))
  sig <- sqrt(log(var/(mean^2) +1))
  return(list(m=m, sig=sig))
}
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

alpha.true <- log(c(1, 5, 20)) #Birds/km2
Het.true <- c(0.1, 1) #CV around D.true to simulate heterogeneity in density
sample.true <- c(15, 20, 25, 30) #Number of grids to sample
r.true <- log(c(0.9, 0.97, 1.03, 1.1)) #Trend magnitudes
yrs <- c(16,31) #Number of years of sampling
#yrs <- c(15,30) #Length of time for monitoring
sample.freq <- c(1, 3, 5)
ann.var <- c(0.1, 0.5) #Annual variation around trend line
pts.per.grid <- 11 #Number of points in a grid

n.iter <- 100 #number of iterations

### Detection parameters
cutoff <- 145 #Distance in meters
area.circle <- pi*(cutoff/1000)^2 #area of point count circle in km2
nG <- 10 #Number of distance classes
breaks <- seq(0, cutoff, length.out=nG + 1) #Break points for distance classes
area.band <- pi*breaks[-1]^2 - pi*breaks[-(nG+1)]^2
area.prop <- area.band/sum(area.band)

theta.true <- 50 #True half-normal shape parameter. Gives overall detection of 0.54 with 250m cutoff

int.true <- (theta.true^2*(1 - exp(-breaks[-1]^2/(2*theta.true^2))) -
               theta.true^2*(1 - exp(-breaks[-(nG+1)]^2/(2*theta.true^2)))) #integral of half-normal function
p.true <- 2*pi*int.true/area.band #Complete half-normal calculations
pi.true <- p.true*area.prop #Correct for availability
overall.p.true <- sum(pi.true) #Overall detection probability
n.det.xtra <- 200 # number of additional detections for estimating p

### Create scenarios ###
output.table <- expand.grid(Alpha=alpha.true,
                            Trend=exp(r.true),
                            Heterogeneity=Het.true,
                            Annual.var=ann.var,
                            N.samp=sample.true,
                            Years=yrs,
                            Samp.freq=sample.freq,
                            Iteration=1:n.iter,
                            r.est=NA,
                            r.SD=NA,
                            r.CV=NA,
                            r.bias=NA,
                            r.pctbias=NA,
                            r.coverage=NA,
                            r.f=NA)

for(a in 1:length(alpha.true)) { #Density loop  1
  alpha.a <- alpha.true[a]
  for(h in 1:length(Het.true)) { #Heterogeneity loop 3
    Het.h <- Het.true[h]
    for(v in 1:length(ann.var)) { # 3
      var.v <- ann.var[v]
      for(r in 1:length(r.true)) { # 1
        r.r <- r.true[r]
        for(n in 1:length(sample.true)) { #Sample size loop 3
          samp.n <- sample.true[n]
          for(t in 1:length(yrs)) { #1
            yrs.t <- yrs[t]
            for(f in 1:length(sample.freq)) {
              samp.freq <- sample.freq[f]
              for(i in 1:n.iter) { #iteration loop
                pt1=proc.time()
                out.row <- which(output.table$Alpha==alpha.a &
                                   output.table$Heterogeneity==Het.h &
                                   output.table$N.samp==samp.n &
                                   output.table$Years==yrs.t &
                                   output.table$Trend==exp(r.r) &
                                   output.table$Annual.var==var.v &
                                   output.table$Samp.freq==samp.freq &
                                   output.table$Iteration==i)
                
                n.pt <- samp.n*pts.per.grid #Number of points
                grid.pt.ind <- rep(1:samp.n, each = pts.per.grid) #Index mapping point to grid
                
                lam.mu <- exp(alpha.a + r.r*((1:yrs.t)-1)) #Mean density per year
                mean.i <- lognorm.params(lam.mu, (lam.mu[1]*var.v)^2)$m #lognormal means
                sd.i <- lognorm.params(lam.mu, (lam.mu[1]*var.v)^2)$sig #lognormal SDs
                lam.true <- exp(rnorm(yrs.t, mean.i, sd.i)) #Year-specific abundance, given heterogeneity level
                grid.mn.i <- lognorm.params(lam.true,(lam.true*Het.h)^2)$m #lognormal means for spatial variation
                grid.sd.i <- lognorm.params(lam.true,(lam.true*Het.h)^2)$sig #lognormal SDs for spatial variation
                lambda.true <- matrix(NA,samp.n,yrs.t) #Grid-level means by year
                N.true <- matrix(NA,n.pt,yrs.t) #Point-level abundances by year
                for(k in 1:yrs.t) { #Sample grid-level means and point-level abundance
                  lambda.true[,k] <- exp(rnorm(samp.n, grid.mn.i[k], grid.sd.i[k]))*area.circle
                  N.true[,k] <- rpois(n.pt,lambda.true[grid.pt.ind,k])
                }
                
                y.sum <- matrix(NA,n.pt,yrs.t)
                for(k in 1:yrs.t) {
                  y.sum[,k] <- rbinom(n.pt, N.true[,k], overall.p.true) #Sample number of observations at each point
                }
                y.bin <- array(NA, dim=c(n.pt, yrs.t, nG)) #Create matrix to hold distance class data
                for(k in 1:yrs.t) {
                  for(j in 1:n.pt) {
                    y.bin[j,k,] <- rmultinom(1, N.true[j,k], pi.true)
                  } #close j loop
                }# close k loop
                if(samp.freq>1) {
                  ind.keep <- seq(1, yrs.t, by = samp.freq)
                  y.sum[,-c(ind.keep)] <- NA
                  y.bin[,-c(ind.keep),] <- NA
                }
                
                det.xtra <- rmultinom(1, n.det.xtra, pi.true) #sample extra detections
                det.all <- det.xtra + apply(y.bin,3,sum,na.rm=T) #add point-level detections
                
                for(j in 1:nG) {
                  if(j == 1) {dclass <- rep(j, det.all[j])}
                  if(j > 1) {dclass <- c(dclass, rep(j, det.all[j]))}
                } #close j loop
                
                ### Initial values ###
                #        Z.init <- y.sum
                #        Z.init[which(Z.init>0)] <- 1
                
                data <- list(n.grid=samp.n, n.pt=n.pt, grid.pt.ind=grid.pt.ind,
                             area.circle=area.circle, nG=nG, breaks=breaks,
                             area.band=area.band, area.prop=area.prop, nind=length(dclass),
                             y.sum=y.sum, dclass=dclass, T=yrs.t)
                #params <- c('alpha', 'r','eps.sd','theta.mn', 'pcap', 'p', 'D.all') # Doesn't seem like we need most of this....
                params <- c('r')
                inits <- list(list('theta.mn' = 400, 'N.pt' = y.sum),
                              list('theta.mn' = 400, 'N.pt' = y.sum),
                              list('theta.mn' = 400, 'N.pt' = y.sum))
                #              if(samp.n*yrs.t>=30) {mod.file='trend_sim_HNdot_RE.jags'}
                #              if(samp.n*yrs.t<30) {mod.file='trend_sim_HNdot_noRE.jags'}
                model.out <- jagsUI(data = data, parameters.to.save = params, inits = inits,
                                    model.file = str_c(scripts.loc, 'trend_sim_HNdot_noRE.jags'),
                                    n.chains = 3, n.iter = 1000, n.burnin = 1000*0.5, n.thin = 1,
                                    parallel = TRUE, verbose = FALSE)
                
                # Populate r output row
                output.table$r.est[out.row[1]] <- median(exp(model.out$sims.list$r))
                output.table$r.SD[out.row[1]] <- sd(exp(model.out$sims.list$r))
                output.table$r.CV[out.row[1]] <-
                  (output.table$r.SD[out.row[1]])/output.table$r.est[out.row[1]]
                output.table$r.bias[out.row[1]] <-
                  output.table$r.est[out.row[1]] - output.table$Trend[out.row[1]]
                output.table$r.pctbias[out.row[1]] <-
                  output.table$r.bias[out.row[1]]/output.table$Trend[out.row[1]]
                if(output.table$r.bias[out.row[1]]==0) {output.table$r.pctbias[out.row[1]] <- 0}
                output.table$r.coverage[out.row[1]] <-
                  ifelse(quantile(exp(model.out$sims.list$r),probs=0.025) <=
                           output.table$Trend[out.row[1]] &&
                           quantile(exp(model.out$sims.list$r),probs=0.975) >=
                           output.table$Trend[out.row[1]], 1, 0)
                if(output.table$Trend[out.row[1]]<1) {output.table$r.f[out.row[1]] <-
                  length(which(model.out$sims.list$r<0))/model.out$mcmc.info$n.samples}
                if(output.table$Trend[out.row[1]]>1) {output.table$r.f[out.row[1]] <-
                  length(which(model.out$sims.list$r>0))/model.out$mcmc.info$n.samples}
                
                pt2=proc.time()
                pt2-pt1
                gc()
              } #close i loop
              
              cat('Alpha = ',a,' of ',length(alpha.true), 'Het = ',h,' of ',length(Het.true),
                  ', Annual var = ',v,' of ',length(ann.var),
                  ', Trend = ',r,' of ',length(r.true),', n = ',n,' of ',length(sample.true),
                  ', Frequency = ',f,' of ',length(samp.freq),
                  ', Years = ',t,' of ',length(yrs),'\n',sep='')
            } # close f loop
          } #close t loop
        } #close n loop
        write.csv(output.table, paste(out.loc,'IMBCR trend sim.csv',sep=''), row.names = FALSE)
      } #close r loop
    } #close v loop
  } #close h loop
} #close a loop

output.table$r.pctbias <- output.table$r.pctbias*100
write.csv(output.table, paste(out.loc,'IMBCR trend sim.csv',sep=''), row.names = FALSE)
# If loading from files:
out.table<-read.csv(paste(out.loc,'IMBCR trend sim.csv',sep=''),header=TRUE)

# sum.table <- output.table %>%
#   dplyr::group_by(Alpha, Trend, Heterogeneity, Annual.var, N.samp, Years, Samp.freq) %>%
#   summarise(r.bias.median = median(r.bias),
#             r.pctbias.median = median(r.pctbias),
#             #r.f.median = median(r.f),
#             #r.f.LCI = quantile(r.f, prob = 0.025, type = 8),
#             #r.f.UCI = quantile(r.f, prob = 0.975, type = 8),
#             r.power = (sum(r.f >= 0.975) / n()) * 100,
#             r.coverage = mean(r.coverage) * 100,
#             r.RMSE = sqrt(mean(r.bias^2)), # When first troubleshooting, make sure there are no NAs.
#             r.RMSE.pct = sqrt(mean((r.pctbias/100)^2)) * 100) # When first troubleshooting, make sure there are no NAs.
# write.csv(sum.table, 'IMBCR trend sim summary.csv', row.names = F)
sum.table <- read.csv(str_c(out.loc, 'IMBCR trend sim summary.csv'), header = T, stringsAsFactors = F)


### Plot % r bias
pdf(str_c(out.loc, 'Trend percent bias.pdf'))

offset=seq(-0.4,0.4,length.out=length(sample.true))

for(t in 1:length(yrs)) {
  yrs.t<-yrs[t]
  for(a in 1:length(alpha.true)) {
    alpha.a<-alpha.true[a]
    for(f in 1:length(sample.freq)) {
      s.freq <- sample.freq[f]
      for(r in 1:length(r.true)) {
        r.r<-exp(r.true[r])
        
        sub.tab <- sum.table %>% filter(round(Alpha, digits=4)==round(alpha.a, digits=4) &
                                          Trend==r.r & Years==yrs.t & Samp.freq==s.freq)
        ctr=1
        for(i in 1:length(Het.true)) {
          for(j in 1:length(ann.var)) {
            sub.tab$N.samp[which(sub.tab$Heterogeneity==Het.true[i] &
                                   sub.tab$Annual.var==ann.var[j])] <-
              sub.tab$N.samp[which(sub.tab$Heterogeneity==Het.true[i] &
                                     sub.tab$Annual.var==ann.var[j])] +
              offset[ctr]
            ctr=ctr+1
          }
        }
        sub.tab$Heterogeneity <- as.factor(sub.tab$Heterogeneity)
        sub.tab$Annual.var <- as.factor(sub.tab$Annual.var)
        
        sub.tab <- sub.tab[which(is.na(sub.tab$r.pctbias.median)==FALSE),]
        
        plot.num <- paste('g',r,sep='')
        sf.name <- ifelse(s.freq == 1, 'Yearly for ',
                          paste('Every ', s.freq,' for ',sep=''))
        yr.name <- paste(yrs.t,' yrs',sep='')
        D.name <- paste(exp(alpha.a),'/km2',sep='')
        r.name <- paste(100*(r.r-1),'%',sep='')
        title.name <- paste(str_c(sf.name,yr.name),D.name,r.name,sep=", ")
        eval(parse(text=paste(plot.num,"<-ggplot(sub.tab, aes(x=N.samp, y=r.pctbias.median)) + ",
                              "geom_point(aes(col=Heterogeneity,shape=Annual.var),size=2) + ",
                              "coord_cartesian(ylim=c(min(c(-10, sum.table$r.pctbias.median)),
                                max(c(10, sum.table$r.pctbias.median)))) + ",
                              "geom_hline(yintercept = 0) + ",
                              "theme_bw() + ",
                              "ggtitle(title.name) + ",
                              "labs(x='Number of grids', y='Trend % bias') + ",
                              "theme(axis.title=element_text(size=14),  ",
                              "plot.title=element_text(size=12), ",
                              "axis.title.x=element_text(vjust=00,size=12), ",
                              "axis.title.y=element_text(size=12),  ",
                              "legend.position=c(0.5,0),",
                              "legend.justification=c(0.5,0),",
                              "legend.background=element_blank(),",
                              "legend.key=element_blank(),",
                              "legend.spacing.y=unit(0,'mm'),",
                              "axis.text.x=element_text(size=10, ",
                              "angle = 30, ",
                              "vjust=.5),  ",
                              "axis.text.y=element_text(size=10)) +",
                              "guides(col = guide_legend(ncol = 2,
                              direction='horizontal',
                              keyheight=unit(1,'mm')),
                              shape = guide_legend(ncol = 2,
                              direction='horizontal',
                              keyheight=unit(1,'mm')))",sep="")))
      } #close r loop
      sub.tab<-expand.grid(N.samp=sample.true,Heterogeneity=Het.true,Annual.var=ann.var,r.pctbias.median=NA)
      sub.tab$Heterogeneity<-as.factor(sub.tab$Heterogeneity)
      sub.tab$Annual.var<-as.factor(sub.tab$Annual.var)
      plot.name<-paste('g',length(r.true)+1,sep='')
      eval(parse(text=paste(plot.name, " <- ggplot(sub.tab, aes(x=N.samp, y=r.pctbias.median)) + ",
                            "geom_point(aes(col=Heterogeneity,shape=Annual.var),size=2) + ",
                            "theme_bw()",sep='')))
      plots<-paste("grid.arrange(",paste(paste('g',1:(length(r.true)-1),sep=''),
                                         collapse="+theme(legend.position='hidden'),"),
                   "+theme(legend.position='hidden'),",paste('g',length(r.true),sep=''),
                   ",nrow=2,ncol=2,newpage=TRUE)",sep="")
      eval(parse(text=plots))
    } #close f loop
  } #close a loop
} #close t loop
dev.off()

# ### Plot r bias (Don't really care about this as much.)
# pdf(str_c(out.loc, 'Trend bias.pdf'))
# 
# offset=seq(-1,1,by=0.1)
# 
# for(t in 1:length(yrs)) {
#   yrs.t<-yrs[t]
#   for(a in 1:length(alpha.true)) {
#     alpha.a<-alpha.true[a]
#     for(f in 1:length(sample.freq)) {
#       s.freq <- sample.freq[f]
#       for(r in 1:length(r.true)) {
#         r.r<-exp(r.true[r])
#         
#         sub.tab <- sum.table %>% filter(round(Alpha, digits=4)==round(alpha.a, digits=4) &
#                                           Trend==r.r & Years==yrs.t & Samp.freq==s.freq)
#         ctr=1
#         for(i in 1:length(Het.true)) {
#           for(j in 1:length(ann.var)) {
#             sub.tab$N.samp[intersect(which(sub.tab$Heterogeneity==Het.true[i]),which(sub.tab$Heterogeneity==ann.var[j]))] <- sub.tab$N.samp[intersect(which(sub.tab$Heterogeneity==Het.true[i]),which(sub.tab$Heterogeneity==ann.var[j]))] - offset[ctr]
#             ctr=ctr+1
#           }
#         }
#         sub.tab$Heterogeneity<-as.factor(sub.tab$Heterogeneity)
#         sub.tab$Annual.var<-as.factor(sub.tab$Annual.var)
#         
#         plot.num <- paste('g',r,sep='')
#         sf.name <- ifelse(s.freq == 1, 'Yearly for ',
#                           paste('Every ', s.freq,' for ',sep=''))
#         yr.name <- paste(yrs.t,' yrs',sep='')
#         D.name <- paste(exp(alpha.a),'/km2',sep='')
#         r.name <- paste(100*(r.r-1),'%',sep='')
#         title.name <- paste(str_c(sf.name,yr.name),D.name,r.name,sep=", ")
#         eval(parse(text=paste(plot.num,"<-ggplot(sub.tab, aes(x=N.samp, y=r.bias.median)) + ",
#                               "geom_point(aes(col=Heterogeneity,shape=Annual.var),size=2) + ",
#                               "coord_cartesian(ylim=c(min(sum.table$r.bias.median),
#                                 max(sum.table$r.bias.median))) + ",
#                               "geom_hline(yintercept = 0) + ",
#                               "theme_bw() + ",
#                               "ggtitle(title.name) + ",
#                               "labs(x='Number of grids', y='Trend bias') + ",
#                               "theme(axis.title=element_text(size=14),  ",
#                               "plot.title=element_text(size=12), ",
#                               "axis.title.x=element_text(vjust=00,size=12), ",
#                               "axis.title.y=element_text(size=12),  ",
#                               "legend.position=c(0.5,0),",
#                               "legend.justification=c(0.5,0),",
#                               "legend.background=element_blank(),",
#                               "legend.key=element_blank(),",
#                               "legend.spacing.y=unit(0,'mm'),",
#                               "axis.text.x=element_text(size=10, ",
#                               "angle = 30, ",
#                               "vjust=.5),  ",
#                               "axis.text.y=element_text(size=10)) +",
#                               "guides(col = guide_legend(ncol = 2,
#                               direction='horizontal',
#                               keyheight=unit(1,'mm')),
#                               shape = guide_legend(ncol = 2,
#                               direction='horizontal',
#                               keyheight=unit(1,'mm')))",sep="")))
#       } #close r loop
#       
#       sub.tab<-expand.grid(N.samp=sample.true,Heterogeneity=Het.true,Annual.var=ann.var,r.bias.median=NA)
#       sub.tab$Heterogeneity<-as.factor(sub.tab$Heterogeneity)
#       sub.tab$Annual.var<-as.factor(sub.tab$Annual.var)
#       plot.name<-paste('g',length(r.true)+1,sep='')
#       eval(parse(text=paste(plot.name, " <- ggplot(sub.tab, aes(x=N.samp, y=r.bias.median)) + ",
#                             "geom_point(aes(col=Heterogeneity,shape=Annual.var),size=2) + ",
#                             "theme_bw()",sep='')))
#       plots<-paste("grid.arrange(",paste(paste('g',1:(length(r.true)-1),sep=''),
#                                          collapse="+theme(legend.position='hidden'),"),
#                    "+theme(legend.position='hidden'),",paste('g',length(r.true),sep=''),
#                    ",nrow=2,ncol=2,newpage=TRUE)",sep="")
#       eval(parse(text=plots))
#     } # close f loop
#   } #close a loop
# } #close t loop
# 
# dev.off()

### Plot r coverage
pdf(str_c(out.loc, 'Trend coverage.pdf'))

offset=seq(-0.4,0.4,by=0.1)

for(t in 1:length(yrs)) {
  yrs.t<-yrs[t]
  for(a in 1:length(alpha.true)) {
    alpha.a<-alpha.true[a]
    for(f in 1:length(sample.freq)) {
      s.freq <- sample.freq[f]
      for(r in 1:length(r.true)) {
        r.r<-exp(r.true[r])
        
        sub.tab <- sum.table %>% filter(round(Alpha, digits=4)==round(alpha.a, digits=4) &
                                          Trend==r.r & Years==yrs.t & Samp.freq==s.freq)
        ctr=1
        for(i in 1:length(Het.true)) {
          for(j in 1:length(ann.var)) {
            sub.tab$N.samp[intersect(which(sub.tab$Heterogeneity==Het.true[i]),which(sub.tab$Heterogeneity==ann.var[j]))] <- sub.tab$N.samp[intersect(which(sub.tab$Heterogeneity==Het.true[i]),which(sub.tab$Heterogeneity==ann.var[j]))] - offset[ctr]
            ctr=ctr+1
          }
        }
        sub.tab$Heterogeneity<-as.factor(sub.tab$Heterogeneity)
        sub.tab$Annual.var<-as.factor(sub.tab$Annual.var)
        
        plot.num <- paste('g',r,sep='')
        sf.name <- ifelse(s.freq == 1, 'Yearly for ',
                          paste('Every ', s.freq,' for ',sep=''))
        yr.name <- paste(yrs.t,' yrs',sep='')
        D.name <- paste(exp(alpha.a),'/km2',sep='')
        r.name <- paste(100*(r.r-1),'%',sep='')
        title.name <- paste(str_c(sf.name,yr.name),D.name,r.name,sep=", ")
        eval(parse(text=paste(plot.num,"<-ggplot(sub.tab, aes(x=N.samp, y=r.coverage)) + ",
                              "geom_point(aes(col=Heterogeneity,shape=Annual.var),size=2) + ",
                              "geom_hline(yintercept = 95) + ",
                              "coord_cartesian(ylim=c(-5,100)) + ",
                              "theme_bw() + ",
                              "ggtitle(title.name) + ",
                              "labs(x='Number of grids', y='Trend coverage') + ",
                              "theme(axis.title=element_text(size=14),  ",
                              "plot.title=element_text(size=12), ",
                              "axis.title.x=element_text(vjust=00,size=12), ",
                              "axis.title.y=element_text(size=12),  ",
                              "legend.position=c(0.5,0),",
                              "legend.justification=c(0.5,0),",
                              "legend.background=element_blank(),",
                              "legend.key=element_blank(),",
                              "legend.spacing.y=unit(0,'mm'),",
                              "axis.text.x=element_text(size=10, ",
                              "angle = 30, ",
                              "vjust=.5),  ",
                              "axis.text.y=element_text(size=10)) +",
                              "guides(col = guide_legend(ncol = 2,
                              direction='horizontal',
                              keyheight=unit(1,'mm')),
                              shape = guide_legend(ncol = 2,
                              direction='horizontal',
                              keyheight=unit(1,'mm')))",sep="")))
      } #close r loop
      
      sub.tab<-expand.grid(N.samp=sample.true,Heterogeneity=Het.true,Annual.var=ann.var,r.coverage=NA)
      sub.tab$Heterogeneity<-as.factor(sub.tab$Heterogeneity)
      sub.tab$Annual.var<-as.factor(sub.tab$Annual.var)
      plot.name<-paste('g',length(r.true)+1,sep='')
      eval(parse(text=paste(plot.name, " <- ggplot(sub.tab, aes(x=N.samp, y=r.coverage)) + ",
                            "geom_point(aes(col=Heterogeneity,shape=Annual.var),size=2) + ",
                            "theme_bw()",sep='')))
      plots<-paste("grid.arrange(",paste(paste('g',1:(length(r.true)-1),sep=''),
                                         collapse="+theme(legend.position='hidden'),"),
                   "+theme(legend.position='hidden'),",paste('g',length(r.true),sep=''),
                   ",nrow=2,ncol=2,newpage=TRUE)",sep="")
      eval(parse(text=plots))
    } # close f loop
  } #close a loop
} #close t loop

dev.off()

### Plot power
pdf(str_c(out.loc, 'Trend power.pdf'))

offset=seq(-0.4,0.4,by=0.1)

for(t in 1:length(yrs)) {
  yrs.t<-yrs[t]
  for(a in 1:length(alpha.true)) {
    alpha.a<-alpha.true[a]
    for(f in 1:length(sample.freq)) {
      s.freq <- sample.freq[f]
      for(r in 1:length(r.true)) {
        r.r<-exp(r.true[r])
        
        sub.tab <- sum.table %>% filter(round(Alpha, digits=4)==round(alpha.a, digits=4) &
                                          Trend==r.r & Years==yrs.t & Samp.freq==s.freq)
        ctr=1
        for(i in 1:length(Het.true)) {
          for(j in 1:length(ann.var)) {
            sub.tab$N.samp[intersect(which(sub.tab$Heterogeneity==Het.true[i]),
                                     which(sub.tab$Heterogeneity==ann.var[j]))] <-
              sub.tab$N.samp[intersect(which(sub.tab$Heterogeneity==Het.true[i]),
                                       which(sub.tab$Heterogeneity==ann.var[j]))] - offset[ctr]
            ctr=ctr+1
          }
        }
        sub.tab$Heterogeneity<-as.factor(sub.tab$Heterogeneity)
        sub.tab$Annual.var<-as.factor(sub.tab$Annual.var)
        
        plot.num <- paste('g',r,sep='')
        sf.name <- ifelse(s.freq == 1, 'Yearly for ',
                          paste('Every ', s.freq,' for ',sep=''))
        yr.name <- paste(yrs.t,' yrs',sep='')
        D.name <- paste(exp(alpha.a),'/km2',sep='')
        r.name <- paste(100*(r.r-1),'%',sep='')
        title.name <- paste(str_c(sf.name,yr.name),D.name,r.name,sep=", ")
        eval(parse(text=paste(plot.num,"<-ggplot(sub.tab, aes(x=N.samp, y=r.power)) + ",
                              "geom_point(aes(col=Heterogeneity,shape=Annual.var),size=2) + ",
                              "geom_hline(yintercept = 80) + ",
                              "coord_cartesian(ylim=c(-5,100)) + ",
                              "theme_bw() + ",
                              "ggtitle(title.name) + ",
                              "labs(x='Number of grids', y='Trend power') + ",
                              "theme(axis.title=element_text(size=14),  ",
                              "plot.title=element_text(size=12), ",
                              "axis.title.x=element_text(vjust=00,size=12), ",
                              "axis.title.y=element_text(size=12),  ",
                              "legend.position=c(0.5,0),",
                              "legend.justification=c(0.5,0),",
                              "legend.background=element_blank(),",
                              "legend.key=element_blank(),",
                              "legend.spacing.y=unit(0,'mm'),",
                              "axis.text.x=element_text(size=10, ",
                              "angle = 30, ",
                              "vjust=.5),  ",
                              "axis.text.y=element_text(size=10)) +",
                              "guides(col = guide_legend(ncol = 2,
                              direction='horizontal',
                              keyheight=unit(1,'mm')),
                              shape = guide_legend(ncol = 2,
                              direction='horizontal',
                              keyheight=unit(1,'mm')))",sep="")))
      } #close r loop
      
      sub.tab<-expand.grid(N.samp=sample.true,Heterogeneity=Het.true,Annual.var=ann.var,r.coverage=NA)
      sub.tab$Heterogeneity<-as.factor(sub.tab$Heterogeneity)
      sub.tab$Annual.var<-as.factor(sub.tab$Annual.var)
      plot.name<-paste('g',length(r.true)+1,sep='')
      eval(parse(text=paste(plot.name, " <- ggplot(sub.tab, aes(x=N.samp, y=r.coverage)) + ",
                            "geom_point(aes(col=Heterogeneity,shape=Annual.var),size=2) + ",
                            "theme_bw()",sep='')))
      plots<-paste("grid.arrange(",paste(paste('g',1:(length(r.true)-1),sep=''),
                                         collapse="+theme(legend.position='hidden'),"),
                   "+theme(legend.position='hidden'),",paste('g',length(r.true),sep=''),
                   ",nrow=2,ncol=2,newpage=TRUE)",sep="")
      eval(parse(text=plots))
    } # close f loop
  } #close a loop
} #close t loop

dev.off()

# ### Plot r RMSE (***Less interpretable than percent RMSE***)
# pdf(str_c(out.loc, 'Trend RMSE.pdf'))
# 
# offset=seq(-0.4,0.4,by=0.1)
# 
# for(t in 1:length(yrs)) {
#   yrs.t<-yrs[t]
#   for(a in 1:length(alpha.true)) {
#     alpha.a<-alpha.true[a]
#     for(f in 1:length(sample.freq)) {
#       s.freq <- sample.freq[f]
#       for(r in 1:length(r.true)) {
#         r.r<-exp(r.true[r])
#         
#         sub.tab <- sum.table %>% filter(round(Alpha, digits=4)==round(alpha.a, digits=4) &
#                                           Trend==r.r & Years==yrs.t & Samp.freq==s.freq)
#         ctr=1
#         for(i in 1:length(Het.true)) {
#           for(j in 1:length(ann.var)) {
#             sub.tab$N.samp[intersect(which(sub.tab$Heterogeneity==Het.true[i]),which(sub.tab$Heterogeneity==ann.var[j]))] <- sub.tab$N.samp[intersect(which(sub.tab$Heterogeneity==Het.true[i]),which(sub.tab$Heterogeneity==ann.var[j]))] - offset[ctr]
#             ctr=ctr+1
#           }
#         }
#         sub.tab$Heterogeneity<-as.factor(sub.tab$Heterogeneity)
#         sub.tab$Annual.var<-as.factor(sub.tab$Annual.var)
#         
#         plot.num <- paste('g',r,sep='')
#         sf.name <- ifelse(s.freq == 1, 'Yearly for ',
#                           paste('Every ', s.freq,' for ',sep=''))
#         yr.name <- paste(yrs.t,' yrs',sep='')
#         D.name <- paste(exp(alpha.a),'/km2',sep='')
#         r.name <- paste(100*(r.r-1),'%',sep='')
#         title.name <- paste(str_c(sf.name,yr.name),D.name,r.name,sep=", ")
#         eval(parse(text=paste(plot.num,"<-ggplot(sub.tab, aes(x=N.samp, y=r.RMSE)) + ",
#                               "geom_point(aes(col=Heterogeneity,shape=Annual.var),size=2) + ",
#                               "coord_cartesian(ylim = c(0,
#                                 0.05*(max(sum.table$r.RMSE)-min(sum.table$r.RMSE))+
#                                 max(sum.table$r.RMSE))) + ",
#                               "theme_bw() + ",
#                               "ggtitle(title.name) + ",
#                               "labs(x='Number of grids', y='Trend RMSE') + ",
#                               "theme(axis.title=element_text(size=14),  ",
#                               "plot.title=element_text(size=12), ",
#                               "axis.title.x=element_text(vjust=00,size=12), ",
#                               "axis.title.y=element_text(size=12),  ",
#                               "legend.position=c(0.5,1),",
#                               "legend.justification=c(0.5,1),",
#                               "legend.background=element_blank(),",
#                               "legend.key=element_blank(),",
#                               "legend.spacing.y=unit(0,'mm'),",
#                               "axis.text.x=element_text(size=10, ",
#                               "angle = 30, ",
#                               "vjust=.5),  ",
#                               "axis.text.y=element_text(size=10)) +",
#                               "guides(col = guide_legend(ncol = 2,
#                               direction='horizontal',
#                               keyheight=unit(1,'mm')),
#                               shape = guide_legend(ncol = 2,
#                               direction='horizontal',
#                               keyheight=unit(1,'mm')))",sep="")))
#       } #close r loop
#       
#       sub.tab<-expand.grid(N.samp=sample.true,Heterogeneity=Het.true,Annual.var=ann.var,r.RMSE=NA)
#       sub.tab$Heterogeneity<-as.factor(sub.tab$Heterogeneity)
#       sub.tab$Annual.var<-as.factor(sub.tab$Annual.var)
#       plot.name<-paste('g',length(r.true)+1,sep='')
#       eval(parse(text=paste(plot.name, " <- ggplot(sub.tab, aes(x=N.samp, y=r.RMSE)) + ",
#                             "geom_point(aes(col=Heterogeneity,shape=Annual.var),size=2) + ",
#                             "theme_bw()",sep='')))
#       plots<-paste("grid.arrange(",paste(paste('g',1:(length(r.true)-1),sep=''),
#                                          collapse="+theme(legend.position='hidden'),"),
#                    "+theme(legend.position='hidden'),",paste('g',length(r.true),sep=''),
#                    ",nrow=2,ncol=2,newpage=TRUE)",sep="")
#       eval(parse(text=plots))
#     } # close f loop
#   } #close a loop
# } #close t loop
# 
# dev.off()


### Plot r RMSE %
pdf(str_c(out.loc, 'Trend percent RMSE.pdf'))

offset=seq(-0.4,0.4,by=0.1)

for(t in 1:length(yrs)) {
  yrs.t<-yrs[t]
  for(a in 1:length(alpha.true)) {
    alpha.a<-alpha.true[a]
    for(f in 1:length(sample.freq)) {
      s.freq <- sample.freq[f]
      for(r in 1:length(r.true)) {
        r.r<-exp(r.true[r])
        
        sub.tab <- sum.table %>% filter(round(Alpha, digits=4)==round(alpha.a, digits=4) &
                                          Trend==r.r & Years==yrs.t & Samp.freq==s.freq)
        ctr=1
        for(i in 1:length(Het.true)) {
          for(j in 1:length(ann.var)) {
            sub.tab$N.samp[intersect(which(sub.tab$Heterogeneity==Het.true[i]),which(sub.tab$Heterogeneity==ann.var[j]))] <- sub.tab$N.samp[intersect(which(sub.tab$Heterogeneity==Het.true[i]),which(sub.tab$Heterogeneity==ann.var[j]))] - offset[ctr]
            ctr=ctr+1
          }
        }
        sub.tab$Heterogeneity<-as.factor(sub.tab$Heterogeneity)
        sub.tab$Annual.var<-as.factor(sub.tab$Annual.var)
        
        plot.num <- paste('g',r,sep='')
        sf.name <- ifelse(s.freq == 1, 'Yearly for ',
                          paste('Every ', s.freq,' for ',sep=''))
        yr.name <- paste(yrs.t,' yrs',sep='')
        D.name <- paste(exp(alpha.a),'/km2',sep='')
        r.name <- paste(100*(r.r-1),'%',sep='')
        title.name <- paste(str_c(sf.name,yr.name),D.name,r.name,sep=", ")
        eval(parse(text=paste(plot.num,"<-ggplot(sub.tab, aes(x=N.samp, y=r.RMSE.pct)) + ",
                              "geom_point(aes(col=Heterogeneity,shape=Annual.var),size=2) + ",
                              "coord_cartesian(ylim = c(0, 100)) + ",
                              "theme_bw() + ",
                              "ggtitle(title.name) + ",
                              "labs(x='Number of grids', y='Trend % RMSE ') + ",
                              "theme(axis.title=element_text(size=14),  ",
                              "plot.title=element_text(size=12), ",
                              "axis.title.x=element_text(vjust=00,size=12), ",
                              "axis.title.y=element_text(size=12),  ",
                              "legend.position=c(0.5,1),",
                              "legend.justification=c(0.5,1),",
                              "legend.background=element_blank(),",
                              "legend.key=element_blank(),",
                              "legend.spacing.y=unit(0,'mm'),",
                              "axis.text.x=element_text(size=10, ",
                              "angle = 30, ",
                              "vjust=.5),  ",
                              "axis.text.y=element_text(size=10)) +",
                              "guides(col = guide_legend(ncol = 2,
                              direction='horizontal',
                              keyheight=unit(1,'mm')),
                              shape = guide_legend(ncol = 2,
                              direction='horizontal',
                              keyheight=unit(1,'mm')))",sep="")))
      } #close r loop
      
      sub.tab<-expand.grid(N.samp=sample.true,Heterogeneity=Het.true,Annual.var=ann.var,r.RMSE.pct=NA)
      sub.tab$Heterogeneity<-as.factor(sub.tab$Heterogeneity)
      sub.tab$Annual.var<-as.factor(sub.tab$Annual.var)
      plot.name<-paste('g',length(r.true)+1,sep='')
      eval(parse(text=paste(plot.name, " <- ggplot(sub.tab, aes(x=N.samp, y=r.RMSE.pct)) + ",
                            "geom_point(aes(col=Heterogeneity,shape=Annual.var),size=2) + ",
                            "theme_bw()",sep='')))
      plots<-paste("grid.arrange(",paste(paste('g',1:(length(r.true)-1),sep=''),
                                         collapse="+theme(legend.position='hidden'),"),
                   "+theme(legend.position='hidden'),",paste('g',length(r.true),sep=''),
                   ",nrow=2,ncol=2,newpage=TRUE)",sep="")
      eval(parse(text=plots))
    } # close f loop
  } #close a loop
} #close t loop

dev.off()

