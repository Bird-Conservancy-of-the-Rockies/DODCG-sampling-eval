library(stringr)
library(dplyr)

# Compile IMBCR Bayesian estimates for ADC
rm(list=ls())
setwd("/home/rstudio04@RMBO.local/IMBCR_density")

sp.list <- str_sub(list.files()[which(str_sub(list.files(), -3, -1) == "zip")], 1, 4)

out.tab <- data.frame(species = sp.list, sigma.md = NA, sigma.lo = NA, sigma.hi = NA,
                      stringsAsFactors = F)

for(i in 1:length(sp.list)) {
  
  species=sp.list[i]
  zip.file=paste(species,'.zip',sep='')
  file.lst <- unzip(zip.file,unzip='/usr/bin/unzip',list=TRUE)
  data.file <- paste(species,' IMBCR detection output 2008_2020 trunc_90 bins_10.rdata',sep='')
  unzip(zipfile = zip.file,files = data.file)
  load(data.file)
  file.remove(data.file)
  out.tab[i, c("sigma.md", "sigma.lo", "sigma.hi")] <- det.output$hn.dot["theta.mn", c("50%", "2.5%", "97.5%")] # Need to replace theta with sigma in future analyses.
}

write.csv(out.tab, 'IMBCR spp det param estimates.csv',row.names=F)
