setwd("C:/Users/Quresh.Latif/files/projects/power_analyses/DOD")

library(tidyverse)

dat <- foreign::read.dbf("C:/Users/Quresh.Latif/files/GIS/DOD/Surv_points_GCS83.dbf", as.is = T)

sum.habs <- dat %>%
  dplyr::group_by(HabType) %>%
  summarise(Count = n()) %>%
  mutate(Prop = round((Count / sum(Count)) * 100))

write.csv(sum.habs, "Habitat_summary.csv", row.names = F)
