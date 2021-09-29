setwd("C:/Users/Quresh.Latif/files/projects/power_analyses/DOD")

library(tidyverse)
library(BCRDataAPI)
strata <- "WY-DOD-CG"
dat_species <- read.csv("Species of interest_habitat.csv", header = T, stringsAsFactors = F) %>%
  select(Common.Name:Habitat)
spp <- dat_species$BirdCode

# Check Birds of the World Species accounts #
hab <- dat_species$Habitat
names(hab) <- dat_species$BirdCode
hab["AMGO"] <- "multiple"
hab["BLGR"] <- "multiple"
hab["BGGN"] <- "multiple"
hab["CONI"] <- "multiple"
hab["SPTO"] <- "shrubland"
dat_species$Habitat <- hab
rm(hab)

# Average density estimates #
dat_density <- read.csv("Camp Guernsey density estimates.csv", header = T, stringsAsFactors = F)
dat_species <- dat_species %>%
  dplyr::left_join(
    dat_density %>%
      filter(Year == 2013) %>%
      select(Species, TotalD, CV_D, n) %>%
      rename(D_2013 = TotalD, CV_2013 = CV_D, n_2013 = n),
    by = c("BirdCode" = "Species")
  ) %>%
  dplyr::left_join(
    dat_density %>%
      filter(Year == 2015) %>%
      select(Species, TotalD, CV_D, n) %>%
      rename(D_2015 = TotalD, CV_2015 = CV_D, n_2015 = n),
    by = c("BirdCode" = "Species")
  ) %>%
  dplyr::left_join(
    dat_density %>%
      filter(Year == 2020) %>%
      select(Species, TotalD, CV_D, n) %>%
      rename(D_2020 = TotalD, CV_2020 = CV_D, n_2020 = n),
    by = c("BirdCode" = "Species")
  )

# Trend estimates #
dat_trend <- read.csv("WY trend estimates through 2020.csv", header = T, stringsAsFactors = F)


# Get habitat specialization indices #
dat_hab_special <- read.csv("Species_habs_CampG.csv", header = T, stringsAsFactors = F) %>%
  rename(Spp = X)
dat_species <- dat_species %>%
  dplyr::left_join(dat_hab_special, by = c("BirdCode" = "Spp"))

dat_species <- dat_species %>%
  mutate(., D_mean = rowMeans(select(., D_2013, D_2015, D_2020))) %>%
  mutate(., CV_mean = rowMeans(select(., CV_2013, CV_2015, CV_2020))) %>%
  mutate(., n = rowSums(select(., n_2013, n_2015, n_2020))) %>%
  rowwise() %>% mutate(D_min = min(D_2013, D_2015, D_2020)) %>%
  mutate(CV_max = max(CV_2013, CV_2015, CV_2020)) %>%
  select(Common.Name, Habitat, n, wetland_riparian:Grassland, D_mean, D_min, CV_mean, CV_max, BirdCode) %>%
  arrange(Common.Name)
dat_species <- dat_species %>%
  select(Common.Name:Grassland) %>%
  bind_cols(
    dat_species %>%
      select(D_mean:CV_max) %>%
      mutate_all(function(x) round(x, digits = 2) %>% as.character()) %>%
      mutate_all(function(x) ifelse(is.na(x), "--", x)),
    dat_species %>%
      select(BirdCode)
  ) %>%
  mutate(n = ifelse(is.na(n), 0, n))



write.csv(dat_species, "Species_densities_habitat_summary.csv", row.names = F)
