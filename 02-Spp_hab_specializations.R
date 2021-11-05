setwd("C:/Users/Quresh.Latif/files/projects/power_analyses/DOD")

library(tidyverse)
library(BCRDataAPI)
dat_species <- read.csv("Species of interest_habitat.csv", header = T, stringsAsFactors = F) %>%
  select(Common.Name:Habitat)
spp <- dat_species$BirdCode
SampDesign <- c("IMBCR", "GRTS")

############################################
# Specialization indices for Camp Guernsey #
############################################

# Get all species detections recorded at Camp G.
strata <- "WY-DOD-CG"

BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('analysis.api.bcr.eco')
BCRDataAPI::add_columns(c('TransectNum|str',
                          'Point|int',
                          'Year|int',
                          'BirdCode|str',
                          'Species|str',
                          'CL_count|int'
))

BCRDataAPI::filter_on(c(str_c('BirdCode in ', str_c(spp, collapse = ",")),
                        str_c('Stratum in ', str_c(strata, collapse = ",")),
                        'ninetynine = 0',
                        'eightyeight = 0',
                        'How <> F',
                        'Sex <> J',
                        'Migrant = 0',
                        'TimePeriod > -1'))
BCRDataAPI::group_by(c('TransectNum', 'Point', 'Year', 'BirdCode', 'Species', 'CL_count'))
grab <- BCRDataAPI::get_data(interpolate_effort = T) %>%
  dplyr::group_by(TransectNum, Point, Year, BirdCode, Species) %>%
  summarize(CL_count = sum(CL_count, na.rm = T)) %>%
  ungroup()

# Get GIS-based hab designations #
dat_hab <- foreign::read.dbf("C:/Users/Quresh.Latif/files/GIS/DOD/Surv_points_GCS83.dbf", as.is = T)
grab <- grab %>%
  dplyr::left_join(
    dat_hab %>% select(TransectNu, Point, HabType) %>%
      rename(TransectNum = TransectNu),
    by = c("TransectNum", "Point")
  )
grab <- grab %>%
  filter(!is.na(HabType))

hab.codes <- list(wetland_riparian = c("wetland_riparian"),
                  Canyon = c("canyon/woodland"),
                  Shrubland = c("shrubland", "shrubland/woodland"),
                  Woodland = c("woodland", "canyon/woodland", "shrubland/woodland"),
                  Grassland = c("grassland"))

rows <- spp
cols <- names(hab.codes)
out <- matrix("", nrow = length(rows), ncol = length(cols),
              dimnames = list(rows, cols))

for(i in 1:nrow(out)) for(j in 1:ncol(out)) {
  if(spp[i] %in% grab$BirdCode) {
    if(sum(grab %>% filter(BirdCode == spp[i]) %>% pull(CL_count)) >= 1) {
      nhab <- grab %>%
        filter(HabType %in% hab.codes[[j]]) %>%
        select(TransectNum, Point, Year) %>%
        distinct() %>% nrow
      nnon <- grab %>%
        filter(!HabType %in% hab.codes[[j]]) %>%
        select(TransectNum, Point, Year) %>%
        distinct() %>% nrow
      
      count_hab <- grab %>%
        filter(HabType %in% hab.codes[[j]] &
                 BirdCode == spp[i]) %>%
        pull(CL_count) %>% sum
      count_non <- grab %>%
        filter(!HabType %in% hab.codes[[j]] &
                 BirdCode == spp[i]) %>%
        pull(CL_count) %>% sum
      
      RA_hab <- count_hab / nhab
      RA_non <- count_non / nnon
      
      out[i, j] <- str_c(round(RA_hab / (RA_hab + RA_non), digits = 2),
                         " (",
                         count_hab,
                         ")")
    }
  } else {
    out[i, j] <- "0 (0)"
  }
}

write.csv(out, "Species_habs_CampG.csv", row.names = T)

## Tabulate grid cells and points by GIS habitat categories ##
rows <- c(names(hab.codes), "Total")
cols <- c("Grid cells", "Points")
out <- matrix(NA, nrow = length(rows), ncol = length(cols),
              dimnames = list(rows, cols))

for(i in 1:length(hab.codes)) {
  out[i, "Grid cells"] <- grab %>%
    select(TransectNum, HabType) %>%
    filter(HabType %in% hab.codes[[i]]) %>%
    pull(TransectNum) %>% unique() %>%
    length()
  out[i, "Points"] <- grab %>%
    select(TransectNum, Point, HabType) %>%
    filter(HabType %in% hab.codes[[i]]) %>%
    select(TransectNum, Point) %>% distinct() %>%
    nrow()
}
out["Total", "Grid cells"] <- grab %>%
  pull(TransectNum) %>% unique() %>%
  length()
out["Total", "Points"] <- grab %>%
  select(TransectNum, Point) %>% distinct() %>%
  nrow()

# Add proportion areas #
veg <- foreign::read.dbf("C:/Users/Quresh.Latif/files/GIS/DOD/Vegetation_communites_utm13.dbf", as.is = T)

out <- cbind(out, rep(NA, nrow(out)))
dimnames(out)[[2]][ncol(out)] <- "Percent area"
out["wetland_riparian", "Percent area"] <- round((sum(veg$Acres[which(veg$HabitatTyp == "wetland_riparian")]) / sum(veg$Acres)) * 100, digits = 1)
out["Canyon", "Percent area"] <- round((sum(veg$Acres[which(veg$HabitatTyp == "canyon/woodland")]) / sum(veg$Acres)) * 100, digits = 1)
out["Shrubland", "Percent area"] <- round((sum(veg$Acres[which(veg$HabitatTyp %in% c("shrubland", "shrubland/woodland"))]) / sum(veg$Acres)) * 100, digits = 1)
out["Woodland", "Percent area"] <- round((sum(veg$Acres[which(veg$HabitatTyp %in% c("woodland", "shrubland/woodland"))]) / sum(veg$Acres)) * 100, digits = 1)
out["Grassland", "Percent area"] <- round((sum(veg$Acres[which(veg$HabitatTyp == "grassland")]) / sum(veg$Acres)) * 100, digits = 1)
out["Total", "Percent area"] <- 100

write.csv(out, "Hab_summary_CG_GIS.csv", row.names = T)

###############################################################
# Cross-reference "vegetation communities" with habitat types #
###############################################################

foreign::read.dbf("C:/Users/Quresh.Latif/files/GIS/DOD/Vegetation_communites_utm13.dbf", as.is = T) %>%
  select(featureNam, HabitatTyp, Acres) %>% distinct() %>%
  write.csv("Veg_communities_habitat_type_crossRef.csv", row.names = F)

#####################################
# Specialization indices for BCR 17 #
#####################################

# hab.codes <- list(Wetland = c("WE", "HW"),
#                   Riparian = "RI",
#                   Shrubland = c("SH", "SA", "DS"),
#                   Woodland = c("PJ", "PP", "MC", "LP", "DW", "OA"),
#                   Grassland = c("GR", "4D"))
# 

# # Get all species detections recorded in BCR 17
# BCRDataAPI::reset_api()
# BCRDataAPI::set_api_server('analysis.api.bcr.eco')
# BCRDataAPI::add_columns(c('TransectNum|str',
#                           'Point|int',
#                           'Year|int',
#                           'BirdCode|str',
#                           'Species|str',
#                           'CL_count|int'
# ))
# 
# BCRDataAPI::filter_on(c(str_c('SelectionMethod in ', str_c(SampDesign, collapse = ",")),
#                         str_c('BirdCode in ', str_c(spp, collapse = ",")),
#                         'bcr = 17',
#                         'ninetynine = 0',
#                         'eightyeight = 0',
#                         'How <> F',
#                         'Sex <> J',
#                         'Migrant = 0',
#                         'TimePeriod > -1'))
# BCRDataAPI::group_by(c('TransectNum', 'Point', 'bcr', 'Year', 'BirdCode', 'Species', 'CL_count'))
# grab <- BCRDataAPI::get_data(interpolate_effort = T) %>%
#   dplyr::group_by(TransectNum, Point, Year, BirdCode, Species) %>%
#   summarize(CL_count = sum(CL_count, na.rm = T)) %>%
#   ungroup()
# 
# # Get primary habitat designations #
# BCRDataAPI::reset_api()
# BCRDataAPI::set_api_server('analysis.api.bcr.eco')
# BCRDataAPI::add_columns(c('TransectNum|str',
#                           'Point|int',
#                           'Year|int',
#                           'primaryHabitat|str',
#                           'HabitatCommonName|str'
# ))
# 
# BCRDataAPI::filter_on(c(str_c('SelectionMethod in ', str_c(SampDesign, collapse = ",")),
#                         'bcr = 17'))
# BCRDataAPI::group_by(c('TransectNum', 'Point', 'Year', 'primaryHabitat', 'HabitatCommonName'))
# grab_phab <- BCRDataAPI::get_data()
# 
# grab <- grab %>%
#   dplyr::left_join(grab_phab, by = c("TransectNum", "Point", "Year")) %>%
#   filter(!is.na(primaryHabitat))
# 
# # spp_d10 <- grab %>%
# #   dplyr::group_by(BirdCode) %>%
# #   summarise(sumCount = sum(CL_count)) %>%
# #   filter(sumCount >= 10) %>%
# #   pull(BirdCode)
# # All focal species were detected >= 10 times.
# 
# rows <- spp
# cols <- names(hab.codes)
# out <- matrix(NA, nrow = length(rows), ncol = length(cols),
#               dimnames = list(rows, cols))
# 
# for(i in 1:nrow(out)) for(j in 1:ncol(out)) {
#   if(spp[i] %in% grab$BirdCode) {
#     nhab <- grab %>%
#       filter(primaryHabitat %in% hab.codes[[j]]) %>%
#       select(TransectNum, Point, Year) %>%
#       distinct() %>% nrow
#     nnon <- grab %>%
#       filter(!primaryHabitat %in% hab.codes[[j]]) %>%
#       select(TransectNum, Point, Year) %>%
#       distinct() %>% nrow
#     
#     count_hab <- grab %>%
#       filter(primaryHabitat %in% hab.codes[[j]] &
#                BirdCode == spp[i]) %>%
#       pull(CL_count) %>% sum
#     count_non <- grab %>%
#       filter(!primaryHabitat %in% hab.codes[[j]] &
#                BirdCode == spp[i]) %>%
#       pull(CL_count) %>% sum
#     
#     RA_hab <- count_hab / nhab
#     RA_non <- count_non / nnon
#     
#     out[i, j] <- RA_hab / (RA_hab + RA_non)
#   }
# }
# 
# write.csv(out, "Species_habs_BCR17.csv", row.names = T)

