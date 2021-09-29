setwd("C:/Users/Quresh.Latif/files/projects/power_analyses/DOD")

library(tidyverse)
library(BCRDataAPI)
strata <- "WY-DOD-CG"

BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('analysis.api.bcr.eco')
BCRDataAPI::add_columns(c('TransectNum|str',
                          'Stratum|str',
                          'Point|int',
                          'easting|int',
                          'northing|int',
                          'zone|int'
))

BCRDataAPI::filter_on(c(str_c('Stratum in ', str_c(strata, collapse = ","))))
BCRDataAPI::group_by(c('TransectNum', 'Stratum', 'Point', 'easting', 'northing', 'zone'))
grab <- BCRDataAPI::get_data()

write.csv(grab, "IMBCR_point_coords.csv", row.names = F)
