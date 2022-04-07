#
# Title: Cleaning of ARU presence/absence data
# Created: March 28th, 2022
# Last Updated: April 7th, 2022
# Author: Brandon Allen
# Objective: Clean and filter the WildTrax data for cooccurrence analyses.
# Keywords: Notes, Occurrence summaries, Calling Intensity, Landcover summaries, Climate summaries
#

######### 
# Notes # 
#########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# 1) Each recording represents a single observaton.
# 2) For definitions of the labels used to distinguish the types of ARU data, see the BU_acoustic_recording_analysis_protocol-v10
# 3) There are new labels for Wildtrax, we are going to use the 1m and 3m species recordings, but only the detections in the first minute.
# 4) New lookup for the vegetation and soil classes (v2020) are implemented
# 5) This analysis includes only data included in WildTrax. This means some information from other BU projects may not be included. 
# 6) Data was downloaded from all projects with amphibian detections on December 2, 2021
#
########################
# Occurrence summaries # 
########################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##############################
# Occurrence standardization # 
##############################

# Clear memory
rm(list=ls())
gc()

# Load libraries
library(ggplot2)
library(mefa4)
library(opticut)

# Using the original data sets, create one merged document
column.lookup <- read.csv("data/lookup/column-id_lookup_v2022.csv")

occ.files <- list.files(path = "data/base/species/", full.names = TRUE)
occ.files <- occ.files[!(occ.files %in% c("data/base/species/english_column_definitions.csv",
                                          "data/base/species/APPENDED_REPORT_2022.csv"))]

species.data <- NA
for (x in 1:length(occ.files)) {
  
  # Load file
  temp.data <- read.csv(occ.files[x])
  
  # Standardize columns
  temp.col <- column.lookup$ColID[!(column.lookup$ColID %in% colnames(temp.data))]
  new.data <- data.frame(matrix(data = NA, nrow = nrow(temp.data), ncol = length(temp.col), dimnames = list(NULL, temp.col)))
  temp.data <- cbind(temp.data, new.data)
  temp.data <- temp.data[, column.lookup$ColID]
  
  species.data <- rbind(species.data, temp.data)
  print(x)
  
  rm(temp.col, new.data,temp.data)
  
}

write.csv(species.data, file = "data/processed/species/APPENDED_REPORT_2022.csv", row.names = FALSE)
rm(x, species.data, occ.files, column.lookup)

#####################
# Project Selection # 
#####################

# Load amphibian data 
species.data <- read.csv("data/processed/species/APPENDED_REPORT_2022.csv")

# Keep projects based on the lookup table
project.lookup <- read.csv("data/lookup/project-lookup_v2022.csv")
project.lookup <- project.lookup$Project[project.lookup$Inclusion]
species.data <- species.data[species.data$project %in% project.lookup, ]

#########################
# Survey Standarization # 
#########################

# Standardize the recording intervals as characters
species.data$min0_start <- as.character(species.data$min0_start)
species.data$min1_start <- as.character(species.data$min1_start)
species.data$min2_start <- as.character(species.data$min2_start)

# Remove "Heavy" rain, wind, and noise
species.data <- species.data[!(species.data[, "rain"] %in% "Heavy"), ]
species.data <- species.data[!(species.data[, "wind"] %in% "Heavy"), ]
species.data <- species.data[!(species.data[, "industry_noise"] %in% "Heavy"), ]
species.data <- species.data[!(species.data[, "noise"] %in% "Heavy"), ]

# Keep only Confident and Confirmed Calls
species.data <- species.data[species.data$confidence %in% c("Confident", "Confirmed"), ]

# Take the 1 and 3 minute recording intervals
species.data <- species.data[species.data$method %in% c("1m 1SPM", "1m 2SPM", "3m 1SPM", "3m 2SPM"),]
species.data$maxdur <- as.integer(sapply(strsplit(as.character(species.data$method), "m"), "[[", 1))

# Split the characters of the start minutes to capture the first detection time
f <- function(v) {
  v <- strsplit(as.character(v), ",")
  v <- sapply(v, function(z) if (length(z) < 1) NA else z[1])
  v <- as.numeric(v)
  v
}
species.data$min0_start_num <- f(species.data$min0_start)
species.data$min1_start_num <- f(species.data$min1_start)
species.data$min2_start_num <- f(species.data$min2_start)

species.data$Start <- strptime(paste(species.data$recording_date, species.data$recording_time),  "%Y-%m-%d %H:%M:%S")

# Standardize the common and scientific names for each recording
tx <- nonDuplicated(species.data[,c("species_code", "scientific_name","species_english_name")], species_code, TRUE)
tx$m1 <- tx$species_code

# Match the species lookup table to create the presence/absence information for each station
species.data$Spp <- as.factor(tx$m1)[match(species.data$species_code, tx$species_code)]

###############################
# Survey time standardization # 
###############################

# Create information for the time of year, time of year cuts, time of day, and subset base on the recording time (11:30 - 2:30)
species.data$ToY <- species.data$Start$yday
species.data$ToY2 <- species.data$ToY * species.data$ToY
species.data$Year <- species.data$Start$year + 1900
species.data$site_stn <- paste0(species.data$location, "_", species.data$Year)
species.data$ToYc <- as.integer(cut(species.data$ToY, c(0, 105, 120, 140, 150, 160, 170, 180, 365)))
species.data$replicate <- as.character(species.data$Start)
species.data$replicate <- gsub("[[:punct:]]", "", species.data$replicate)
species.data$replicate <- gsub("[[:space:]]", "", species.data$replicate)
species.data$visit <- paste0("ABMISM::", species.data$site_stn, "::", species.data$replicate)

species.data$ToD <- species.data$Start$hour + species.data$Start$min / 60
species.data$ToDx <- round(species.data$ToD, 0)
species.data$ToDc <- as.factor(ifelse(species.data$ToDx %in% c(0, 1, 2, 24), "Midnight", "Morning"))

# Remove columns that are not used
species.data <- species.data[, !(colnames(species.data) %in% c("min3_voc", "min3_start", "min3_freq_range", "min3_tag_duration",
                                                               "min4_voc", "min4_start", "min4_freq_range", "min4_tag_duration",
                                                               "min5_voc", "min5_start", "min5_freq_range", "min5_tag_duration",
                                                               "min6_voc", "min6_start", "min6_freq_range", "min6_tag_duration",
                                                               "min7_voc", "min7_start", "min7_freq_range", "min7_tag_duration",
                                                               "min8_voc", "min8_start", "min8_freq_range", "min8_tag_duration",
                                                               "min9_voc", "min9_start", "min9_freq_range", "min9_tag_duration"))]

# Store long form
species.long <- species.data

#################################
# Vocalization standardization  # 
#################################

# Standardize the output
species.occ <- species.data[, c("project", "organization", "location", "Year", "site_stn", "latitude", "longitude", "hourly_weather_station_distance", "hourly_temp", "hourly_precipitation_mm", "ToY", "ToY2")]
species.occ <- species.occ[!duplicated(species.occ), ]
colnames(species.occ) <- c("Project", "Organization", "Location", "Year", "Site_Location", "Lat", "Long", "Weather_Station_Distance",
                                  "Hourly_Temp", "Hourly_Precipitation", "ToY", "ToY2")

# Add factor to account for stations that have been visited multiple times
visit.template <- data.frame(Site = names(table(species.occ$Location)), 
                             visits = as.numeric(table(species.occ$Location)))

species.occ["visit"] <- visit.template$visits[match(species.occ$Location, visit.template$Site)]
species.occ["Site_Loc_ToY"] <- paste0(species.occ$Site_Location, "_", species.occ$ToY)
species.occ <- species.occ[!duplicated(species.occ$Site_Loc_ToY), ]

# Create species by site matrix

for (spp in sort(tx$species_code)) { 
            
            site.subset <- species.data[species.data$species_code == spp, ]
            site.list <- unique(paste0(site.subset$site_stn, "_", site.subset$ToY))
            species.occ[[spp]] <- ifelse(species.occ$Site_Loc_ToY %in% site.list, 1, 0)
            print(spp)
            
}

#
# Removal of known incorrect species identifications for amphibians
#

# WETO
species.occ[species.occ$Site_Loc_ToY %in% c("1578-NE_2017_123",
                                                "34-NE_2015_156",
                                                "1362-NE_2016_179",
                                                "8-SW_2018_157"), "WETO"] <- 0

# PLSP
species.occ[species.occ$Site_Loc_ToY %in% c("183-NW_2017_195"), "PLSP"] <- 0

# NLFR
species.occ[species.occ$Site_Loc_ToY %in% c("592-SW_2018_190",
                                                "W479-CL_2018_169"), "NLFR"] <- 0

# CSFR
species.occ[species.occ$Site_Loc_ToY %in% c("OG-EI-750-41-1_2017_148",
                                                "1311-NE_2018_164"), "CSFR"] <- 0

# CATO
species.occ[species.occ$Site_Loc_ToY %in% c("1644-NW_2016_146",
                                                "1645-NE_2016_144",
                                                "497-NW_2016_133",
                                                "555-NW_2017_169",
                                                "588-SW_2017_139",
                                                "653-NE_2017_186"), "CATO"] <- 0

# WOFR
# No corrections

# BCFR
# No corrections

# GPTO
# No corrections

# Save results
species.wide <- species.occ
save(species.long, species.wide, file = paste0("data/processed/species/species-abundance_2022-04-07.Rdata"))

rm(list=ls())
gc()
