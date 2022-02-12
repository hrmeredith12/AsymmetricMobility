###############################################################################
# Ronan Corgel
# Asymmetry Mobility Project
# Prepare Namibia Data
# February 2022
# Steps:
# 1. Import administrative data
# 2. Import CDR files
# 3. Separate and rewrite dates
# 4. Add in other administrative details
# 5. Save data
###############################################################################

# Remove objects previously stored in R
rm(list = ls())

# Set Directories
setwd('/Users/rcorgel/OneDrive/Projects/asymmetry-project') 

# Load libraries
library('geosphere')  # haverstine distance function
library('tidyverse')  # joining
library('lubridate')  # parsing dates
library('reshape')    # for melting data
library('assertr')    # test characteristics of a dataset
library("tidyr")

# 1. Import administrative data
# The CDRs just have an adm2 ID number - no name or adm1 ID. 
# This dataframe has characteristics of each adm2 unit compiled previously from a number of sources. 
# The information will be merged with the CDRs later. You can ignore the variables on urbanicity for now (urb...)
load("raw/NAM.ID.link.RData")  # This dataframe links the admin2 unit names with the IDs 

# 2. Import CDR files, add in relevant trip info (origin and destination details) 
load("raw/NAM_between_day_mobility_updated_joinedAdm2.RData") # This dataframe has the daily trip counts between each admin2 unit
trip.data <- mydata.joined
trip.data.long <- reshape2::melt(trip.data, id.vars = c("i", "j")) ## To get the dataframe in long form
colnames(trip.data.long) <- c("start.adm2.code", "end.adm2.code", "date", "trip.count")

# 3. Separate and rewrite dates
dates.sep <- strsplit(as.character(trip.data.long$date), "/", perl=TRUE)               # separate date string 
trip.data.long$d <- as.integer(sapply(dates.sep, function(x) x[1]))                    # define day of date
trip.data.long$m <- as.integer(sapply(dates.sep, function(x) x[2]))                    # month of date
trip.data.long$y <- as.integer(sapply(dates.sep, function(x) x[3]))                    # year of date
trip.data.long$date <- as.Date(with(trip.data.long, paste(y, m, d,sep="-")), "%Y-%m-%d")

# 4. Add in other administrative details
trip.data.long <- left_join(trip.data.long, county.file.ordered, by = c("start.adm2.code" = "ID_2"))
trip.data.long <- left_join(trip.data.long, county.file.ordered, by = c("end.adm2.code" = "ID_2"))

colnames(trip.data.long) <- c('start.adm2.code', 'end.adm2.code', 'date', 'trip.count', 'd', 'm', 'y',   
                              'start.adm1.name', 'start.adm2.name', 'X_start', 'Y_start', 'pop.start', 
                              'urb.2.5.start', 'urb.10.start', 'urb.19.start', 'start.adm1.code', 
                              'end.adm1.name', 'end.adm2.name', 'X_end', 'Y_end', 'pop.end', 
                              'urb.2.5.end', 'urb.10.end', 'urb.19.end', 'end.adm1.code')

# 5. Save files as .RData or .csv files for use in other scripts
save(trip.data.long, file = "tmp/NAM_trips_data_long.RData")
