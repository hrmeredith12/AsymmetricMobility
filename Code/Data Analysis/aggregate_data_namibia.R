###############################################################################
# Ronan Corgel
# Asymmetry Project
# Aggregate Namibia Data
# February 2022
# Steps:
# 1. Load data from previous file
# 2. Calculate average proportion by month and week and overall
# 3. Create trip matrices and trip proportion matrices (Origin-Destination matrix)
# 4. Save files as .RData files for use in other scripts
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
library('asymmetry')  # asymmetry matrices

# 1. Load data from previous file
load('tmp/NAM_trips_data_long.RData')

# 2. Calculate average proportion overall and by month and week
# Overall
adm2.trip.overall <- trip.data.long %>%                                         # calculate number of trips between each origin/destination pair for each month/year
  group_by(start.adm2.code, end.adm2.code) %>%
  mutate(adm2.single.trip.sum = sum(trip.count, na.rm = TRUE)) %>%
  distinct(start.adm2.code, end.adm2.code, .keep_all=TRUE) 

# Calculate the total trips out of each origin each day (excluding stays)
adm2.trip.overall <- adm2.trip.overall %>% group_by(start.adm2.code, date) %>%
  mutate(trip.total = ifelse(start.adm2.code == end.adm2.code, 0, 
                             sum(adm2.single.trip.sum[which(start.adm2.code != end.adm2.code)], na.rm = T)))

# Calculate trip proportion of a specific trip between a given origin and destination
adm2.trip.overall$overall.trip.proportion = adm2.trip.overall$adm2.single.trip.sum/adm2.trip.overall$trip.total  

# Month
adm2.trip.month <- trip.data.long %>%                                           # calculate number of trips between each origin/destination pair for each month/year
  group_by(start.adm2.code, end.adm2.code, m, y) %>%
  mutate(adm2.single.trip.sum = sum(trip.count, na.rm = TRUE)) %>%
  distinct(start.adm2.code, end.adm2.code, m, y, .keep_all=TRUE) 

# Calculate the total trips out of each origin each day (excluding stays)
adm2.trip.month <- adm2.trip.month %>% group_by(start.adm2.code, date) %>%
  mutate(trip.total = ifelse(start.adm2.code == end.adm2.code, 0, 
                             sum(adm2.single.trip.sum[which(start.adm2.code != end.adm2.code)], na.rm = T)))

# Calculate trip proportion of a specific trip between a given origin and destination
adm2.trip.month$monthly.trip.proportion = adm2.trip.month$adm2.single.trip.sum/adm2.trip.month$trip.total  

# Calculate average proportion by month and average number of trips by month
adm2.trip.month.avg <- adm2.trip.month %>%                                      # calculate average trip proportion each month between origin/destination
  group_by(start.adm2.code, end.adm2.code) %>%
  mutate(avg.monthly.trip.proportion = mean(monthly.trip.proportion, na.rm = TRUE)) %>%
  mutate(var.monthly.trip.proportion = var(monthly.trip.proportion, na.rm = TRUE)) %>%
  mutate(adm2.single.trip.avg = ceiling(mean(adm2.single.trip.sum, na.rm = TRUE))) %>%
  mutate(adm2.single.trip.var = var(adm2.single.trip.sum, na.rm = TRUE)) %>%
  distinct(start.adm2.code, end.adm2.code, .keep_all = TRUE)

# Week
trip.data.long$year_week <- floor_date(trip.data.long$date, "1 week")
adm2.trip.week <- trip.data.long %>%                                            # calculate number of trips between each origin/destination pair for each month/year
  group_by(start.adm2.code, end.adm2.code, year_week) %>%
  mutate(adm2.single.trip.sum = sum(trip.count, na.rm = TRUE)) %>%
  distinct(start.adm2.code, end.adm2.code, year_week, .keep_all=TRUE) 

# Calculate the total trips out of each origin each day (excluding stays)
adm2.trip.week <- adm2.trip.week %>% group_by(start.adm2.code, date) %>%
  mutate(trip.total = ifelse(start.adm2.code == end.adm2.code, 0, 
                             sum(adm2.single.trip.sum[which(start.adm2.code != end.adm2.code)], na.rm = T)))

# Calculate trip proportion of a specific trip between a given origin and destination
adm2.trip.week$weekly.trip.proportion = adm2.trip.week$adm2.single.trip.sum/adm2.trip.week$trip.total  

adm2.trip.week.avg <- adm2.trip.week %>%                                        # calculate average trip proportion each month between origin/destination
  group_by(start.adm2.code, end.adm2.code) %>%
  mutate(avg.weekly.trip.proportion = mean(weekly.trip.proportion, na.rm = TRUE)) %>%
  mutate(var.weekly.trip.proportion = var(weekly.trip.proportion, na.rm = TRUE)) %>%
  mutate(adm2.single.trip.avg = ceiling(mean(adm2.single.trip.sum, na.rm = TRUE))) %>%
  mutate(adm2.single.trip.var = var(adm2.single.trip.sum, na.rm = TRUE)) %>%
  distinct(start.adm2.code, end.adm2.code, .keep_all = TRUE)

# Day Check (Single Day)
# Calculate the total trips out of each origin each day
adm2.trip.day <- trip.data.long %>% group_by(start.adm2.code, date) %>%
  mutate(trip.total = ifelse(start.adm2.code == end.adm2.code, 0, 
                                   sum(trip.count[which(start.adm2.code != end.adm2.code)], na.rm = T)))

# Calculate trip proportion of a specific trip between a given origin and destination
adm2.trip.day$daily.trip.proportion = adm2.trip.day$trip.count/adm2.trip.day$trip.total

# Filter to specific day
day.test <- adm2.trip.day %>% filter(date == "2010-10-02")

# Create Matrix
prop.day.summary <- day.test[,c('start.adm2.code', 'end.adm2.code', 'daily.trip.proportion')]
# Reshape to wide
M.day.prop <- reshape::cast(prop.day.summary, start.adm2.code ~ end.adm2.code)            
# Label rows with district numbers
rownames(M.day.prop) <- M.day.prop$start.adm2.code                           
# Get rid of the first column
M.day.prop <- M.day.prop[ ,-1]
class(M.day.prop) <- 'data.frame'
M.day.prop <- as.matrix(M.day.prop)
names(dimnames(M.day.prop)) <- c('origin', 'destination')
# Replace NAs with 0
M.day.prop[is.na(M.day.prop)] <- 0
M.day.prop[is.infinite(M.day.prop)] <- 0
# Skew Symmetric Matrix
SS.day.prop <- skewsymmetry(M.day.prop)

# creates a color palette from red to blue
my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 299)
col_breaks = c(seq(-50000,-.001,length=100),  # negative values are red
               seq(-.001,0.01,length=100),   # zeroes are white
               seq(0.01,50000,length=100))  # positive values are blue

# Make heat maps of asymmetric matrices
pdf(file='visuals/hmap_single_day.pdf')
hmap(SS.day.prop, col = my_palette, ylab = "origin", xlab = "destination", main = "Day")
dev.off()

# 3. Create trip matrices of trip proportion (Origin-Destination matrix)
# Overall
# Proportions Matrix
prop.overall.summary <- adm2.trip.overall[,c('start.adm2.code', 'end.adm2.code', 'overall.trip.proportion')]
# Reshape to wide
M.overall.prop <- reshape::cast(prop.overall.summary, start.adm2.code ~ end.adm2.code)            
# Label rows with district numbers
rownames(M.overall.prop) <- M.overall.prop$start.adm2.code                           
# Get rid of the first column
M.overall.prop <- M.overall.prop[ ,-1]
class(M.overall.prop) <- 'data.frame'
M.overall.prop <- as.matrix(M.overall.prop)
names(dimnames(M.overall.prop)) <- c('origin', 'destination')
# Replace NAs with 0
M.overall.prop[is.na(M.overall.prop)] <- 0
M.overall.prop[is.infinite(M.overall.prop)] <- 0

# Trips Matrix
trip.overall.summary <- adm2.trip.overall[,c('start.adm2.code', 'end.adm2.code', 'adm2.single.trip.sum')]
# Reshape to wide
M.overall.trips <- reshape::cast(trip.overall.summary, start.adm2.code ~ end.adm2.code)            
# Label rows with district numbers
rownames(M.overall.trips) <- M.overall.trips$start.adm2.code                           
# Get rid of the first column
M.overall.trips <- M.overall.trips[ ,-1]
class(M.overall.trips) <- 'data.frame'
M.overall.trips <- as.matrix(M.overall.trips)
names(dimnames(M.overall.trips)) <- c('origin', 'destination')
# Replace NAs with 0
M.overall.trips[is.na(M.overall.trips)] <- 0

# Month
# Proportions Matrix
prop.month.avg.summary <- adm2.trip.month.avg[,c('start.adm2.code', 'end.adm2.code', 'avg.monthly.trip.proportion')]
# Reshape to wide
M.monthly.avg.prop <- reshape::cast(prop.month.avg.summary, start.adm2.code ~ end.adm2.code)            
# Label rows with district numbers
rownames(M.monthly.avg.prop) <- M.monthly.avg.prop$start.adm2.code                           
# Get rid of the first column
M.monthly.avg.prop <- M.monthly.avg.prop[ ,-1]
class(M.monthly.avg.prop) <- 'data.frame'
M.monthly.avg.prop <- as.matrix(M.monthly.avg.prop)
names(dimnames(M.monthly.avg.prop)) <- c('origin', 'destination')
# Replace NAs with 0
M.monthly.avg.prop[is.na(M.monthly.avg.prop)] <- 0
M.monthly.avg.prop[is.infinite(M.monthly.avg.prop)] <- 0

# Trips Matrix
trip.month.avg.summary <- adm2.trip.month.avg[,c('start.adm2.code', 'end.adm2.code', 'adm2.single.trip.avg')]
# Reshape to wide
M.monthly.avg <- reshape::cast(trip.month.avg.summary, start.adm2.code ~ end.adm2.code)            
# Label rows with district numbers
rownames(M.monthly.avg) <- M.monthly.avg$start.adm2.code                           
# Get rid of the first column
M.monthly.avg <- M.monthly.avg[ ,-1]
class(M.monthly.avg) <- 'data.frame'
M.monthly.avg <- as.matrix(M.monthly.avg)
names(dimnames(M.monthly.avg)) <- c('origin', 'destination')
# Replace NAs with 0
M.monthly.avg[is.na(M.monthly.avg)] <- 0

# Week
# Proportions Matrix
prop.week.avg.summary <- adm2.trip.week.avg[,c('start.adm2.code', 'end.adm2.code', 'avg.weekly.trip.proportion')]
# Reshape to wide
M.weekly.avg.prop <- reshape::cast(prop.week.avg.summary, start.adm2.code ~ end.adm2.code)            
# Label rows with district numbers
rownames(M.weekly.avg.prop) <- M.weekly.avg.prop$start.adm2.code                           
# Get rid of the first column
M.weekly.avg.prop <- M.weekly.avg.prop[ ,-1]
class(M.weekly.avg.prop) <- 'data.frame'
M.weekly.avg.prop <- as.matrix(M.weekly.avg.prop)
names(dimnames(M.weekly.avg.prop)) <- c('origin', 'destination')
# Replace NAs with 0
M.weekly.avg.prop[is.na(M.weekly.avg.prop)] <- 0
M.weekly.avg.prop[is.infinite(M.weekly.avg.prop)] <- 0

# Trips Matrix
trip.week.avg.summary <- adm2.trip.week.avg[,c('start.adm2.code', 'end.adm2.code', 'adm2.single.trip.avg')]
# Reshape to wide
M.weekly.avg <- reshape::cast(trip.week.avg.summary, start.adm2.code ~ end.adm2.code)            
# Label rows with district numbers
rownames(M.weekly.avg) <- M.weekly.avg$start.adm2.code                           
# Get rid of the first column
M.weekly.avg <- M.weekly.avg[ ,-1]
class(M.weekly.avg) <- 'data.frame'
M.weekly.avg <- as.matrix(M.weekly.avg)
names(dimnames(M.weekly.avg)) <- c('origin', 'destination')
# Replace NAs with 0
M.weekly.avg[is.na(M.weekly.avg)] <- 0

# 4. Save files as .RData files for use in other scripts
save(M.overall.prop, file = 'output/NAM_overall_prop_ODMatrix.RData')
save(M.monthly.avg.prop, file = 'output/NAM_monthly_prop_ODMatrix.RData')
save(M.weekly.avg.prop, file = 'output/NAM_weekly_prop_ODMatrix.RData')
save(M.overall.trips, file = 'output/NAM_overall_trips_ODMatrix.RData')
save(M.monthly.avg, file = 'output/NAM_monthly_trips_ODMatrix.RData')
save(M.weekly.avg, file = 'output/NAM_weekly_trips_ODMatrix.RData')
# Save monthly data for next file
save(adm2.trip.month, file = 'tmp/adm2.trip.month.RData')
