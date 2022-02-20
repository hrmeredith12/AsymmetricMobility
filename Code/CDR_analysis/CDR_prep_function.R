## Prepare data for asymmetry paper plots ###
## Input: Daily CDRs, adm details
## Output: Daily CDRs, formatted for passing to asymmetry functions
## Last updated: Feb 17, 2022
## Updated by: Hannah Meredith

library("geosphere")
library("dplyr")       # joining
library("lubridate")  # parsing dates
library("reshape")    # for melting data
library("tidyr")

## Namibia #####

load("../../Data/NAM/NAM_between_day_mobility_updated_joinedAdm2.RData")
load("../../Data/NAM/NAM.ID.link.RData")

NAM <- mydata.joined
NAM <- NAM %>%
  pivot_longer(cols = contains("/"),
               names_to = "trip.date",
               values_to = "trips")
dates.sep <- strsplit(as.character(NAM$trip.date), "/", perl=TRUE)          # separate date string
NAM$d <- as.integer(sapply(dates.sep, function(x) x[1]))                    # define day of date
NAM$m <- as.integer(sapply(dates.sep, function(x) x[2]))                    # month of date
NAM$y <- as.integer(sapply(dates.sep, function(x) x[3]))                    # year of date
NAM$trip.date <- as.Date(with(NAM, paste(y, m, d,sep="-")), "%Y-%m-%d")

NAM <- left_join(NAM, county.file.ordered[-c(5:8)], by = c("i" = "ID_2"))
NAM <- left_join(NAM, county.file.ordered[-c(5:8)], by = c("j" = "ID_2"))

colnames(NAM) <- c("start.adm2.code", "end.adm2.code", "trip.date", "trip.count", "d", "m", "y",
                   "start.adm1.name", "start.adm2.name", "X_start", "Y_start", "start.adm1.code",
                   "end.adm1.name", "end.adm2.name", "X_end", "Y_end", "end.adm1.code")
NAM$link.distance.km <- distHaversine(matrix(c(NAM$X_start, NAM$Y_start), ncol =2),   # distance in km
                                   matrix(c(NAM$X_end, NAM$Y_end), ncol = 2))/1000

NAM <- NAM[ , c("start.adm1.name", "start.adm2.name", "start.adm1.code", "start.adm2.code", "end.adm1.name", "end.adm2.name", "end.adm1.code", "end.adm2.code",
                "d", "m", "y", "trip.date", "trip.count", "X_start", "Y_start", "X_end", "Y_end", "link.distance.km")]

# save(NAM, file = "../../Data/NAM/NAM.daily.trips.RData")

NAM.daily.avg <- NAM %>%
  group_by(start.adm1.name, start.adm2.name, start.adm1.code, start.adm2.code, end.adm1.name, end.adm2.name, end.adm1.code, end.adm2.code) %>%
  summarise(trip.count.avg = mean(trip.count, na.rm = T)) %>%
  ungroup()

NAM.daily.avg <- NAM.daily.avg[order(NAM.daily.avg$start.adm2.code, NAM.daily.avg$end.adm2.code),]
# save(NAM.daily.avg, file = "../../Data/NAM/NAM.daily.trips.avg.RData")

NAM.daily.avg.matrix <-  NAM.daily.avg %>% select(start.adm2.code, end.adm2.code, trip.count.avg) %>%
  pivot_wider(names_from = end.adm2.code,
              values_from = trip.count.avg) 
NAM.daily.avg.matrix <- as.matrix(NAM.daily.avg.matrix)
rownames(NAM.daily.avg.matrix) <- NAM.daily.avg.matrix[,1]
NAM.daily.avg.matrix <- NAM.daily.avg.matrix[,-1]
NAM.daily.avg.matrix[is.na(NAM.daily.avg.matrix)] <- 0
NAM.daily.avg.matrix.order <- NAM.daily.avg.matrix[order(as.integer(rownames(NAM.daily.avg.matrix))), order(as.integer(colnames(NAM.daily.avg.matrix)))]
# save(NAM.daily.avg.matrix.order, file = "../../Data/NAM/NAM.daily.avg.matrix.RData")





## Kenya #####

#1. Import location details
fid <- read.csv("../../Data/KEN/KEN_adm2_coords_pop_urbanicity.csv")[c('fid','SP_ID', 'ADM2_NAME')]
map.file <- read.csv("../../Data/KEN/KEN_adm2_pop_urb_buildingFootprints_coords.csv")#read.csv('final/KEN_adm2_pop_urb_2_5_10_19_ppp_coord.csv')
map.file <- left_join(map.file, fid, by = c('SP_ID', 'ADM2_NAME'))

adm2.name<- read.csv('../../Data/KEN/KEN_District_pop_coords_fid.csv')[c('fid','ADM2_NAME')] #  use different file because spelling of # 63 and 69 are different #details <- read.csv('KEN_adm2_coords_pop_urbanicity.csv')[c(1,5)]  
adm2.name$ADM2_NAME <- toupper(adm2.name$ADM2_NAME)
adm2.name <- adm2.name[order(adm2.name$fid),]

adm.details <- map.file[c('fid', 'ADM2_NAME', 'ADM1_NAME', 'X_coord', 'Y_coord', 'pop2010sum', 'build_urb')]#1,9,18,19)]

adm1.codes <- read.csv('../../Data/KEN/KEN_adm1_coords.csv')[ , c('ID_1','NAME_1')]
# adm1.codes$NAME_1 <- as.character(adm1.codes$NAME_1)
adm.details <- left_join(adm.details, adm1.codes, by = c("ADM1_NAME" = "NAME_1"))
colnames(adm.details) <- c("adm2_ID", 'adm2_name', "adm1_name", "X_coord", "Y_coord", 'pop2010', 'build_urb', "adm1_ID")
adm.details <- adm.details [ , c("adm1_ID", "adm1_name", "adm2_ID", "adm2_name", "X_coord", "Y_coord", 'pop2010', 'build_urb')]
adm.details$adm2_name_CAPS <- toupper(adm.details$adm2_name)
adm.details$adm2_name_CAPS[adm.details$adm2_name_CAPS == "NANDI NORTH"] <- "NANDI"
adm.details$adm2_name_CAPS[adm.details$adm2_name_CAPS == "BUTERE MUMIAS"] <- "BUTERE/MUMIAS"


# import trip data and wrangle dates
trip.data <- read.csv('../../Data/KEN/KEN_entrances_per_day.csv', header = TRUE, check.names = FALSE)
trip.long <- melt(trip.data, id.vars = c("origin", "destination"))
colnames(trip.long) <- c("i", "j", "date", "trip.count")
trip.long <- left_join(trip.long, adm.details[, c("adm1_ID", "adm1_name", "adm2_ID", "adm2_name", "X_coord", "Y_coord",  
                                                  "adm2_name_CAPS")], by = c("i" = "adm2_name_CAPS")) 
trip.long <- left_join(trip.long, adm.details[, c("adm1_ID", "adm1_name", "adm2_ID", "adm2_name", "X_coord", "Y_coord",  
                                                  "adm2_name_CAPS")], by = c("j" = "adm2_name_CAPS"))#[-c(1,2)]
colnames(trip.long) <- c("start.adm2.name.Case","end.adm2.name.Case",
                         "trip.date", "trip.count",
                         "start.adm1.code", "start.adm1.name", "start.adm2.code", "start.adm2.name", "X_start", "Y_start", 
                         "end.adm1.code", "end.adm1.name", "end.adm2.code","end.adm2.name", "X_end", "Y_end")

trip.long <- subset(trip.long, start.adm2.name.Case != " ")
trip.long <- subset(trip.long, end.adm2.name.Case != " ")

# rewrite dates
dates.sep <- strsplit(as.character(trip.long$trip.date), "-", perl=TRUE)   # separate date string 
trip.long$d <- as.integer(sapply(dates.sep, function(x) x[3]))                   # define day of day
trip.long$m <- as.integer(sapply(dates.sep, function(x) x[2]))                    # month of date
trip.long$y <- as.integer(sapply(dates.sep, function(x) x[1]))                    # year of date

# import stay data and wrangle dates
stay.data <- read.csv('../../Data/KEN/KEN_stays_per_day.csv', header = TRUE, check.names = FALSE)
stay.data$destination <- stay.data$origin
stay.long <- melt(stay.data, id.vars = c("origin", "destination"))
colnames(stay.long) <- c("i", "j", "date", "trip.count")
stay.long <- left_join(stay.long, adm.details[, c("adm1_ID", "adm1_name", "adm2_ID", "adm2_name", "X_coord", "Y_coord", 
                                                  "adm2_name_CAPS")], by = c("i" = "adm2_name_CAPS")) 
stay.long <- left_join(stay.long, adm.details[, c("adm1_ID", "adm1_name", "adm2_ID", "adm2_name", "X_coord", "Y_coord", 
                                                  "adm2_name_CAPS")], by = c("j" = "adm2_name_CAPS"))#[-c(1,2)]
colnames(stay.long) <- c("start.adm2.name.Case","end.adm2.name.Case",
                         "trip.date", "trip.count",
                         "start.adm1.code", "start.adm1.name", "start.adm2.code", "start.adm2.name", "X_start", "Y_start", 
                         "end.adm1.code", "end.adm1.name", "end.adm2.code","end.adm2.name", "X_end", "Y_end")
stay.long <- subset(stay.long, start.adm2.name.Case != " ")

dates.sep <- strsplit(as.character(stay.long$trip.date), "/", perl=TRUE)   # separate date string 
stay.long$d <- as.integer(sapply(dates.sep, function(x) x[2]))                    # define day of day
stay.long$m <- as.integer(sapply(dates.sep, function(x) x[1]))                    # month of date
stay.long$y <- as.integer(sapply(dates.sep, function(x) x[3]))                    # year of date
stay.long$trip.date <- as.Date(with(stay.long, paste(y, m, d,sep="-")), "%Y-%m-%d")

trip.data.long <- rbind.data.frame(stay.long, trip.long)

trip.data.long <- trip.data.long[order(trip.data.long[,'y'],trip.data.long[,'m'], trip.data.long[,'d'],trip.data.long[ ,'start.adm2.code'], trip.data.long[ ,'end.adm2.code']), ]
trip.data.long$link.distance.km <- distHaversine(matrix(c(trip.data.long$X_start, trip.data.long$Y_start), ncol =2),   # distance in km
                                      matrix(c(trip.data.long$X_end, trip.data.long$Y_end), ncol = 2))/1000


KEN <- trip.data.long[ , c("start.adm1.name", "start.adm2.name", "start.adm1.code", "start.adm2.code", "end.adm1.name", "end.adm2.name", "end.adm1.code", "end.adm2.code",
                                     "d", "m", "y", "trip.date", "trip.count", "X_start", "Y_start", "X_end", "Y_end", "link.distance.km")]

KEN <- subset(KEN, !trip.date %in% as.Date(c("2008-06-01", "2009-03-04"))) ## Dates seemed to have aggregrated trip counts that are not in line with counts seen on other dates

# save(KEN, file = "../../Data/KEN/KEN.daily.trips.RData")

KEN.daily.avg <- KEN %>%
  group_by(start.adm1.name, start.adm2.name, start.adm1.code, start.adm2.code, end.adm1.name, end.adm2.name, end.adm1.code, end.adm2.code) %>%
  summarise(trip.count.avg = mean(trip.count, na.rm = T)) %>%
  ungroup()

KEN.daily.avg <- KEN.daily.avg[order(KEN.daily.avg$start.adm2.code, KEN.daily.avg$end.adm2.code),]

# save(KEN.daily.avg, file = "../../Data/KEN/KEN.daily.trips.avg.RData")

KEN.daily.avg.matrix <-  KEN.daily.avg %>% select(start.adm2.code, end.adm2.code, trip.count.avg) %>%
  pivot_wider(names_from = end.adm2.code,
              values_from = trip.count.avg) 
KEN.daily.avg.matrix <- as.matrix(KEN.daily.avg.matrix)
rownames(KEN.daily.avg.matrix) <- KEN.daily.avg.matrix[,1]
KEN.daily.avg.matrix <- KEN.daily.avg.matrix[,-1]
KEN.daily.avg.matrix[is.na(KEN.daily.avg.matrix)] <- 0
KEN.daily.avg.matrix.order <- KEN.daily.avg.matrix[order(as.integer(rownames(KEN.daily.avg.matrix))), order(as.integer(colnames(KEN.daily.avg.matrix)))]

# save(KEN.daily.avg.matrix.order, file = "../../Data/KEN/KEN.daily.avg.matrix.RData")
