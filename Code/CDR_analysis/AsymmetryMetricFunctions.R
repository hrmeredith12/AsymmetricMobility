## Calculating different metrics for network symmetry

library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(scales)

#### Entropy Index ####

Entropy.Index <- function(df, location.name, trip.source){
  EI <- df %>%
    group_by(trip.date)%>%
    mutate(links = n(),
           total.trips = sum(trip.count, na.rm = T),
           Z.l = trip.count/total.trips,
           Z.l.term = ifelse(Z.l == 0, 0, Z.l*log(Z.l)/log(links)),
           entropy.idx = -sum(Z.l.term),
           entropy.date = as.Date(paste("2010", m, d,sep="-"), "%Y-%m-%d")) %>%
    distinct(entropy.date, y, m, entropy.idx)
  
  EI$location.name <- rep(location.name, nrow = nrow(EI))
  EI$trip.source <- rep(trip.source, nrow = nrow(EI))
  EI$EI_mean <- mean(EI$entropy.idx, na.rm = T)
  EI$EI_sd <- sd(EI$entropy.idx, na.rm = T)
  return(EI)
}

#### Node Symmetry Index (NSI) #### 

Node.symmetry.index <- function(df, location.name, trip.source){
  NSI <- df %>%
    group_by(trip.date, start.adm2.code)%>%
    mutate(outgoing.trips = sum(trip.count[start.adm2.code != end.adm2.code], na.rm = T))

  NSI <- NSI %>%
    group_by(trip.date, end.adm2.code)%>%
    mutate(incoming.trips = sum(trip.count[start.adm2.code != end.adm2.code], na.rm = T)) %>% ungroup()
  
  NSI$NSI <- (NSI$incoming.trips - NSI$outgoing.trips)/(NSI$incoming.trips + NSI$outgoing.trips)
  NSI <- subset(NSI, start.adm2.code == end.adm2.code) # since each origin has unique NSI, regardless of destinations, just select one row/location
  
  NSI.overall <- NSI %>%
    group_by(start.adm2.code)%>%
    summarise(NSI.yearly.avg = mean(NSI, na.rm = T),
              NSI.yearly.sd = sd(NSI, na.rm = T))%>%
    distinct(start.adm2.code, .keep_all = T)
  
  NSI.overall$location.name <- rep(location.name, nrow = nrow(NSI.overall))
  NSI.overall$trip.source <- rep(trip.source, nrow = nrow(NSI.overall))
  
  return(NSI.overall)
}

## Mapping NSI mean and variation
library(ggplot2)
library(maptools)
library(rgeos)
library(Cairo)
library(ggmap)
library(scales)
library(RColorBrewer)
library(sf)
library(sp)
library(rgdal)
library(readxl)
library(dplyr)

NAM.NSI.map <- function(NSI_NAM){
  set.seed(8000)
  states.shp <- readOGR("../../Data/NAM/NAM_adm2_joined.shp")
  qgis.details <- read.csv("../../Data/NAM/NAM_adm2Joined_pop_urb_buildingFootprints_coords.csv")#("../prep code for gm data/NAM/NAM_adm2Joined_pop_urb10-19_pph_ppp_coords.csv")
  qgis.details$urb.cat.2 <- ifelse(qgis.details$build_urb < 0.5, 0, 1)
  qgis.details <- qgis.details[ , colnames(qgis.details) %in% c("ID_1", "NAME_1", "ID_2", "NAME_2", "build_urb", "urb.cat.2", "X_coord", "Y_coord")]
  
  load("../../Data/NAM/NAM_adm2Joined_monthly_trips_details.RData")
  adm2.data.cdr <- distinct(adm2.trip.month.summary[, c(1:6)])
  
  key.adm2 <- left_join(qgis.details, adm2.data.cdr[ , c("start.adm1.name", "start.adm1.code", "start.adm2.name", "start.adm2.code")], by = c("NAME_1" = "start.adm1.name", "NAME_2" = "start.adm2.name"))
  key.adm2$start.adm1.code[key.adm2$NAME_1 == "Oshana" & is.na(key.adm2$start.adm1.code)] = 11
  key.adm2$start.adm2.code[key.adm2$NAME_2 == "Okatana"] = 82
  key.adm2$start.adm2.code[key.adm2$NAME_2 == "Ompundja"] = 84
  key.adm2 <- left_join(key.adm2, NSI_NAM, by = c("start.adm2.code"))
  
  states.shp.2 <- merge(states.shp, key.adm2, by.x = c("NAME_1", "NAME_2", "ID_1", "ID_2", "build_urb"), by.y = c("NAME_1", "NAME_2", "ID_1", "ID_2", "build_urb"))
  states.shp.f <- fortify(states.shp.2, region = "start.adm2.code", name = "NAME_2")
  class(states.shp.f)
  
  ## merge with coefficients and reorder
  merge.shp.coef <- merge(states.shp.f, key.adm2, by.x= "id",  by.y = "start.adm2.code", all.x = TRUE)
  final.plot <- merge.shp.coef[order(merge.shp.coef$order), ]
  
  ## plot
  NSI.map <- ggplot() +
    geom_polygon(data = final.plot,
                 aes(x = long, y = lat, group = group, fill = NSI.yearly.avg),  # could have fill depend on variable (i.e. population density)
                 color = "black", size = 0.25)+
    coord_map() + 
    scale_fill_distiller(name = "Node Strength Index", 
                         palette = "YlOrRd", 
                         breaks = c(-0.015, -0.01, -0.005, 0, 0.005, 0.01, 0.015),#pretty_breaks(n = 3),
                         labels = c('-0.015', '-0.01', '-0.005', '0', '0.005', '0.01', '0.015'),
                         limits = c(min = -0.015, max = 0.015)) +
    # labs(title = "Average Daily NSI")+
    theme(legend.position = c(0.85, 0.45),
          axis.ticks = element_blank(), 
          axis.title = element_blank(), 
          axis.text =  element_blank(), 
          panel.background = element_blank())
  
  NSI.var.map <- ggplot() +
    geom_polygon(data = final.plot,
                 aes(x = long, y = lat, group = group, fill = NSI.yearly.sd),  # could have fill depend on variable (i.e. population density)
                 color = "black", size = 0.25)+
    coord_map() + 
    scale_fill_distiller(name = "Variance in NSI", 
                         palette = "YlGnBu", 
                         breaks = c(0, 0.05, 0.1, 0.2, 0.3, 0.4),
                         labels = c('0', '0.05', '0.1', '0.2', '0.3', '0.4'),
                         limits = c(min = 0, max = 0.4)) +
    labs(title = "Variance in Daily NSI")+
    theme(legend.position = c(0.85, 0.45),
          axis.ticks = element_blank(), 
          axis.title = element_blank(), 
          axis.text =  element_blank(), 
          panel.background = element_blank())
  return(NSI.map)
}
KEN.NSI.map <- function(NSI_KEN){
  set.seed(8000)
  
  states.shp <- readOGR("../../Data/KEN/Kenya_adm2_09.shp") # import shape file
  region.shp <- readOGR("../../Data/KEN/KEN_adm1.shp")
  
  #import city locations and names
  adm2.data <- read.csv("../../Data/KEN/KEN_adm2_pop_urb_buildingFootprints_coords.csv") 
  names(adm2.data)[1] <- "fid"
  adm1.data <- read.csv('../../Data/KEN/KEN_adm1.csv')[, 5:6]
  
  region.data <- as.data.frame(coordinates(region.shp))
  colnames(region.data) <- c('Longitude', 'Latitude')
  region.data$NAME <- region.shp$NAME_1
  region.data <- left_join(region.data, adm1.data, by = c("NAME" = "NAME_1"))
  
  ## input data to plot on map
  adm1.data$NAME_1 <- factor(adm1.data$NAME_1, 
                             levels = c("Central", "Coast", "Eastern", "Nairobi","North-Eastern","Nyanza","Rift Valley", "Western"),
                             labels = c("Central", "Coast", "Eastern", "Nairobi","North Eastern","Nyanza","Rift Valley", "Western"))
  adm2.data <- left_join(adm2.data, adm1.data, by = c("ADM1_NAME" = "NAME_1"))
  adm2 <- adm2.data[, c("X_coord", "Y_coord", "ID_1", "fid", "ADM1_NAME", "ADM2_NAME")]
  
  num.states <- length(states.shp$ADM1_NAME)
  mydata<-data.frame(Name_1 = states.shp$ADM1_NAME, id = states.shp$ADM1_CODE, Name_2 = states.shp$ADM2_NAME, id.2 = states.shp$ADM2_CODE, popDensity = rnorm(num.states, 55,20))
  
  # need to reorder files because ID_2/NAME_2 of the NAM_adm2_centroids.. don't match NAM_adm2.shp
  mydata.2 <- merge(mydata, adm2, by.x = "Name_2", by.y = "ADM2_NAME")
  mydata.2 <- left_join(mydata.2, NSI_KEN, by = c("Name_2" = "start.adm2.name"))
  
  # fortify shape file to get into dataframe
  states.shp.f <- fortify(states.shp, region = "ADM2_CODE", name = "ADM2_NAME")
  class(states.shp.f)
  
  ## merge with coefficients and reorder
  merge.shp.coef <- merge(states.shp.f, mydata.2, by.x= "id",  by.y = "id.2", all.x = TRUE)
  final.plot <- merge.shp.coef[order(merge.shp.coef$order), ]
  
  ## plot
  ggplot() +
    geom_polygon(data = final.plot,
                 aes(x = long, y = lat, group = group, fill = NSI.yearly.avg),  
                 color = "black", size = 0.25)+
    coord_map() + 
    scale_fill_distiller(name = "Node Strength Index", 
                         palette = "YlOrRd", 
                         breaks = c(-0.05, -0.01, 0, 0.01, 0.05),
                         labels = c('-0.05','-0.01', '0', '0.01', '0.05'),
                         limits = c(min = -0.17, max = 0.03)) +
    # labs(title = "Kenya urbanicity")+
    theme(legend.position = c(1, 0.45),
          axis.ticks = element_blank(), 
          axis.title = element_blank(), 
          axis.text =  element_blank(), 
          panel.background = element_blank())
}



