## Calculating different metrics for network symmetry
library(asymmetry)
library(Cairo)
library(dplyr)
library(ggmap)
library(ggplot2)
library(lubridate)
library(maptools)
library(RColorBrewer)
library(readxl)
library(rgdal)
library(rgeos)
library(scales)
library(sf)
library(sp)
library(tidyr)

#### Entropy Index ####
## Measures extent to which the trips are distributed evenly across all links in the network (0 = one node dominates, 1 = trips are evenly distributed)

Entropy.Index <- function(df, location.name, trip.source){
  EI <- df %>%
    group_by(trip.date)%>%
    mutate(links = n(),   # number of OD pairs
           total.trips = sum(trip.count, na.rm = T),   
           Z.l = trip.count/total.trips,  # proportion of journeys on link l in relation to the total number of journeys in the network
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
  NSI$total.trips = NSI$incoming.trips + NSI$outgoing.trips
  
  # ggplot(NSI, aes(outgoing.trips, incoming.trips)) +
  #   geom_point()+
  #   geom_abline(intercept = 0, slope = 1) +
  #   scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
  #   scale_x_continuous(labels = function(x) format(x, scientific = TRUE))+
  #   facet_wrap(~start.adm2.name, nrow = 9, scales = "free")+
  #   theme(axis.text = element_text(size = 6))
  
  
  NSI.overall <- NSI %>%
    group_by(start.adm2.code, start.adm2.name, pop.start, urb.start)%>%
    summarise(total.trips.avg = mean(total.trips, na.rm = T),
              total.trips.sd = sd(total.trips, na.rm = T),
              incoming.trips.avg = mean(incoming.trips, na.rm = T),
              incoming.trips.sd = sd(incoming.trips, na.rm = T),
              outgoing.trips.avg = mean(outgoing.trips, na.rm = T),
              outgoing.trips.sd = sd(outgoing.trips, na.rm = T),
              NSI.yearly.avg = mean(NSI, na.rm = T),
              NSI.yearly.sd = sd(NSI, na.rm = T))%>%
    distinct(start.adm2.code, .keep_all = T)
  
  # ggplot(NSI.overall, aes(outgoing.trips.avg, incoming.trips.avg)) +
  #   geom_point()+
  #   geom_abline(intercept = 0, slope = 1) +
  #   scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
  #   scale_x_continuous(labels = function(x) format(x, scientific = TRUE))
  
  NSI.overall$location.name <- rep(location.name, nrow = nrow(NSI.overall))
  NSI.overall$trip.source <- rep(trip.source, nrow = nrow(NSI.overall))
  
  return(NSI.overall)
}

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
  
  states.shp.2 <- merge(states.shp, key.adm2, by.x = c("NAME_1", "NAME_2", "ID_1", "ID_2","pop2010ppp", "build_urb"), by.y = c("NAME_1", "NAME_2", "ID_1", "ID_2", "pop.start", "urb.start"))
  # Key_code_name_NAM <- states.shp.2@data[, c("NAME_1", "NAME_2", "start.adm1.code", "start.adm2.code","pop2010ppp", "build_urb")]  ## used to link code with names in Basic gravity model function
  # colnames(Key_code_name_NAM) <- c("start.adm1.name", "start.adm2.name", "start.adm1.code", "start.adm2.code", "start.pop", "urb.pop")
  # save(Key_code_name_NAM, file = "../../Data/NAM/Key_code_name_NAM.RData")
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
    scale_fill_distiller(name = "Average NSI", 
                         palette = "YlOrRd", 
                         breaks = c(-0.015, -0.01, -0.005, 0, 0.005, 0.01, 0.015),#pretty_breaks(n = 3),
                         labels = c('-0.015', '-0.01', '-0.005', '0', '0.005', '0.01', '0.015'),
                         limits = c(min = -0.015, max = 0.015)) +
    # guides(fill = guide_colourbar(barwidth = 1, barheight = 15))+
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
    theme(legend.position = c(0.85, 0.45),
          axis.ticks = element_blank(), 
          axis.title = element_blank(), 
          axis.text =  element_blank(), 
          panel.background = element_blank())
  maps <- list(NSI.map, NSI.var.map)
  return(maps)
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
  adm2 <- adm2.data[, c("X_coord", "Y_coord", "ID_1", "fid", "ADM1_NAME", "ADM2_NAME", "build_urb", "pop2010sum")]
  adm2 <- adm2[order(adm2$ID_1),]
  adm2$ID_2 <- seq(1:nrow(adm2))
  
  num.states <- length(states.shp$ADM1_NAME)
  mydata<-data.frame(Name_1 = states.shp$ADM1_NAME, id = states.shp$ADM1_CODE, Name_2 = states.shp$ADM2_NAME, id.2 = states.shp$ADM2_CODE, popDensity = rnorm(num.states, 55,20))
 
  mydata.2 <- merge(mydata, adm2, by.x = "Name_2", by.y = "ADM2_NAME")
  
  # Key_code_name_KEN <- mydata.2[, c("Name_1", "Name_2", "ID_1", "ID_2", "pop2010sum", "build_urb")]  ## used to link code with names in Basic gravity model function
  # colnames(Key_code_name_KEN) <- c("start.adm1.name", "start.adm2.name", "start.adm1.code", "start.adm2.code", "pop.start", "urb.start")
  # save(Key_code_name_KEN, file = "../../Data/KEN/Key_code_name_KEN.RData")

  mydata.2 <- left_join(mydata.2, NSI_KEN, by = c("Name_2" = "start.adm2.name"))
  
  
  # fortify shape file to get into dataframe
  states.shp.f <- fortify(states.shp, region = "ADM2_CODE", name = "ADM2_NAME")
  class(states.shp.f)
  
  ## merge with coefficients and reorder
  merge.shp.coef <- merge(states.shp.f, mydata.2, by.x= "id",  by.y = "id.2", all.x = TRUE)
  final.plot <- merge.shp.coef[order(merge.shp.coef$order), ]
  
  
  NSI.map <- ggplot() +
    geom_polygon(data = final.plot,
                 aes(x = long, y = lat, group = group, fill = NSI.yearly.avg),  # could have fill depend on variable (i.e. population density)
                 color = "black", size = 0.25)+
    coord_map() + 
    scale_fill_distiller(name = "Average NSI", 
                         palette = "YlOrRd", 
                         breaks = c(-0.015, -0.01, -0.005, 0, 0.005, 0.01, 0.015),#pretty_breaks(n = 3),
                         labels = c('-0.015', '-0.01', '-0.005', '0', '0.005', '0.01', '0.015'),
                         limits = c(min = -0.021, max = 0.021)) +
    # guides(fill = guide_colourbar(barwidth = 1, barheight = 15))+
    theme(#legend.position = c(0.85, 0.45),
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
    theme(#legend.position = c(0.85, 0.45),
          axis.ticks = element_blank(), 
          axis.title = element_blank(), 
          axis.text =  element_blank(), 
          panel.background = element_blank())
  maps <- list(NSI.map, NSI.var.map)
  return(maps)
}

create.NSI.figure <- function(NSI.obs, NSI.all, location){
  
  income.outgoing.trips <- ggplot(NSI.obs, aes(incoming.trips.avg, outgoing.trips.avg))+
    geom_point()+
    geom_errorbar(aes(ymin = incoming.trips.avg - incoming.trips.sd, ymax = incoming.trips.avg + incoming.trips.sd))+
    geom_errorbar(aes(xmin = outgoing.trips.avg - outgoing.trips.sd, xmax = outgoing.trips.avg + outgoing.trips.sd))+
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap(~location.name, nrow = 2, scales = "free")+
    scale_y_log10()+
    scale_x_log10()+
    labs(x = "Average daily incoming trips", y = "Average daily outgoing trips")+
    theme(panel.background = element_blank(),
          strip.background =element_rect(fill="white"),
          strip.text = element_text(colour = 'white'),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5))
  
  NSI.pop <- ggplot(NSI.obs, aes(pop.start, NSI.yearly.avg))+
    geom_point()+
    geom_errorbar(aes(ymin = NSI.yearly.avg - NSI.yearly.sd, ymax = NSI.yearly.avg + NSI.yearly.sd))+
    scale_x_log10()+
    facet_wrap(~location.name, nrow = 2, scales = "free")+
    labs(x = "Population size", y = "Average daily NSI")+
    theme(legend.position = c(0.8, 0.9),
          panel.background = element_blank(),
          strip.background =element_rect(fill="white"),
          strip.text = element_text(colour = 'white'),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5))
  
  NSI.urb <- ggplot(NSI.obs, aes(urb.start, NSI.yearly.avg))+
    geom_point()+
    geom_errorbar(aes(ymin = NSI.yearly.avg - NSI.yearly.sd, ymax = NSI.yearly.avg + NSI.yearly.sd))+
    facet_wrap(~location.name, nrow = 2, scales = "free")+
    labs(x = "Urbanicity", y = "Average daily NSI")+
    theme(legend.position = c(0.8, 0.9),
          panel.background = element_blank(),
          strip.background =element_rect(fill="white"),
          strip.text = element_text(colour = 'white'),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5))
  
  if(location == "Namibia"){
    NSI_map = NAM.NSI.map(NSI.obs)
  } else{
    NSI_map = KEN.NSI.map(NSI.obs)
  }
  
 NSI.compare.model <- ggplot(NSI.all, aes(NSI.yearly.avg, fill = trip.source))+
    geom_density(alpha = 0.5)+
    labs(x = "Averge Daily NSI", y = "Density", fill = "Trip count source")+
    lims(x = c(-0.05,0.05))+
    scale_y_sqrt()+
    theme(legend.position = c(0.8, 0.8),
          panel.background = element_blank(),
          strip.background =element_rect(fill="white"),
          strip.text = element_text(colour = 'white'))
  
  NSI_plots <- ggarrange(income.outgoing.trips, NSI.pop, NSI.urb, NSI_map[[1]], NSI_map[[2]], NSI.compare.model,
                            ncol = 3, nrow = 2,
                         labels = c("A", "B", "C", "D", "E", "F"))
  
  return(NSI_plots)
}


#### Transitive Symmetry Index ####

transitive.symmetry <- function(df, location.name, trip.source){
  trip.dates <- df %>% select(trip.date, m, y) %>% distinct(trip.date, m, y)
  asym.df <- rep()
  
  ## Calculate skew symmetry for each date. Make sure matrices are ordered correctly, otherwise they will not be ordered correctly in heatmap/transitive matrix
  for(i in 1:nrow(trip.dates)){  
    print(trip.dates$trip.date[i])
    
    # turn trips from date[i] into matrix
    trips.m <- subset(df, trip.date == trip.dates$trip.date[i])
    trips.m <- trips.m %>% select(start.adm2.code, end.adm2.code,trip.count) %>%
      pivot_wider(names_from = end.adm2.code, values_from = trip.count) 
    trips.m <- as.matrix(trips.m)
    rownames(trips.m) <- trips.m[,1]
    trips.m <- trips.m[,-1]
    trips.m[is.na(trips.m)] <- 0
    trips.m.order <- trips.m[order(as.integer(rownames(trips.m))), order(as.integer(colnames(trips.m)))]
    
    # calculate skew-symmetry matrix for date[i]
    trips.ss <- skewsymmetry(trips.m.order)
    asym.matrx <- trips.ss$A
    asym.df.i <- reshape2::melt(asym.matrx, value.name = "skew.symmetry", varnames=c('start.adm2.code', 'end.adm2.code'))
    asym.df.i$trip.date <- rep(trip.dates$trip.date[i], nrow = nrow(asym.df.i))
    
    # collect all skew_symmetry matrices. Used to calculate average skew-symmetries next 
    asym.df <- rbind(asym.df, asym.df.i)
  }
  
  # Average skew-symmetry across all dates and turn into matrix
  asym.overall <- asym.df %>% group_by(start.adm2.code, end.adm2.code)%>%
    summarise(skew.sym.yearly.avg = mean(skew.symmetry, na.rm = T),
              skew.sym.yearly.sd = sd(skew.symmetry, na.rm = T))%>%
    distinct(start.adm2.code, end.adm2.code, .keep_all = T)
  asym.wide <- asym.overall %>% select(start.adm2.code, end.adm2.code, skew.sym.yearly.avg) %>%
    pivot_wider(names_from = end.adm2.code, values_from = skew.sym.yearly.avg) 
  asym.wide <- as.matrix(asym.wide)
  rownames(asym.wide) <- asym.wide[,1]
  asym.wide <- asym.wide[,-1]
  asym.wide.order <- asym.wide[order(as.integer(rownames(asym.wide))), order(as.integer(colnames(asym.wide)))]
  
  # Order skew-symmetry matrix by number of negative or positive cells per row (origin)
  bin <- asym.wide.order * 0
  bin[asym.wide.order > 0] <- 1
  bin[asym.wide.order < 0] <- -1
  rsbin <- rowSums(bin)
  asym.matrx.order <- asym.wide.order[order(rsbin), order(rsbin)]
  
  # Only interested in top half of the matrix to characterize transitive nature of asymmetry
  asym.matrx.order[lower.tri(asym.matrx.order)] <- NA
  asym.long <- reshape2::melt(asym.matrx.order, value.name = "skew.symmetry", varnames=c('start.adm2.code', 'end.adm2.code'))
  asym.plot <- subset(asym.long, !is.na(skew.symmetry) & start.adm2.code != end.adm2.code)
  asym.plot$pos.neg <- ifelse(asym.plot$skew.symmetry > 0, 1, 
                              ifelse(asym.plot$skew.symmetry < 0, -1, 0)) 
  asym.plot$sqrt.skew <- sqrt(abs(asym.plot$skew.symmetry))*asym.plot$pos.neg
  asym.plot$location.name <- rep(location.name, nrow = nrow(asym.plot))
  asym.plot$trip.source <- rep(trip.source, nrow = nrow(asym.plot))
  return(asym.plot)
}

