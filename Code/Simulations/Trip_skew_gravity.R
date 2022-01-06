## Simulation to explore basic gravity model's ability to recreate asymmetry patterns
## Set-up remains the same (16 locations set up in a grid system with defined population size)
## Compare impact different skew scenarios

## Load libraries
library(asymmetry)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(png)
library(reshape2)
library(dplyr)
library(DescTools)

grav.model.skew <- function(N, D, theta, omega.1, omega.2, gamma, skew, title){
  
  trip.prop <- trips <- matrix(NA, nrow = length(N), ncol = length(N))
  for (i in 1:length(N)) {
    for (j in 1:length(N)) {
      trips[i,j] <- ifelse(
        i == j, 0,
        exp(log(theta) + (log(skew[i,j]) + omega.1*log(N[i]) + omega.2*log(N[j]) - gamma*log(D[i,j])))
      )
    }
    trip.prop[i,] <- (trips[i,]/sum(trips[i,])) ## trip proportion. 
  }
  colnames(trip.prop) <- seq(1,length(N),1)
  rownames(trip.prop) <- seq(1,length(N),1)
  plots <- asymmetry.plots(trip.prop, title)
  return(plots)
}

## Plotting asymmetry map and trip count/prop map
asymmetry.plots <- function(trips, title){
  trips.symmetry <- symmetry.metric(trips) # caluclate measures of symmetry
  trips.ss <- skewsymmetry(trips); summary(trips) # Calculate skew-symmetry SSQ 
  
  # creates a color palette from red to blue
  my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 299)
  col_breaks = c(seq(-4000,-.001,length=100),  # negative values are red
                 seq(-.001,0.001,length=100),  # zeroes are white
                 seq(0.001,4000,length=100))   # positive values are blue
  
  # heatmap for asymmetric portion 
  png("hmap.png", width = 400, height = 500)
  hmap(trips.ss, col = my_palette, xlab = "Destination", ylab = "Origin")
  dev.off()
  
  img1 <- readPNG("hmap.png")
  asym.map <- ggplot() + background_image(img1) + ggtitle(title)
  
  # tile plot of trips
  trips.df <- melt(trips, value.name = "trips", varnames = c('origin', 'destination'))
  
  trip.plot <- ggplot(trips.df, aes(x = destination, y = origin)) +
    geom_tile(aes(fill = trips))+
    xlab("Destination") +
    ylab("Origin ") +
    guides(fill = guide_legend(title = "Trips"))+
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11))+
    ggtitle(title)
  
  list(asym.map, trip.plot)
}

## Measuring symmetry
symmetry.metric <- function(trip.matrix) {
  trips <- trip.matrix
  trips.t <- t(trip.matrix)
  As = (abs(trips) + abs(trips.t))/2
  Aa = (abs(trips) - abs(trips.t))/2
  sym.metric = norm(Aa, type = "m")/norm(As, type = "m")
  return(sym.metric)
}

#### Generate trip details
# dummy.trips <- cbind.data.frame(origin = seq(1,16,1),  
#                                 x.coord = c(0, 25, 50, 75, 0, 25, 50, 75, 0, 25, 50, 75, 0, 25, 50, 75),
#                                 y.coord = c(0, 0, 0, 0, 25, 25, 25, 25, 50, 50, 50, 50, 75, 75, 75, 75),
#                                 pop = rep(10000, 16))
# dummy.trips$varied.pop <- runif(n = nrow(dummy.trips), 500, 500000)

# save(dummy.trips, file = "dummy.trip.df.RData")

load("dummy.trip.df.RData") # used a saved dataframe to be consistent between runs due to random number generator

D <- matrix(NA, nrow=dim(dummy.trips)[1], ncol=dim(dummy.trips)[1])
for (i in 1:nrow(dummy.trips)){
  for (j in 1:nrow(dummy.trips)){
    x.dist <- abs(dummy.trips$x.coord[i] - dummy.trips$x.coord[j])
    y.dist <- abs(dummy.trips$y.coord[i] - dummy.trips$y.coord[j])
    D[i,j] <- sqrt(x.dist^2 + y.dist^2)
  }
}

colnames(D) <- rownames(D) <- dummy.trips$origin


trip.details <- melt(D, value.name = "distance", varnames = c('origin', 'destination'))
trip.details <- left_join(trip.details, dummy.trips[ , c('origin','varied.pop')], by = c('origin'))
trip.details <- left_join(trip.details, dummy.trips[ , c('origin','varied.pop')], by = c('destination' = 'origin'))
colnames(trip.details) <- c('origin', 'destination', 'distance', 'pop.start', 'pop.end')

## Defining different skew scenarios
no_skew <- matrix(1, nrow=dim(dummy.trips)[1], ncol=dim(dummy.trips)[1]) 
skew_2mixed_low <- no_skew;  skew_2mixed_low[1,2] <- skew_2mixed_low[8,2] <- 0.5 ## all locations have no skew, except two routes with reduced skew
skew_2mixed_high <- no_skew;  skew_2mixed_high[1,2] <- skew_2mixed_high[8,2] <- 2 ## all locations have no skew, except two routes with increased skew
skew_all_mixed <- matrix(runif(dim(dummy.trips)[1]*dim(dummy.trips)[1], 1/5, 2), nrow=dim(dummy.trips)[1], ncol=dim(dummy.trips)[1]) # Random selection of skew between 0.5 - 2
skew_row <- no_skew; skew_row[5,] <- 2 # all trips made from a certain location have increased skew (2)
skew_col <- no_skew; skew_col[,5] <- 2 # all trips made to a certain location have increased skew (2)
skew_shortDist <- no_skew * ifelse(D < quantile(D)[2], 2,                              # Trips are skewed towards shorter trips; longer trips are penalized
                                     ifelse(D > quantile(D)[2] & D < quantile(D)[3], 1.5, 
                                            ifelse(D > quantile(D)[3] & D < quantile(D)[4], 0.75, 
                                                   0.5)))
skew_longDist <- no_skew * ifelse(D < quantile(D)[2], 0.5,                            # Trips are skewed towards longer trips; shorter trips are penalized
                                      ifelse(D > quantile(D)[2] & D < quantile(D)[3], 0.75, 
                                             ifelse(D > quantile(D)[3] & D < quantile(D)[4], 1.5, 
                                                    2)))
skew_smallPop<- no_skew * ifelse(dummy.trips$varied.pop < quantile(dummy.trips$varied.pop)[2], 2,   # Trips are skewed towards smaller populations; trips to larger populations are penalized
                                    ifelse(dummy.trips$varied.pop > quantile(dummy.trips$varied.pop)[2] & dummy.trips$varied.pop < quantile(dummy.trips$varied.pop)[3], 1.5, 
                                           ifelse(dummy.trips$varied.pop > quantile(dummy.trips$varied.pop)[3] & dummy.trips$varied.pop < quantile(dummy.trips$varied.pop)[4], 0.75, 
                                                  0.5))) 
skew_largePop <- no_skew * ifelse(dummy.trips$varied.pop < quantile(dummy.trips$varied.pop)[2], 0.5, # Trips are skewed towards larger populations; trips to smaller populations are penalized
                                     ifelse(dummy.trips$varied.pop > quantile(dummy.trips$varied.pop)[2] & dummy.trips$varied.pop < quantile(dummy.trips$varied.pop)[3], 0.75, 
                                            ifelse(dummy.trips$varied.pop > quantile(dummy.trips$varied.pop)[3] & dummy.trips$varied.pop < quantile(dummy.trips$varied.pop)[4], 1.5, 
                                                   2)))

## Common parameters
N = dummy.trips$varied.pop
theta = 1 
omega.1 = 1
omega.2 = 1
gamma = 1

plots_noSkew <- grav.model.skew(N,  D, theta, omega.1, omega.2, gamma, skew = no_skew, title = "    No Skew")
plots_skew_2mixed_low <- grav.model.skew(N,  D, theta, omega.1, omega.2, gamma, skew = skew_2mixed_low, title = "    Two routes skewed low")
plots_skew_2mixed_high <- grav.model.skew(N,  D, theta, omega.1, omega.2, gamma, skew = skew_2mixed_high, title = "    Two routes skewed high")
plots_skew_all_mixed <- grav.model.skew(N,  D, theta, omega.1, omega.2, gamma, skew = skew_all_mixed, title = "    All routes randomly skewed")
plots_skew_row <- grav.model.skew(N, D, theta, omega.1, omega.2, gamma, skew = skew_row, title = "    One origin (#5) skewed high for all")
plots_skew_col <- grav.model.skew(N, D, theta, omega.1, omega.2, gamma, skew = skew_col, title = "    One destination (#5) skewed high for all")
plots_skew_longDist <- grav.model.skew(N, D, theta, omega.1, omega.2, gamma, skew = skew_longDist, title = "    Skew favors longer distance trips")
plots_skew_shortDist <- grav.model.skew(N, D, theta, omega.1, omega.2, gamma, skew = skew_shortDist, title = "    Skew favors shorter distance trips")
plots_skew_largePop <- grav.model.skew(N, D, theta, omega.1, omega.2, gamma, skew = skew_largePop, title = "    Skew favors larger population sizes")
plots_skew_smallPop <- grav.model.skew(N, D, theta, omega.1, omega.2, gamma, skew = skew_smallPop, title = "    Skew favors smaller population sizes")


asym.heat.maps <- ggarrange(plots_noSkew[[1]], plots_skew_2mixed_low[[1]], plots_skew_2mixed_high[[1]], plots_skew_all_mixed[[1]], plots_skew_row[[1]],
                            plots_skew_col[[1]], plots_skew_longDist[[1]], plots_skew_shortDist[[1]], plots_skew_largePop[[1]], plots_skew_smallPop[[1]], 
                            ncol = 5, nrow = 2, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"))

prop.OD.maps <- ggarrange(plots_noSkew[[2]], plots_skew_2mixed_low[[2]], plots_skew_2mixed_high[[2]], plots_skew_all_mixed[[2]], plots_skew_row[[2]],
                          plots_skew_col[[2]], plots_skew_longDist[[2]], plots_skew_shortDist[[2]], plots_skew_largePop[[2]], plots_skew_smallPop[[2]], 
                          ncol = 5, nrow = 2, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"))