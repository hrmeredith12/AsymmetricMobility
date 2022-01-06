## Simulation to explore basic gravity model's ability to recreate asymmetry patterns
## Set-up remains the same (16 locations set up in a grid system with defined population size)
## Compare impact different clustering scenarios. Does geographical layout matter?

## Load libraries
library(asymmetry)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(png)
library(reshape2)
library(dplyr)
library(DescTools)

grav.model.skew.cluster <- function(N, D, theta, omega.1, omega.2, gamma, skew, title){
  
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

grav.model.param.cluster <- function(N, D, theta, omega.1, omega.2, gamma, title){
  
  trip.prop <- trips <- matrix(NA, nrow = length(N), ncol = length(N))
  for (i in 1:length(N)) {
    for (j in 1:length(N)) {
      
      trips[i,j] <- ifelse(
        i == j,
        0,
        exp(log(theta) + (omega.1[i,j]*log(N[i]) + omega.2[i,j]*log(N[j]) - gamma[i,j]*log(D[i,j])))
      )
      
    }
    trip.prop[i,] <- (trips[i,]/sum(trips[i,]))
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
dummy.trips$row_cluster <- ifelse(dummy.trips$origin %in% c(1:4), 1,
                              ifelse(dummy.trips$origin %in% c(5:8), 2,
                                     ifelse(dummy.trips$origin %in% c(9:12), 3, 4)))
dummy.trips$col_cluster <- ifelse(dummy.trips$origin %in% c(1,5,9,13), 1,
                                  ifelse(dummy.trips$origin %in% c(2,6,10,14), 2,
                                         ifelse(dummy.trips$origin %in% c(3,7,11,15), 3, 4)))
dummy.trips$quad_cluster <- ifelse(dummy.trips$origin %in% c(1,2,5,6), 1,
                                  ifelse(dummy.trips$origin %in% c(3,4,7,8), 2,
                                         ifelse(dummy.trips$origin %in% c(9,10,13,14), 3, 4)))


D <- row_cluster <- col_cluster <- quad_cluster <- matrix(NA, nrow=dim(dummy.trips)[1], ncol=dim(dummy.trips)[1])
for (i in 1:nrow(dummy.trips)){
  for (j in 1:nrow(dummy.trips)){
    x.dist <- abs(dummy.trips$x.coord[i] - dummy.trips$x.coord[j])
    y.dist <- abs(dummy.trips$y.coord[i] - dummy.trips$y.coord[j])
    D[i,j] <- sqrt(x.dist^2 + y.dist^2)
    row_cluster[i,j] <- ifelse(dummy.trips$row_cluster[i] == dummy.trips$row_cluster[j], 1, 0) 
    col_cluster[i,j] <- ifelse(dummy.trips$col_cluster[i] == dummy.trips$col_cluster[j], 1, 0) 
    quad_cluster[i,j] <- ifelse(dummy.trips$quad_cluster[i] == dummy.trips$quad_cluster[j], 1, 0) 
  }
}

colnames(D) <- rownames(D) <- colnames(row_cluster) <- rownames(row_cluster) <- colnames(col_cluster) <- rownames(col_cluster) <- colnames(quad_cluster) <- rownames(quad_cluster) <- dummy.trips$origin


trip.details <- melt(D, value.name = "distance", varnames = c('origin', 'destination'))
trip.details <- left_join(trip.details, dummy.trips[ , c('origin','varied.pop')], by = c('origin'))
trip.details <- left_join(trip.details, dummy.trips[ , c('origin','varied.pop')], by = c('destination' = 'origin'))
colnames(trip.details) <- c('origin', 'destination', 'distance', 'pop.start', 'pop.end')

# 1. Account for clustering with Skew factor
ones <- matrix(1, nrow=dim(dummy.trips)[1], ncol=dim(dummy.trips)[1]) 
skew_row_cluster <- ones * ifelse(row_cluster == 1, 2, 0.5)
skew_row_cluster_capital <- skew_row_cluster; skew_row_cluster_capital[ ,dummy.trips$origin[which.max(dummy.trips$varied.pop)]] <- 2
skew_col_cluster <- ones * ifelse(col_cluster == 1, 2, 0.5)
skew_col_cluster_capital <- skew_col_cluster; skew_col_cluster_capital[ ,dummy.trips$origin[which.max(dummy.trips$varied.pop)]] <- 2
skew_quad_cluster <- ones * ifelse(quad_cluster == 1, 2, 0.5)
skew_quad_cluster_capital <- skew_quad_cluster; skew_quad_cluster_capital[ ,dummy.trips$origin[which.max(dummy.trips$varied.pop)]] <- 2

## Common parameters
N = dummy.trips$varied.pop
theta = 1 
omega.1 = 1
omega.2 = 1
gamma = 1

plots_row_skew_cluster <- grav.model.skew.cluster(N, D, theta, omega.1, omega.2, gamma, skew = skew_row_cluster, title = "    Trips are clustered by row")
plots_row_skew_cluster_capital <- grav.model.skew.cluster(N, D, theta, omega.1, omega.2, gamma, skew = skew_row_cluster_capital, title = "    Trips are clustered by row + capital")
plots_col_skew_cluster <- grav.model.skew.cluster(N, D, theta, omega.1, omega.2, gamma, skew = skew_col_cluster, title = "    Trips are clustered by column")
plots_col_skew_cluster_capital <- grav.model.skew.cluster(N, D, theta, omega.1, omega.2, gamma, skew = skew_col_cluster_capital, title = "    Trips are clustered by column + capital")
plots_quad_skew_cluster <- grav.model.skew.cluster(N, D, theta, omega.1, omega.2, gamma, skew = skew_quad_cluster, title = "    Trips are clustered by quadrants")
plots_quad_skew_cluster_capital <- grav.model.skew.cluster(N, D, theta, omega.1, omega.2, gamma, skew = skew_quad_cluster_capital, title = "    Trips are clustered by quad + capital")


asym.heat.maps_skew <- ggarrange(plots_row_skew_cluster[[1]], plots_row_skew_cluster_capital[[1]], plots_col_skew_cluster[[1]], plots_col_skew_cluster_capital[[1]], 
                            plots_quad_skew_cluster[[1]], plots_quad_skew_cluster_capital[[1]], 
                            ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E", "F"))



#2. Account for clustering with parameters

omega.1_row_cluster <- omega.2_row_cluster <- ones * ifelse(row_cluster == 1, 1.2, 0.8)
gamma_row_cluster <- ones * ifelse(row_cluster == 1, 2, 0.9)

omega.1_row_cluster_capital <- omega.2_row_cluster_capital <- omega.1_row_cluster; 
omega.1_row_cluster_capital[ ,dummy.trips$origin[which.max(dummy.trips$varied.pop)]] <- omega.2_row_cluster_capital[ ,dummy.trips$origin[which.max(dummy.trips$varied.pop)]] <- 1.2
gamma_row_cluster_capital <- gamma_row_cluster
gamma_row_cluster_capital[ ,dummy.trips$origin[which.max(dummy.trips$varied.pop)]] <- 2

omega.1_col_cluster <- omega.2_col_cluster <- ones * ifelse(col_cluster == 1, 1.2, 0.8)
gamma_col_cluster <- ones * ifelse(col_cluster == 1, 2, 0.9)

omega.1_col_cluster_capital <- omega.2_col_cluster_capital <- omega.1_col_cluster; 
omega.1_col_cluster_capital[ ,dummy.trips$origin[which.max(dummy.trips$varied.pop)]] <- omega.2_col_cluster_capital[ ,dummy.trips$origin[which.max(dummy.trips$varied.pop)]] <- 1.2
gamma_col_cluster_capital <- gamma_col_cluster
gamma_col_cluster_capital[ ,dummy.trips$origin[which.max(dummy.trips$varied.pop)]] <- 2

omega.1_quad_cluster <- omega.2_quad_cluster <- ones * ifelse(quad_cluster == 1, 1.2, 0.8)
gamma_quad_cluster <- ones * ifelse(quad_cluster == 1, 2, 0.9)

omega.1_quad_cluster_capital <- omega.2_quad_cluster_capital <- omega.1_quad_cluster; 
omega.1_quad_cluster_capital[ ,dummy.trips$origin[which.max(dummy.trips$varied.pop)]] <- omega.2_quad_cluster_capital[ ,dummy.trips$origin[which.max(dummy.trips$varied.pop)]] <- 1.2
gamma_quad_cluster_capital <- gamma_quad_cluster
gamma_quad_cluster_capital[ ,dummy.trips$origin[which.max(dummy.trips$varied.pop)]] <- 2

plots_row_cluster <- grav.model.param.cluster(N, D, theta, omega.1_row_cluster, omega.2_row_cluster, gamma_row_cluster, title = "    Trips are clustered by row")
plots_row_cluster_capital <- grav.model.param.cluster(N, D, theta, omega.1_row_cluster_capital, omega.2_row_cluster_capital, gamma_row_cluster_capital, title = "    Trips are clustered by row + capital")
plots_col_cluster <- grav.model.param.cluster(N, D, theta, omega.1_col_cluster, omega.2_col_cluster, gamma_col_cluster, title = "    Trips are clustered by column")
plots_col_cluster_capital <- grav.model.param.cluster(N, D, theta, omega.1_col_cluster_capital, omega.2_col_cluster_capital, gamma_col_cluster_capital, title = "    Trips are clustered by column + capital")
plots_quad_cluster <- grav.model.param.cluster(N, D, theta, omega.1_quad_cluster, omega.2_quad_cluster, gamma_quad_cluster, title = "    Trips are clustered by quadrants")
plots_quad_cluster_capital <- grav.model.param.cluster(N, D, theta, omega.1_quad_cluster_capital, omega.2_quad_cluster_capital, gamma_quad_cluster_capital, title = "    Trips are clustered by quad + capital")


asym.heat.maps <- ggarrange(plots_row_cluster[[1]], plots_row_cluster_capital[[1]], plots_col_cluster[[1]], plots_col_cluster_capital[[1]], 
                            plots_quad_cluster[[1]], plots_quad_cluster_capital[[1]], 
                            ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E", "F"))
prop.OD.maps <- ggarrange(plots_row_cluster[[2]], plots_row_cluster_capital[[2]], plots_col_cluster[[2]], plots_col_cluster_capital[[2]], 
                            plots_quad_cluster[[2]], plots_quad_cluster_capital[[2]], 
                            ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E", "F"))
