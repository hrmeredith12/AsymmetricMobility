## Simulation to explore gravity model's ability to recreate asymmetry patterns
library(asymmetry)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(ggplot2)

## Basic gravity model
grav.model.basic <- function(N_orig, N_dest, D, theta, alpha, beta, gamma, title){
  
  trip.prop <- trips <- matrix(NA, nrow = length(N_orig), ncol = length(N_dest))
  for (i in 1:length(N_orig)) {
    for (j in 1:length(N_dest)) {
      
      trips[i,j] <- ifelse(
        i == j, 0,
        exp(log(theta) + (alpha*log(N_dest[i]) + beta*log(N_orig[j]) - gamma*log(D[i,j])))
      )
      
    }
    trip.prop[i,] <- (trips[i,]/sum(trips[i,])) ## trip proportion. 
  }
  
  colnames(trip.prop) <- seq(1,16,1)
  rownames(trip.prop) <- seq(1,16,1)
  plots <- asymmetry.plots(trip.prop, title)
  return(plots)
}

## Plotting asymmetry map and trip count/prop map
asymmetry.plots <- function(trips, title){
  
  trips.ss <- skewsymmetry(trips); summary(trips) # Calculate skew-symmetry SSQ 
  
  # creates a color palette from red to blue
  my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 299)
  col_breaks = c(seq(-4000,-.001,length=100),  # negative values are red
                 seq(-.001,0.001,length=100),  # zeroes are white
                 seq(0.001,4000,length=100))   # positive values are blue
  
  # heatmap for asymmetric portion 
  png("hmap.png", width = 400, height = 500)
  hmap(trips.ss, col = my_palette, 
       xlab = "Destination", ylab = "Origin")
  dev.off()
  
  img1 <- readPNG("hmap.png")
  asym.map <- ggplot() + 
    background_image(img1)+ 
    ggtitle(title)
  
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


#### Generate trip details  ### Commented out for now - keeps randomly assigned populations the same between runs
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

## Run simulations and generate plots

#### 1. Basic trips - all parameters = 1 ####

basic_trips <- grav.model.basic(N_orig = dummy.trips$varied.pop,
                                N_dest = dummy.trips$varied.pop,
                                D = D,
                                theta = 1, 
                                alpha = 1, 
                                beta = 1, 
                                gamma = 1,
                                title = "Basic model, Varied Pop, all param = 1")
#### 2. Basic trips - all parameters = 1, except gamma = 2 ####
basic_trips_gamma2 <- grav.model.basic(N_orig = dummy.trips$varied.pop,
                                N_dest = dummy.trips$varied.pop,
                                D = D,
                                theta = 1, 
                                alpha = 1, 
                                beta = 1, 
                                gamma = 2,
                                title = "Basic model, Varied Pop, gamma = 2")



## Combine plots

asym.heat.maps <- ggarrange(basic_trips[[1]], basic_trips_gamma2[[1]],  
                            ncol = 2, nrow = 1, labels = c("A", "B"))

prop.OD.maps <-  ggarrange(basic_trips[[2]], basic_trips_gamma2[[2]],  
                           ncol = 2, nrow = 1, labels = c("A", "B"))
