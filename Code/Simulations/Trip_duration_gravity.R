## Script to simulate effect of trip duration on asymmetry
## Based off of John R Gile's function "sim.gravity.duration" available on github
## Last updated Jan 4, 2022 by Hannah Meredith


##' EXtra information for simulating connectivity values using gravity model with trip duration
##'
##' This function simulates a connectivity matrix supplied model parameters in a gravity model formula that incorporates trip duration by using a conditional dispersal kernel (\eqn{f(d_ij | \lambda_ij)}) in 
##' the denominator. The gravity model still uses a Gamma distribution as the dispersal kernel, but this is scaled by the probability \eqn{Pr(\lambda_ij | d_ij)} according to Bayes theorem. If a 
##' vector of shape (\code{s}) and rate (\code{r}) parameters is supplied, the function will simulate route specific dispersal kernels based on the origin location (\eqn{i}). A null model (where all model parameters = 1) 
##' can be simulated by supplying only population sizes (\code{N}) and pairwise distances (\code{D}).
##' \deqn{
##'     \theta * ( N_i^\omega_1 N_j^\omega_2 / f(d_ij | \lambda_ij) )
##' }
##' 
##' @param N vector of population sizes
##' @param D matrix of distances among all \eqn{ij} pairs
##' @param theta scalar giving the proportionality constant of gravity formula (default = 1)
##' @param omega.1 scalar giving exponential scaling of origin population size (default = 1)
##' @param omega.2 scalar giving exponential scaling of destination population size (default = 1)
##' @param gamma scalar giving the dispersal kernel paramater (default = 1)
##' @param lambda matrix of trip duration decay parameters for each \eqn{ij} route
##' @param alpha model fitting parameter for the ECDF of lambda (default = 1)
##' @param counts logical indicating whether or not to return a count variable by scaling the connectivity matrix by origin population size (\eqn{N_i}) (default = FALSE)
##' 
##' @return a matrix with values between 0 and 1 (if \code{counts = FALSE}) or positive integers (if \code{counts = TRUE})
##' 
##' @author John Giles


## Load libraries
library(asymmetry)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(png)
library(reshape2)
library(dplyr)
library(DescTools)

### Functions


## This function estimates trips for each pair of locations. 
## It passes the trip matrix to another function "asymmetry.plots" and then returns the plots 
sim.gravity.duration <- function(N, D, theta, omega.1, omega.2, gamma, alpha, lambda, counts=FALSE, title) {
  
  if (!(identical(length(N), dim(D)[1], dim(D)[2], dim(lambda)[1], dim(lambda)[2]))) stop('Check dimensions of input data (N, D, or lambda)')
  if (!(length(c(theta, omega.1, omega.2, gamma)) == 4)) stop('theta, omega, and zeta parameters must be scalars')
  
  # Initialize simulation matrices
  n.districts <- length(N)
  diag(lambda) <- 0
  
  if (length(alpha) == 1) alpha <- rep(alpha, n.districts)
  
  x <- f.d <- f.d.lambda <- matrix(NA, n.districts, n.districts)
  
  # Simulate gravity model with duration
  for (i in 1:n.districts) {
    message(paste('Origin:', i, sep=' '))
    
    for (j in 1:n.districts) {
      
      # Gravity model
      if (i == j) {    # if considering "stays", set trip count to 0
        x[i,j] <- 0
      } else {
        f.d[i,j] <- (D[i,j]^gamma)
        f.d.lambda[i,j] <- (f.d[i,j] * (1 - sum(lambda[,] <= lambda[i,j])/(n.districts^2))^alpha[i]) + 1e-6 # Conditional dispersal kernel
        x[i,j] <- exp(log(theta) + (omega.1*log(N[i]) + omega.2*log(N[j]) - log(f.d.lambda[i,j])))
      }   
    }
    
    x[i,] <- (x[i,]/sum(x[i,]))
    if (counts == TRUE) x[i,] <- round(x[i,]*N[i])
  }
  
  dimnames(x) <- list(origin=dimnames(D)[[1]], destination=dimnames(D)[[2]])
  colnames(x) <- rownames(x) <- colnames(D)
  
  plots <- asymmetry.plots(x, title)
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
  hmap(trips.ss, col = my_palette, 
       xlab = "Destination", ylab = "Origin")
  dev.off()
  
  img1 <- readPNG("hmap.png")
  asym.map <- ggplot() + 
    background_image(img1) +
    # This ensures that the image leaves some space at the edges
    # theme(plot.margin = margin(t=1, l=1, r=1, b=1, unit = "cm"))+
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

## Matrices of different lambda values (trip duration)
lambda_0.1 <- matrix(0.1, nrow=dim(dummy.trips)[1], ncol=dim(dummy.trips)[1]) # mean trip duration = 10 days; trip duration decay rate = 0.1
lambda_1 <- matrix(1, nrow=dim(dummy.trips)[1], ncol=dim(dummy.trips)[1]) # mean trip duration = 1 day; trip duration decay rate = 0
lambda_2mixed <- lambda_1;  lambda_2mixed[1,2] <- lambda_2mixed[8,2] <- 0.1 ## all locations have mean trip duration of 1 day, except two routes with trip dur = 10 days
lambda_mixed <- matrix(runif(dim(dummy.trips)[1]*dim(dummy.trips)[1], 1/10, 1), nrow=dim(dummy.trips)[1], ncol=dim(dummy.trips)[1]) # Random selection of trip duration between 1-10 days for all trips
lambda_row <- lambda_1; lambda_row[5,] <- 0.1 # all trips made from a certain location have trip duration = 10 days
lambda_col <- lambda_1; lambda_col[,5] <- 0.1 # all trips made to a certain location have trip duration = 10 days
lambda_longDist <- lambda_1 * ifelse(D < quantile(D)[2], 1,
                          ifelse(D > quantile(D)[2] & D < quantile(D)[3], 0.5, 
                          ifelse(D > quantile(D)[3] & D < quantile(D)[4], 0.2, 
                                 0.1)))
lambda_shortDist <- lambda_1 * ifelse(D < quantile(D)[2], 0.1,
                                     ifelse(D > quantile(D)[2] & D < quantile(D)[3], 0.2, 
                                            ifelse(D > quantile(D)[3] & D < quantile(D)[4], 0.5, 
                                                   1)))
lambda_largePop<- lambda_1 * ifelse(dummy.trips$varied.pop < quantile(dummy.trips$varied.pop)[2], 1,  # 1-2 days
                                     ifelse(dummy.trips$varied.pop > quantile(dummy.trips$varied.pop)[2] & dummy.trips$varied.pop < quantile(dummy.trips$varied.pop)[3], 1/2, # 2-3 days
                                            ifelse(dummy.trips$varied.pop > quantile(dummy.trips$varied.pop)[3] & dummy.trips$varied.pop < quantile(dummy.trips$varied.pop)[4], 1/5, # 3-5 days 
                                                   1/10))) # 5-10 days
lambda_smallPop <- lambda_1 * ifelse(dummy.trips$varied.pop < quantile(dummy.trips$varied.pop)[2], 1/10,
                                      ifelse(dummy.trips$varied.pop > quantile(dummy.trips$varied.pop)[2] & dummy.trips$varied.pop < quantile(dummy.trips$varied.pop)[3], 1/5, 
                                             ifelse(dummy.trips$varied.pop > quantile(dummy.trips$varied.pop)[3] & dummy.trips$varied.pop < quantile(dummy.trips$varied.pop)[4], 1/2, 
                                                    1)))

## Define gravity model parameters

## Estimate trip counts with gravity model for different trip duration scenarios
N = dummy.trips$varied.pop 
D = D
theta=1
omega.1=1
omega.2=1
gamma=1
alpha=1

## Run simulations and generate plots

# plots_lambda_0.1 <- sim.gravity.duration(N, D, theta, omega.1, omega.2, gamma, alpha, lambda = lambda_0.1, counts = FALSE, title = "    All trip durations = 10 days \n")
# plots_lambda_1 <- sim.gravity.duration(N, D, theta, omega.1, omega.2, gamma, alpha, lambda = lambda_1, counts = FALSE, title = "    All trip durations = 1 day")
plots_lambda_2mixed <- sim.gravity.duration(N, D, theta, omega.1, omega.2, gamma, alpha, lambda = lambda_2mixed, counts = FALSE, title = "    All trip durations = 1 day, \nexcept for two routes = 10 days")
plots_lambda_mixed <- sim.gravity.duration(N, D, theta, omega.1, omega.2, gamma, alpha, lambda = lambda_mixed, counts = FALSE, title = "    Random trip durations")
plots_lambda_row <- sim.gravity.duration(N, D, theta, omega.1, omega.2, gamma, alpha, lambda = lambda_row, counts = FALSE, title ="    All trip durations = 1 day, \nexcept 1 origin (#5) that has 10 day trips")
plots_lambda_col <- sim.gravity.duration(N, D, theta, omega.1, omega.2, gamma, alpha, lambda = lambda_col, counts = FALSE, title = "    All trip durations = 1 day, \nexcept 1 destination (#5) that has 10 day trips")
plots_lambda_longDist <- sim.gravity.duration(N, D, theta, omega.1, omega.2, gamma, alpha, lambda = lambda_longDist, counts = FALSE, title = "    Trip duration increases with distance")
plots_lambda_shortDist <- sim.gravity.duration(N, D, theta, omega.1, omega.2, gamma, alpha, lambda = lambda_shortDist, counts = FALSE, title = "    Trip duration decreases with distance")
plots_lambda_largePop <- sim.gravity.duration(N, D, theta, omega.1, omega.2, gamma, alpha, lambda = lambda_largePop, counts = FALSE, title = "    Trip duration increases with population")
plots_lambda_smallPop <- sim.gravity.duration(N, D, theta, omega.1, omega.2, gamma, alpha, lambda = lambda_smallPop, counts = FALSE, title = "    Trip duration decreases with population")


asym.heat.maps <- ggarrange(plots_lambda_2mixed[[1]], plots_lambda_mixed[[1]], plots_lambda_largePop[[1]], plots_lambda_smallPop[[1]], plots_lambda_row[[1]],
                            plots_lambda_col[[1]], plots_lambda_1[[1]], plots_lambda_longDist[[1]], plots_lambda_shortDist[[1]],  
                            ncol = 5, nrow = 2, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"))

prop.OD.maps <- ggarrange(plots_lambda_2mixed[[2]], plots_lambda_mixed[[2]], plots_lambda_largePop[[2]], plots_lambda_smallPop[[2]], plots_lambda_row[[2]],
                          plots_lambda_col[[2]], plots_lambda_1[[2]], plots_lambda_longDist[[2]], plots_lambda_shortDist[[2]],  
                          ncol = 5, nrow = 2, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"))
