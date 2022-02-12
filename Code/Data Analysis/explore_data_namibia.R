###############################################################################
# Ronan Corgel
# Asymmetry Project
# Explore Namibia Data
# February 2022
# Steps:
# 1. Load data from previous files
# 2. Visualize matrices and calculate Gini coefficients
# 2.1. Create heat maps and calculate Gini for those heat maps
# 2.2. Calculate Gini by origin and plot distributions
# 2.3. Calclate and visualize distribution for Gini by route
# 2.4. Break up data by month, calculate Gini coefficients over time, and visualize
###############################################################################

# Remove objects previously stored in R
rm(list = ls())

# Set Directories
setwd('/Users/rcorgel/OneDrive/Projects/asymmetry-project') 

# Load libraries
library('ggplot2')
library('asymmetry')
library('RColorBrewer')
library('dplyr')
library('DescTools')
library('reshape')

# 1. Load data from previous files
load('output/NAM_overall_prop_ODMatrix.RData')
load('output/NAM_monthly_prop_ODMatrix.RData')
load('output/NAM_weekly_prop_ODMatrix.RData')
load('output/NAM_overall_trips_ODMatrix.RData')
load('output/NAM_monthly_trips_ODMatrix.RData')
load('output/NAM_weekly_trips_ODMatrix.RData')
load('tmp/NAM_trips_data_long.RData')
load('tmp/adm2.trip.month.RData')

# 2. Visualize matrices and calculate Gini coefficients 
# 2.1. Create heat maps and calculate Gini for those heat maps
SS.monthly.avg <- skewsymmetry(M.monthly.avg)
SS.weekly.avg <- skewsymmetry(M.weekly.avg)
SS.overall.trips <- skewsymmetry(M.overall.trips)
SS.monthly.avg.prop <- skewsymmetry(M.monthly.avg.prop)
SS.weekly.avg.prop <- skewsymmetry(M.weekly.avg.prop)
SS.overall.prop <- skewsymmetry(M.overall.prop)

# creates a color palette from red to blue
my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 299)
col_breaks = c(seq(-50000,-.001,length=100),  # negative values are red
               seq(-.001,0.01,length=100),   # zeroes are white
               seq(0.01,50000,length=100))  # positive values are blue

# Make heat maps of asymmetric matrices
pdf(file='visuals/hmap_month_nam.pdf')
hmap(SS.monthly.avg, col = my_palette, ylab = "origin", xlab = "destination", main = "Monthly Average")
dev.off()
pdf(file='visuals/hmap_week_nam.pdf')
hmap(SS.weekly.avg, col = my_palette, ylab = "origin", xlab = "destination", main = "Weekly Average")
dev.off()
pdf(file='visuals/hmap_overall_nam.pdf')
hmap(SS.overall.trips, col = my_palette, ylab = "origin", xlab = "destination", main = "Overall Trips")
dev.off()
pdf(file='visuals/hmap_month_prop_nam.pdf')
hmap(SS.monthly.avg.prop, col = my_palette, ylab = "origin", xlab = "destination", main = "Monthly Average Proportion")
dev.off()
pdf(file='visuals/hmap_week_prop_nam.pdf')
hmap(SS.weekly.avg.prop, col = my_palette, ylab = "origin", xlab = "destination", main = "Weekly Average Proportion")
dev.off()
pdf(file='visuals/hmap_overall_prop_nam.pdf')
hmap(SS.overall.prop, col = my_palette, ylab = "origin", xlab = "destination", main = "Overall Proportion")
dev.off()

# Overall Gini
Gini(abs(SS.monthly.avg[["A"]]))
Gini(abs(SS.weekly.avg[["A"]]))
Gini(abs(SS.overall.trips[["A"]]))
Gini(abs(SS.monthly.avg.prop[["A"]]))
Gini(abs(SS.weekly.avg.prop[["A"]]))
Gini(abs(SS.overall.prop[["A"]]))

# 2.2. Calculate Gini by origin and plot distributions
# Origin Gini
A.monthly.avg <- as.data.frame(SS.monthly.avg[["A"]])
A.monthly.avg$gini.month.avg <- apply(A.monthly.avg, 1, function(x) Gini(abs(x)))
A.monthly.avg <- A.monthly.avg %>% mutate(id = 1:n())

A.weekly.avg <- as.data.frame(SS.weekly.avg[["A"]])
A.weekly.avg$gini.week.avg <- apply(A.weekly.avg, 1, function(x) Gini(abs(x)))
A.weekly.avg <- A.weekly.avg %>% mutate(id = 1:n())

A.overall.trips <- as.data.frame(SS.overall.trips[["A"]])
A.overall.trips$gini.overall.avg <- apply(A.overall.trips, 1, function(x) Gini(abs(x)))
A.overall.trips <- A.overall.trips %>% mutate(id = 1:n())

A.monthly.prop <- as.data.frame(SS.monthly.avg.prop[["A"]])
A.monthly.prop$gini.month.prop <- apply(A.monthly.prop, 1, function(x) Gini(abs(x)))
A.monthly.prop <- A.monthly.prop %>% mutate(id = 1:n())

A.weekly.prop <- as.data.frame(SS.weekly.avg.prop[["A"]])
A.weekly.prop$gini.week.prop <- apply(A.weekly.prop, 1, function(x) Gini(abs(x)))
A.weekly.prop <- A.weekly.prop %>% mutate(id = 1:n())

A.overall.prop <- as.data.frame(SS.overall.prop[["A"]])
A.overall.prop$gini.overall.prop <- apply(A.overall.prop, 1, function(x) Gini(abs(x)))
A.overall.prop <- A.overall.prop %>% mutate(id = 1:n())

# Make single data frame for visualization of origin Gini distributions
gini.origin.month <- A.monthly.avg %>% select('id', 'gini.month.avg')
gini.origin.month <- left_join(gini.origin.month, A.monthly.prop[, c('id', 'gini.month.prop')], by = 'id')
gini.origin.month.long <- reshape2::melt(gini.origin.month, id.vars = c("id"))

gini.origin.week <- A.weekly.avg[, c('id', 'gini.week.avg')]
gini.origin.week <- left_join(gini.origin.week, A.weekly.prop[, c('id', 'gini.week.prop')], by = 'id')
gini.origin.week.long <- reshape2::melt(gini.origin.week, id.vars = c("id"))

gini.origin.overall <- A.overall.trips[, c('id', 'gini.overall.avg')]
gini.origin.overall <- left_join(gini.origin.overall, A.overall.prop[, c('id', 'gini.overall.prop')], by = 'id')
gini.origin.overall <-  dplyr::rename(gini.origin.overall, gini.overall.trips = gini.overall.avg)
gini.origin.overall.long <- reshape2::melt(gini.origin.overall, id.vars = c("id"))

# Visualize origin Gini distributions
ggplot(gini.origin.month.long, aes(x=value, fill = variable)) + 
  geom_density(alpha=.2) + theme_bw() + labs(title = 'Origin Gini Distribution (Month) \nAverage Trip vs. Trip Proportion')
ggsave('visuals/month_dist_origin.png')  

ggplot(gini.origin.week.long, aes(x=value, fill = variable)) + 
  geom_density(alpha=.2) + theme_bw() + labs(title = 'Origin Gini Distribution (Week) \nAverage Trip vs. Trip Proportion')
ggsave('visuals/week_dist_origin.png')  

ggplot(gini.origin.overall.long, aes(x=value, fill = variable)) + 
  geom_density(alpha=.2) + theme_bw() + labs(title = 'Origin Gini Distribution (Overall) \nTrips vs. Trip Proportion')
ggsave('visuals/overall_dist_origin.png')  

# Compare by level of aggregation
gini.origin.avg <- A.monthly.avg %>% select('id', 'gini.month.avg')
gini.origin.avg <- left_join(gini.origin.avg, A.weekly.avg[, c('id', 'gini.week.avg')], by = 'id')
gini.origin.avg <- left_join(gini.origin.avg, A.overall.trips[, c('id', 'gini.overall.avg')], by = 'id')
gini.origin.avg <-  dplyr::rename(gini.origin.avg, gini.overall.trips = gini.overall.avg)
gini.origin.avg.long <- reshape2::melt(gini.origin.avg, id.vars = c("id"))

gini.origin.prop <- A.monthly.prop %>% select('id', 'gini.month.prop')
gini.origin.prop <- left_join(gini.origin.prop, A.weekly.prop[, c('id', 'gini.week.prop')], by = 'id')
gini.origin.prop <- left_join(gini.origin.prop, A.overall.prop[, c('id', 'gini.overall.prop')], by = 'id')
gini.origin.prop.long <- reshape2::melt(gini.origin.prop, id.vars = c("id"))

# Visualize origin Gini distributions
ggplot(gini.origin.avg.long, aes(x=value, fill = variable)) + 
  geom_density(alpha=.2) + theme_bw() + labs(title = 'Origin Gini Distribution (Average Trip) \nOverall vs. Month vs. Week')
ggsave('visuals/level_dist_origin_avg.png')  

ggplot(gini.origin.prop.long, aes(x=value, fill = variable)) + 
  geom_density(alpha=.2) + theme_bw() + labs(title = 'Origin Gini Distribution (Average Trip Proportion) \nOverall vs. Month vs. Week')
ggsave('visuals/level_dist_origin_prop.png')  

# 2.3. Calclate and visualize distribution for Gini by route
# Calculate the total trips out of each origin each day
trip.data.long <- trip.data.long %>% group_by(start.adm2.code, date) %>%
  mutate(trip.total = ifelse(start.adm2.code == end.adm2.code, 0, 
                             sum(trip.count[which(start.adm2.code != end.adm2.code)], na.rm = T)))

# Calculate trip proportion of a specific trip between a given origin and destination
trip.data.long$trip.prop = trip.data.long$trip.count/trip.data.long$trip.total  

# Calculate trips in, since we currently have trips out
trip.data.long.in <- trip.data.long[, c('start.adm2.name', 'end.adm2.name', 'trip.count', 'trip.prop', 'date')]
trip.data.long.in <- trip.data.long.in %>% dplyr::rename(start.adm2.name = end.adm2.name, 
                                                         end.adm2.name = start.adm2.name,
                                                         trip.count.in = trip.count,
                                                         trip.prop.in = trip.prop)
# Merge trips in to main data
trip.data.long <- dplyr::left_join(trip.data.long, trip.data.long.in, by = c('start.adm2.name', 'end.adm2.name', 'date'))
# Make daily net trip variable
trip.data.long$net.trips <- (trip.data.long$trip.count - trip.data.long$trip.count.in) / 2
trip.data.long$net.trips.prop <- (trip.data.long$trip.prop - trip.data.long$trip.prop.in) / 2
# Make pair ID variable
trip.data.long.pair <- trip.data.long
trip.data.long.pair <- trip.data.long.pair %>%
  mutate(ordered.pair = paste(pmin(start.adm2.code, end.adm2.code), pmax(start.adm2.code, end.adm2.code), sep = "::")) %>%
  distinct(date, ordered.pair, .keep_all = TRUE)
# Make name pair variable
trip.data.long.pair <- trip.data.long.pair %>% mutate(ordered.pair.name = paste(start.adm2.name, end.adm2.name, sep = "::"))

# Gini within Route
trip.data.long.pair <- trip.data.long.pair %>% group_by(start.adm2.name, end.adm2.name) %>%
  mutate(gini.net.trips = Gini(abs(net.trips)),
         gini.net.trips.prop = Gini(abs(net.trips.prop)))

gini.by.route <- trip.data.long.pair %>% distinct(start.adm2.name, end.adm2.name, .keep_all = TRUE)

gini.by.route.sub <- gini.by.route[, c('ordered.pair', 'gini.net.trips', 'gini.net.trips.prop')]
gini.by.route.long <- reshape2::melt(gini.by.route.sub, id.vars = c('ordered.pair'))

# Visualize the distribution
ggplot(data=gini.by.route.long, aes(x= value, fill = variable)) +
  geom_density(alpha = 0.2) + theme_bw() + labs(title = 'Gini by Route \nNet Trips vs. Proportion of Net Trips')
ggsave('visuals/gini_route_dist.png')  

# 2.4. Break up data by month, calculate Gini coefficients over time, and visualize
# Trips
overall.gini.month <- c()
y <- c()
m <- c()
d <- c()
origin.gini <- data.frame(gini= numeric(0), id = numeric(0), month.date = as.Date(character()))
for (i in 2010:2014) {
  for (j in 1:12) {
    sub <- adm2.trip.month %>% filter(y == i & m == j)
    if (dim(sub)[1] == 0) {
      print('skip, empty dataset')
    }
    else {
      trip.month.summary <- sub[,c('start.adm2.code', 'end.adm2.code', 'adm2.single.trip.sum')]
      # Reshape to wide
      M.monthly <- reshape::cast(trip.month.summary, start.adm2.code ~ end.adm2.code)            
      # Label rows with district numbers
      rownames(M.monthly) <- M.monthly$start.adm2.code                           
      # Get rid of the first column
      M.monthly <- M.monthly[ ,-1]
      class(M.monthly) <- 'data.frame'
      M.monthly <- as.matrix(M.monthly)
      names(dimnames(M.monthly)) <- c('origin', 'destination')
      # Replace NAs with 0
      M.monthly[is.na(M.monthly)] <- 0
      assign(paste("M.month",i, j, sep = "."), M.monthly)
      SS.monthly <- skewsymmetry(M.monthly)
      assign(paste("SS.month",i, j, sep = "."), SS.monthly)
      A.monthly <- as.data.frame(SS.monthly[["A"]])
      A.monthly$gini <- apply(A.monthly, 1, function(x) Gini(abs(x)))
      assign(paste("A.month",i, j, sep = "."), A.monthly)
      A.monthly <- A.monthly %>% mutate(id = as.character(rownames(A.monthly)))
      subset <- A.monthly[,c('id', 'gini')]
      subset$month.date <- as.Date(paste(i, j, 1,sep="-"), "%Y-%m-%d")
      origin.gini <- rbind(origin.gini, subset)
      d <- c(d, 1)
      y <- c(y, i)
      m <- c(m, j)
      overall.gini.month <- c(overall.gini.month, Gini(abs(SS.monthly[["A"]])))
    }
  }
}

overall.gini <- data.frame(y, m, d, overall.gini.month)
overall.gini$month.date <- as.Date(with(overall.gini, paste(y,m,d,sep="-")), "%Y-%m-%d")

# Proportions
overall.gini.month.p <- c()
y.p <- c()
m.p <- c()
d.p <- c()
origin.gini.p <- data.frame(gini= numeric(0), id = numeric(0), month.date = as.Date(character()))
for (i in 2010:2014) {
  for (j in 1:12) {
    sub.p <- adm2.trip.month %>% filter(y == i & m == j)
    if (dim(sub.p)[1] == 0) {
      print('skip, empty dataset')
    }
    else {
      trip.month.summary.p <- sub.p[,c('start.adm2.code', 'end.adm2.code', 'monthly.trip.proportion')]
      # Reshape to wide
      M.monthly.p <- reshape::cast(trip.month.summary.p, start.adm2.code ~ end.adm2.code)            
      # Label rows with district numbers
      rownames(M.monthly.p) <- M.monthly.p$start.adm2.code                           
      # Get rid of the first column
      M.monthly.p <- M.monthly.p[ ,-1]
      class(M.monthly.p) <- 'data.frame'
      M.monthly.p <- as.matrix(M.monthly.p)
      names(dimnames(M.monthly.p)) <- c('origin', 'destination')
      # Replace NAs with 0
      M.monthly.p[is.na(M.monthly.p)] <- 0
      M.monthly.p[is.infinite(M.monthly.p)] <- 0
      assign(paste("M.month.p",i, j, sep = "."), M.monthly.p)
      SS.monthly.p <- skewsymmetry(M.monthly.p)
      assign(paste("SS.month.p",i, j, sep = "."), SS.monthly.p)
      A.monthly.p <- as.data.frame(SS.monthly.p[["A"]])
      A.monthly.p$gini <- apply(A.monthly.p, 1, function(x) Gini(abs(x)))
      assign(paste("A.month.p",i, j, sep = "."), A.monthly.p)
      A.monthly.p <- A.monthly.p %>% mutate(id = as.character(rownames(A.monthly.p)))
      subset.p <- A.monthly.p[,c('id', 'gini')]
      subset.p$month.date <- as.Date(paste(i, j, 1,sep="-"), "%Y-%m-%d")
      origin.gini.p <- rbind(origin.gini.p, subset.p)
      d.p <- c(d.p, 1)
      y.p <- c(y.p, i)
      m.p <- c(m.p, j)
      overall.gini.month.p <- c(overall.gini.month.p, Gini(abs(SS.monthly.p[["A"]])))
    }
  }
}

overall.gini.p <- data.frame(y.p, m.p, d.p, overall.gini.month.p)
overall.gini.p$month.date <- as.Date(with(overall.gini.p, paste(y.p,m.p,d.p,sep="-")), "%Y-%m-%d")

# Visualize Gini Over Time
# Overall Trips
ggplot(overall.gini, aes(x = month.date, y = overall.gini.month)) +
  geom_line() +
  theme_bw() +
  labs(title = "Gini Over Time by Month", y = "Overall Gini Coefficient", x = "Month") +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave('visuals/gini_month_over_time.png')

# Overall Proportions
ggplot(overall.gini.p, aes(x = month.date, y = overall.gini.month.p)) +
  geom_line() +
  theme_bw() +
  labs(title = "Gini Over Time by Month (Proportion)", y = "Overall Gini Coefficient", x = "Month") +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave('visuals/gini_month_over_time_prop.png')

# Origin Gini for Trips
origin.gini <- origin.gini %>% group_by(id) %>%
  mutate(gini.percent.change = (gini - mean(gini)) / mean(gini))

ggplot(origin.gini, aes(x=month.date, y=id, fill=gini)) + 
  geom_tile(col = 'white', size = 0) + scale_fill_gradient2(low="red", mid="white", high="blue", midpoint = median(origin.gini$gini)) +
  scale_x_date(date_breaks = "2 months", date_labels = "%m-%d-%y") +
  xlab('') + ylab('Origin Code') + 
  labs(title = 'Gini Over Time by Month \nOrigin Trips') +
  theme(legend.position='bottom', axis.text.y=element_text(size = 4), axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),)  +
  guides(fill = guide_colourbar(title = 'Gini', barwidth = 10, barheight = 0.5, direction = "horizontal", nbin = 10)) 
ggsave('visuals/gini_month_over_time_origin.png')

# Origin Gini for Trips % Change
ggplot(origin.gini, aes(x=month.date, y=id, fill=gini.percent.change)) + 
  geom_tile(col = 'white', size = 0) + scale_fill_gradient2(low="red", mid="white", high="blue", midpoint = 0) +
  scale_x_date(date_breaks = "2 months", date_labels = "%m-%d-%y") +
  xlab('') + ylab('Origin Code') + 
  labs(title = 'Gini Over Time by Month \nOrigin Trips, Percent Change within Origin') +
  theme(legend.position='bottom', axis.text.y=element_text(size = 4), axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),)  +
  guides(fill = guide_colourbar(title = 'Gini % Change', barwidth = 10, barheight = 0.5, direction = "horizontal", nbin = 10)) 
ggsave('visuals/gini_month_over_time_origin_perc.png')

# Origin Gini for Proportions
origin.gini.p <- origin.gini.p %>% group_by(id) %>%
  mutate(gini.percent.change = (gini - mean(gini)) / mean(gini))

ggplot(origin.gini.p, aes(x=month.date, y=id, fill=gini)) + 
  geom_tile(col = 'white', size = 0) + scale_fill_gradient2(low="red", mid="white", high="blue", midpoint = median(origin.gini.p$gini)) +
  scale_x_date(date_breaks = "2 months", date_labels = "%m-%d-%y") +
  xlab('') + ylab('Origin Code') + 
  labs(title = 'Gini Over Time by Month (Proportions) \nOrigin Trips') +
  theme(legend.position='bottom', axis.text.y=element_text(size = 4), axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),)  +
  guides(fill = guide_colourbar(title = 'Gini', barwidth = 10, barheight = 0.5, direction = "horizontal", nbin = 10)) 
ggsave('visuals/gini_month_over_time_origin_prop.png')

# Origin Gini for Proportions % Change
ggplot(origin.gini.p, aes(x=month.date, y=id, fill=gini.percent.change)) + 
  geom_tile(col = 'white', size = 0) + scale_fill_gradient2(low="red", mid="white", high="blue", midpoint = 0) +
  scale_x_date(date_breaks = "2 months", date_labels = "%m-%d-%y") +
  xlab('') + ylab('Origin Code') + 
  labs(title = 'Gini Over Time by Month (Proportions) \nOrigin Trips, Percent Change within Origin') +
  theme(legend.position='bottom', axis.text.y=element_text(size = 4), axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),)  +
  guides(fill = guide_colourbar(title = 'Gini % Change', barwidth = 10, barheight = 0.5, direction = "horizontal", nbin = 10)) 
ggsave('visuals/gini_month_over_time_origin_perc_prop.png')
