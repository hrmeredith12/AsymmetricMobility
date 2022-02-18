## Script that calls functions that calculate different metrics and generate plots featured in Asymmetry paper
## Input: Dataframes with mobility details (trip count, date, origin, destination, etc)
##      - Trip data from CDRs were cleaned/generated in script "CDR_prep_function.R"
##      - Trip data generated from models were generated in script "Basic_gravity_model.R". Note that it takes 1-2 hours to run a basic gravity model. 
## Output: Entropy Index, Node Symmetry Index, Gini Coefficient, Asymmetry Transitive Index
## Last updated: February 18, 2022
## By: Hannah Meredith

library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(scales)

## Source functions
source("AsymmetryMetricFunctions.R")

## Import datafames
load("../../Data/NAM/NAM.daily.trips.RData") ## Daily trips calculated from CDRs
load("../../Data/NAM/GravityModelEstimates/NAM_adm2_daily_estimates_Basic_gm_dist_exp.RData") ## Avg daily trips estimated from gravity model
load("../../Data/NAM/GravityModelEstimates/NAM_adm2_daily_estimates_Basic_gm_dist_pwr.RData") ## Avg daily trips estimated from gravity model

load("../../Data/KEN/KEN.daily.trips.RData") ## Daily trips calculated from CDRs
load("../../Data/KEN/GravityModelEstimates/KEN_adm2_daily_estimates_Basic_gm_dist_exp.RData") ## Avg daily trips estimated from gravity model
load("../../Data/KEN/GravityModelEstimates/KEN_adm2_daily_estimates_Basic_gm_dist_pwr.RData") ## Avg daily trips estimated from gravity model


#### Calculate and compare Entropy Index ####

EI_NAM <- Entropy.Index(NAM, "Namibia", "CDRs")
EI_NAM_gm_basic_exp <- Entropy.Index(NAM.predicted.basic.exp, "Namibia", "Basic_exp")
EI_NAM_gm_basic_pwr <- Entropy.Index(NAM.predicted.basic.pwr, "Namibia", "Basic_pwr")
EI_KEN <- Entropy.Index(KEN, "Kenya", "CDRs")
EI_KEN_gm_basic_exp <- Entropy.Index(KEN.predicted.basic.exp, "Kenya", "Basic_exp")
EI_KEN_gm_basic_pwr <- Entropy.Index(KEN.predicted.basic.pwr, "Kenya", "Basic_pwr")

EI <- rbind(EI_NAM, EI_KEN)

ggplot(EI, aes(entropy.date, entropy.idx, color = as.factor(y)))+
  geom_point()+
  geom_line(aes(entropy.date, EI_mean), color = "black")+
  scale_x_date(labels = date_format("%b"))+
  facet_wrap(~location.name, ncol = 1)+
  labs(x = "Date", y = "Daily Entropy Index", color = "Year")+
  theme(panel.background = element_blank())

EI.table <- rbind.data.frame(EI_NAM[1, c("EI_mean", "EI_sd", "location.name")],   # Entropy index for Namibia: 0.469 +/- 0.0129
                             EI_NAM_gm_basic_exp[1, c("EI_mean", "EI_sd", "location.name")],
                             EI_NAM_gm_basic_pwr[1, c("EI_mean", "EI_sd", "location.name")],
                             EI_KEN[1, c("EI_mean", "EI_sd", "location.name")],   # Entropy index for Kenya: 0.4654 +/- 0.008
                             EI_KEN_gm_basic_exp[1, c("EI_mean", "EI_sd", "location.name")],
                             EI_KEN_gm_basic_pwr[1, c("EI_mean", "EI_sd", "location.name")])   

#### Calculate and compare Node Symmetry Index ####

NSI_NAM <- Node.symmetry.index(NAM, "Namibia", "CDRs")
NSI_NAM_gm_basic_exp <-Node.symmetry.index(NAM.predicted.basic.exp, "Namibia", "Basic_exp")
NSI_NAM_gm_basic_pwr <-Node.symmetry.index(NAM.predicted.basic.pwr, "Namibia", "Basic_pwr")
NSI_KEN <- Node.symmetry.index(KEN, "Kenya", "CDRs")
NSI_KEN_gm_basic_exp <-Node.symmetry.index(KEN.predicted.basic.exp, "Kenya", "Basic_exp")
NSI_KEN_gm_basic_pwr <-Node.symmetry.index(KEN.predicted.basic.pwr, "Kenya", "Basic_pwr")

NSI <- rbind(NSI_NAM, NSI_NAM_gm_basic_exp, NSI_NAM_gm_basic_pwr,
             NSI_KEN, NSI_KEN_gm_basic_exp, NSI_KEN_gm_basic_pwr)

ggplot(NSI, aes(NSI.yearly.avg, fill = trip.source))+
  geom_density(alpha = 0.5)+
  facet_wrap(~location.name, nrow = 2)+
  labs(x = "Averge Daily NSI", y = "Density", fill = "Trip count source")+
  lims(x = c(-0.05,0.05))+
  scale_y_sqrt()

NSI_NAM_map <- NAM.NSI.map(NSI_NAM)
NSI_NAM_gm_basic_exp_map <- NAM.NSI.map(NSI_NAM_gm_basic_exp)
NSI_KEN_map <- KEN.NSI.map(NSI_KEN)



NSI.table <- rbind.data.frame(NSI_NAM[1, c("NSI_mean", "NSI_sd", "location.name")],   # Entropy index for Namibia: 0.469 +/- 0.0129
                             NSI_KEN[1, c("NSI_mean", "NSI_sd", "location.name")])   # Entropy index for Kenya: 0.4654 +/- 0.008
