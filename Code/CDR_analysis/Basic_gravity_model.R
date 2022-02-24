## Basic gravity model functions
## Comment out/uncomment lines relevant to the model being run (exponential or power form of distance kernel)

library(hmob)
library(mobility)
library(DescTools)

fit_basic_gravity <- function (M, D, N = NULL, N_orig = NULL, N_dest = NULL, n_chain = 2,
                                   n_burn = 1000, n_samp = 1000, n_thin = 1, prior = NULL, DIC = FALSE,
                                   parallel = FALSE)
{
  if (all(!is.null(N), is.null(N_orig), is.null(N_dest))) {
    N_dest <- N_orig <- N
  }
  else if (all(is.null(N), !is.null(N_orig), is.null(N_dest))) {
    N_dest <- N_orig
  }
  if (!(identical(dim(M)[1], dim(D)[1], length(N_orig))))
    stop("Dimensions of input data must match")
  if (!(identical(dim(M)[2], dim(D)[2], length(N_dest))))
    stop("Dimensions of input data must match")
  if (!(identical(dimnames(M)[[1]], dimnames(D)[[1]])) | !(identical(dimnames(M)[[1]],
                                                                     names(N_orig)))) {
    stop("Dimension names of input data do not match")
  }
  if (!(identical(dimnames(M)[[2]], dimnames(D)[[2]])) | !(identical(dimnames(M)[[2]],
                                                                     names(N_dest)))) {
    stop("Dimension names of input data do not match")
  }
  message(paste("::Fitting gravity model for", dim(M)[1],
                "origins and", dim(M)[2], "destinations::",
                sep = " "))
  if (!all(unlist(lapply(list(M, N_orig, N_dest), is.integer)))) {
    M[, ] <- as.integer(M)
    N_orig[] <- as.integer(N_orig)
    N_dest[] <- as.integer(N_dest)
  }
  diag(M) <- 0
  if (is.null(prior)) {
    message("Using uniformative priors")
    null_prior <- c(1, 0.5)
    prior <- list(theta = null_prior, omega_1 = null_prior,
                  omega_2 = null_prior, gamma = null_prior)
  }
  else {
    message("Using supplied informative priors")
  }
  jags_data <- list(M = M, 
                    D = D, 
                    N_orig = N_orig, 
                    N_dest = N_dest,
                    prior_theta = prior$theta, 
                    prior_omega_1 = prior$omega_1,
                    prior_omega_2 = prior$omega_2, 
                    prior_gamma = prior$gamma)
  
  jags_model <- "  
    model {
      
      # Poisson likelihood
      for (i in 1:length(N_orig)) {
        for (j in 1:length(N_dest)) {
          M[i,j] ~ dpois(c[i,j])  
        }
      }
      
      # Gravity model
      for (i in 1:length(N_orig)) {
        for (j in 1:length(N_dest)) {
          c[i,j] <- ifelse(
            i == j,
            1e-06,
            exp(log(theta) + (omega_1*log(N_dest[i]) + omega_2*log(N_orig[j]) - log(f_d[i,j] )))
          )
         
          # Exponential distance model:
          f_d[i,j] <- exp(D[i,j]/gamma)

          # # Power distance model:
          # f_d[i,j] <- D[i,j]^gamma
        }
      }
      
      # Priors
      theta ~ dgamma(prior_theta[1], prior_theta[2])
      omega_1 ~ dgamma(prior_omega_1[1], prior_omega_1[2])
      omega_2 ~ dgamma(prior_omega_2[1], prior_omega_2[2])
      gamma ~ dgamma(prior_gamma[1], prior_gamma[2])
      
    }"
  params <- c("omega_1", "omega_2", "theta", "gamma")
  fit_jags(jags_data = jags_data, jags_model = jags_model,
           params = params, n_chain = n_chain, n_burn = n_burn,
           n_samp = n_samp, n_thin = n_thin, DIC = DIC, parallel = parallel)
}

estimate_trips <- function (Key_code_name, dist.fx, M.observed, N, D, theta = 1, omega.1 = 1, omega.2 = 1, gam = 1, model, district.ID, district.label = rank.Id, location.name, counts = FALSE) 
  {
  if (!(identical(length(N), dim(D)[1], dim(D)[2])))
    stop("Check dimensions of input data N and D")
  if (!(length(c(theta, omega.1, omega.2)) == 3))
    stop("theta and omega parameters must be scalars")
  n.districts <- length(N)
  x <- f.d <- matrix(NA, n.districts, n.districts)
  error.fit.origin <- matrix(NA, n.districts, 4)
  
  for (i in 1:n.districts) {
    for (j in 1:n.districts) {
      if (i == j) {
        x[i, j] <- 0
      }
      else {
        if (length(gam) == 1 & dist.fx == "pwr") {
          f.d[i, j] <- D[i,j]^gam
        }
        else if (length(gam) == n.districts & dist.fx == "pwr") {
          f.d[i, j] <- D[i,j]^gam[i]
        }
        else if (length(gam) == 1 & dist.fx == "exp") {
          f.d[i, j] <- exp(D[i,j]/gam)
        }
        else if (length(gam) == n.districts & dist.fx == "exp") {
          f.d[i, j] <- exp(D[i,j]/gam[i])
        }
        else {
          stop("Incorrect length for parameters in Gamma dispersal kernel(s)")
        }
        x[i, j] <- exp(log(theta) + (omega.1 * log(N[i]) + omega.2 * log(N[j]) - log(f.d[i, j])))
      }
    }
  }
  
  colnames(x) <- colnames(D)
  rownames(x) <- rownames(D)
  
  M.predicted <- reshape2::melt(x, 
                                varnames = c("start.adm2.code", "end.adm2.code"),
                                value.name = "trip.count")
  
  M.predicted$y <- rep("2010", nrow = nrow(M.predicted))
  M.predicted$m <- rep("01", nrow = nrow(M.predicted))
  M.predicted$d <- rep("01", nrow = nrow(M.predicted))
  M.predicted$trip.date <- as.Date(paste("2010", M.predicted$m, M.predicted$d, sep="-"), "%Y-%m-%d")
  M.predicted$location.name <- rep(location.name, nrow = nrow(M.predicted))
  M.predicted$trip.source <- rep(model, nrow = nrow(M.predicted))
  M.predicted <- left_join(M.predicted, Key_code_name, by = c("start.adm2.code"))
  M.predicted <- left_join(M.predicted, Key_code_name, by = c("end.adm2.code" = "start.adm2.code"))
  colnames(M.predicted) <- c("start.adm2.code", "end.adm2.code", "trip.count", "y", "m", "d", "trip.date", "location.name",
                             "trip.source",  "start.adm1.name", "start.adm2.name", "start.adm1.code", "pop.start", 'urb.start',
                             "end.adm1.name", "end.adm2.name", "end.adm1.code", 'pop.end', 'urb.end') 
  M.predicted <- M.predicted[, c("start.adm1.name", "start.adm2.name", "start.adm1.code", "start.adm2.code", "end.adm1.name", "end.adm2.name", "end.adm1.code", "end.adm2.code",
                                 "d", "m", "y", "trip.date", "trip.count","pop.start", 'urb.start', 'pop.end', 'urb.end')]
  
  return(M.predicted)
}


### Common parameters for all runs: 
n_chain = 4
n_burn = 4000
n_samp = 10000
n_thin = 5
prior = NULL
DIC = TRUE
parallel = TRUE

t.start = Sys.time(); print(Sys.time())

# ### Namibia adm2 ####
M <- load.obj(1, '../../Data/NAM/NAM.daily.avg.matrix.RData')          	# trip counts
D <- load.obj(1, '../../Data/NAM/NAM_adm2Joined_ppp_distance.RData')    # distance matrix
D <- D[!rownames(D) %in% c("82","84"),!colnames(D) %in% c("82","84")]   # distance matrix, remove districts with no travel recorded
N <- load.obj(1, '../../Data/NAM/NAM_adm2Joined_ppp_pop2010.RData')     # population vector
rownames(N) <- N$ID_2
N <- N[!N$ID_2 %in% c("82","84"), ][2]
N <- t(N) ;
N <- as.vector(N)
names(N) <- colnames(D)
# 
# print(n_chain)
# print(n_burn)
# print(n_samp)
# print(n_thin)
# # Call jags to fit gravity model
# out.NAM <- fit_basic_gravity(M,
#                              D,
#                              N,
#                              n_chain = n_chain,
#                              n_burn = n_burn,
#                              n_samp = n_samp,
#                              n_thin = n_thin,
#                              prior = NULL,
#                              DIC = TRUE,
#                              parallel = TRUE)
# 
# # save(out.NAM, file='../../Data/NAM/GravityModelEstimates/NAM_adm2_daily_estimates_Basic_gm_dist_pwr_chains.RData')
# print('Saved: ../../Data/NAM/GravityModelEstimates/NAM_adm2_daily_estimates_Basic_gm_dist_pwr_chains.RData')
# print(Sys.time() - t.start)
# 
# # adm details to add in to final dataframes
load("../../Data/NAM/Key_code_name_NAM.RData")

#### Uncomment/comment out one of the following chunks, depending on the distance kernel just used for fitting above gravity model
## 1. Namibia, estimate trips using power distance kernel ##
# NAM.basic.pwr <- mobility::summarize_mobility(out.NAM)
# NAM.predicted.basic.pwr <- estimate_trips(Key_code_name_NAM,
#                                           dist.fx = "pwr",
#                                           M,
#                                           N,
#                                           D,
#                                           theta=NAM.basic.pwr['theta', 'Mean'],
#                                           omega.1=NAM.basic.pwr['omega_1', 'Mean'],
#                                           omega.2=NAM.basic.pwr['omega_2', 'Mean'],
#                                           gam = NAM.basic.pwr[rownames(NAM.basic.pwr) %like% 'gamma', 'Mean'],
#                                           model = 'BasicModel-pwr',
#                                           district.ID = districts,
#                                           district.label = rank.Id,
#                                           location = "Namibia",
#                                           counts=TRUE)
# 
# save(NAM.predicted.basic.pwr, file='../../Data/NAM/GravityModelEstimates/NAM_adm2_daily_estimates_Basic_gm_dist_pwr.RData')

# # 2. Namibia, estimate trips using exponential distance kernel ##
# NAM.basic.exp <- mobility::summarize_mobility(out.NAM)
# NAM.predicted.basic.exp <- estimate_trips(Key_code_name_NAM,
#                                           dist.fx = "exp",
#                                           M,
#                                           N,
#                                           D,
#                                           theta=NAM.basic.exp['theta', 'Mean'],
#                                           omega.1=NAM.basic.exp['omega_1', 'Mean'],
#                                           omega.2=NAM.basic.exp['omega_2', 'Mean'],
#                                           gam = NAM.basic.exp[rownames(NAM.basic.exp) %like% 'gamma', 'Mean'],
#                                           model = 'BasicModel-exp',
#                                           district.ID = districts,
#                                           district.label = rank.Id,
#                                           location = "Namibia",
#                                           counts=TRUE)
# # save(NAM.predicted.basic.exp, file='../../Data/NAM/GravityModelEstimates/NAM_adm2_daily_estimates_Basic_gm_dist_exp.RData')

### Kenya adm2 ####
M <- load.obj(1, '../../Data/KEN/KEN.daily.avg.matrix.RData')
D <- load.obj(1, '../../Data/KEN/KEN_adm2_ppp_distance.RData')
N <- load.obj(1, '../../Data/KEN/KEN_adm2_ppp_pop2010.RData')
rownames(N) <- N$fid
N <- t(N$pop2010)
N <- as.vector(N)
names(N) <- colnames(D)

print(n_chain)
print(n_burn)
print(n_samp)
print(n_thin)
out.KEN <- fit_basic_gravity(M,
                             D,
                             N,
                             n_chain = n_chain,
                             n_burn = n_burn,
                             n_samp = n_samp,
                             n_thin = n_thin,
                             prior = NULL,
                             DIC = TRUE,
                             parallel = TRUE)
# save(out.KEN, file='../../Data/KEN/GravityModelEstimates/KEN_adm2_daily_estimates_Basic_gm_dist_pwr_chains.RData')
# print('Saved: ../../Data/KEN/GravityModelEstimates/KEN_adm2_daily_estimates_Basic_gm_dist_pwr_chains.RData')
print(Sys.time() - t.start)

# details to add in adm
load("../../Data/KEN/Key_code_name_KEN.RData")
#### Uncomment/comment out one of the following, depending on the distance kernel just used for fitting above gravity model ####

## 1. Kenya, estimate trips using power distance kernel ###
# KEN.basic.pwr <- mobility::summarize_mobility(out.KEN)
# KEN.predicted.basic.pwr <- estimate_trips(Key_code_name_KEN,
#                                           dist.fx = "pwr",
#                                           M,
#                                           N,
#                                           D,
#                                           theta=KEN.basic.pwr['theta', 'Mean'],
#                                           omega.1=KEN.basic.pwr['omega_1', 'Mean'],
#                                           omega.2=KEN.basic.pwr['omega_2', 'Mean'],
#                                           gam = KEN.basic.pwr[rownames(KEN.basic.pwr) %like% 'gamma', 'Mean'],
#                                           model = 'BasicModel-pwr',
#                                           district.ID = districts,
#                                           district.label = rank.Id,
#                                           location = "Kenya",
#                                           counts=TRUE)
# # save(KEN.predicted.basic.pwr, file='../../Data/KEN/GravityModelEstimates/KEN_adm2_daily_estimates_Basic_gm_dist_pwr.RData')

load("../../Data/KEN/GravityModelEstimates/KEN_adm2_daily_estimates_Basic_gm_dist_exp_chains.RData")

## 2. Kenya, estimate trips using exponential distance kernel ##
KEN.basic.exp <- mobility::summarize_mobility(out.KEN)
KEN.predicted.basic.exp <- estimate_trips(Key_code_name_KEN,
                                          dist.fx = "exp",
                                          M,
                                          N,
                                          D,
                                          theta=KEN.basic.exp['theta', 'Mean'],
                                          omega.1=KEN.basic.exp['omega_1', 'Mean'],
                                          omega.2=KEN.basic.exp['omega_2', 'Mean'],
                                          gam = KEN.basic.exp[rownames(KEN.basic.exp) %like% 'gamma', 'Mean'],
                                          model = 'BasicModel-exp',
                                          district.ID = districts,
                                          district.label = rank.Id,
                                          location = "Kenya",
                                          counts=TRUE)
# save(KEN.predicted.basic.exp, file='../../Data/KEN/GravityModelEstimates/KEN_adm2_daily_estimates_Basic_gm_dist_exp.RData')
