# SPECIES MCMC data-prep -------------------------------------------------------
library(tidyverse)
library(cmdstanr)
library(sf)
setwd("C:/GitHub/GAMYE_Elaboration")


source("functions/neighbours_define_alt.R")
species = "Red Knot"
species_f = "Red_Knot"
sp = species


load("Data/shorebird_full_observation_dataset.Rdata")
load("Data/shorebird_site_map.RData") 
load("data/shorebird_hexagon_grid.RData")



  load(paste0("data/",species_f,"shorebird_data.RData"))
  
  mod.file = "models/GAMYE_iCAR_split_shorebird_two_season.stan"
 
  stan_data[["nknots_year_m1"]] <- stan_data$nknots_year-1
  
  stan_data$year_basis[,stan_data$nknots_year] <- c(1:stan_data$nyears)-mean(1:stan_data$nyears)

  
  
  init_def <- function(){ list(noise_raw = rnorm(stan_data$ncounts,0,0.1),
                               ste_raw = rnorm(stan_data$nsites,0,0.1),
                               STRATA = 0,
                               yeareffect_raw = rnorm(stan_data$nyears,0,0.1),
                               beta_raw_season_1 = rnorm(stan_data$ndays,0,0.1),
                               beta_raw_season_2 = rnorm(stan_data$ndays,0,0.1),
                               sdnoise = 0.2,
                               sdste = 0.1,
                               sdBETA = 1,
                               sdbeta = runif(stan_data$nknots_year_m1,0.01,0.02),
                               sdseason = c(0.1,0.1),
                               sdyear = 0.1,
                               BETA_raw = rnorm(stan_data$nknots_year_m1,0,0.1),
                               beta_raw = matrix(rnorm((stan_data$nknots_year_m1)*stan_data$nstrata,0,0.01),
                                                 nrow = stan_data$nstrata,ncol = stan_data$nknots_year_m1),
                               BETA_raw_lin = rnorm(1,0,0.01),
                               beta_raw_lin = rnorm(stan_data$nstrata,0,0.01),
                               sdbeta_lin = runif(1,0.01,0.07))}
  
  
  

  # Fit model ---------------------------------------------------------------
  
  
  output_dir <- "output/"
  out_base <- paste0(species_f,"_Shorebird_iCAR_split")
  csv_files <- paste0(out_base,"-",1:3,".csv")

  
  
  print(paste("beginning",out_base,Sys.time()))
  
  
  ## compile model
  model <- cmdstan_model(mod.file)
  
  
  stanfit <- model$sample(
    data=stan_data,
    refresh=200,
    chains=3, iter_sampling=2000,
    iter_warmup=2000,
    parallel_chains = 3,
    #pars = parms,
    adapt_delta = 0.95,
    max_treedepth = 14,
    seed = 123,
    init = init_def,
    output_dir = output_dir,
    output_basename = out_base)
  
  
  #stanfit1 <- as_cmdstan_fit(files = paste0(output_dir,csv_files))
  
  
  save(list = c("stanfit","stan_data","csv_files"),
       file = paste0(output_dir,"/",out_base,"_gamye_iCAR.RData"))
  
  


