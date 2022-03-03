


library(tidyverse)
library(bbsBayes)
library(cmdstanr)
library(patchwork)
setwd("C:/GitHub/GAMYE_Elaboration")

filter_low <- FALSE

source("functions/neighbours_define_alt.R")
#species = "Zonotrichia_albicollis"
species = "Cinclus_mexicanus"
# species = "Haliaeetus_leucocephalus"
# species = "Sitta_canadensis"


species_f = species


load(paste0("data/",species_f,"CBC_data.RData"))

stan_data[["N_edges"]] <- NULL
stan_data[["node1"]] <- NULL
stan_data[["node2"]] <- NULL




output_dir <- "output/"
out_base <- paste0(species_f,"_CBC_Non_spatial")
csv_files <- paste0(out_base,"-",1:3,".csv")



ncounts <- stan_data$ncounts
nstrata <- stan_data$nstrata
nknots_year <- stan_data$nknots_year
nsites <- stan_data$nsites
nyears <- stan_data$nyears

print(paste("beginning",out_base,Sys.time()))

mod.file = "models/gamye_CBC_Non_spatial.stan"

## compile model
model <- cmdstan_model(mod.file)


# Initial Values ----------------------------------------------------------


init_def <- function(){ list(noise_raw = rnorm(ncounts,0,0.1),
                             strata_raw = rnorm(nstrata,0,0.1),
                             STRATA = 0,
                             sdstrata = runif(1,0.01,0.1),
                             #eta = 0,
                             yeareffect_raw = matrix(rnorm(nstrata*nyears,0,0.1),nrow = nstrata,ncol = nyears),
                             ste_raw = rnorm(nsites,0,0.1),
                             sdnoise = runif(1,0.01,0.2),
                             sdb = runif(1,0.01,0.1),
                             sdp = runif(1,0.01,0.1),
                             b_raw = rnorm(nstrata,0,0.01),
                             p_raw = rnorm(nstrata,0,0.01),
                             B = 0,
                             P = 0,
                             sdste = runif(1,0.01,0.2),
                             sdbeta = runif(nstrata,0.01,0.1),
                             sdBETA = runif(1,0.01,0.1),
                             sdyear = runif(nstrata,0.01,0.1),
                             BETA_raw = rnorm(nknots_year,0,0.1),
                             beta_raw = matrix(rnorm(nknots_year*nstrata,0,0.01),nrow = nstrata,ncol = nknots_year))}

stanfit <- model$sample(
  data=stan_data,
  refresh=200,
  chains=3, iter_sampling=1000,
  iter_warmup=1000,
  parallel_chains = 3,
  #pars = parms,
  adapt_delta = 0.8,
  max_treedepth = 14,
  seed = 123,
  init = init_def,
  output_dir = output_dir,
  output_basename = out_base)


#stanfit <- as_cmdstan_fit(files = paste0(output_dir,csv_files))


save(list = c("stanfit","stan_data","csv_files"),
     file = paste0(output_dir,"/",out_base,"_gamye_iCAR.RData"))


