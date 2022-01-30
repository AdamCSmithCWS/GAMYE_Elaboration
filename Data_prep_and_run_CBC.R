
library(tidyverse)
library(bbsBayes)
library(cmdstanr)
setwd("C:/Users/adam_/OneDrivedelete/Documents/GitHub/GAMYE_Elaboration")


source("functions/neighbours_define_alt.R")
species = "Zonotrichia_albicollis"
species_f = "species"
data_1 = read.csv(paste0("data/",species,"_modeled_records.csv"))

data_1 <- data_1 %>% 
  mutate(strat = paste(country,state,bcr,sep = "-")) # making a strat value to match BBS strata



# strata_df ---------------------------------------------------------------

strata_map <- bbsBayes::load_map(stratify_by = "bbs_usgs") %>% 
  rename(strat = ST_12) 

strat_df = data_1 %>% 
  select(stratum,strata_vec,
         strat,
         stratum_area_km2,
         bcr,state,country,
         circles_per_stratum,nonzero_circles) %>% 
  distinct() %>% 
  group_by(stratum,strata_vec,
           strat,
           bcr,state,country,
           circles_per_stratum,nonzero_circles) %>% 
  summarise(stratum_area_km2 = sum(stratum_area_km2)) %>% 
    mutate(non_zero = nonzero_circles/circles_per_stratum)

nstrata = max(strat_df$strata_vec)


nonzeroweight <- as.integer(strat_df$non_zero)


# Neighbour relationships -------------------------------------------------

realized_strata_map <- strata_map %>% 
  right_join(.,strat_df,by = "strat")

neighbours = neighbours_define(realized_strata_map,
                  species = species,
                  alt_strat = "strata_vec",
                  plot_dir = "maps/",
                  plot_file = "_strata_map")


N_edges = neighbours$N_edges
node1 = neighbours$node1
node2 = neighbours$node2

# CBC circles -------------------------------------------------------------

sites_df <- data_1 %>% 
  select(circle,circle_vec,
         strat,strata_vec,
         lon,lat) %>% 
  distinct() %>% 
  arrange(circle_vec)
  
  
nsites = max(sites_df$circle_vec)

nsites_strata <- sites_df %>% 
  arrange(strata_vec,
          circle_vec) %>% 
  group_by(strata_vec) %>% 
  summarise(nsites = n())

nsites_strata <- as.integer(nsites_strata$nsites)
maxnsites_strata <- max(nsites_strata)

ste_mat <- matrix(data = 0,
                  nrow = nstrata,
                  ncol = maxnsites_strata)
for(i in 1:nstrata){
  ste_mat[i,1:nsites_strata[i]] <- sites_df[which(sites_df$strata_vec == i),"circle_vec"]
}



# GAM basis ---------------------------------------------------------------

source("Functions/GAM_basis_function_mgcv.R")
years_df <- data.frame(Year = min(data_1$year_vec):max(data_1$year_vec))

nknots = floor(nrow(years_df)/4)

GAM_year <- gam_basis(years_df$Year,
                      nknots = nknots,
                      sm_name = "Year")

nknots_year = GAM_year$nknots_Year
year_basis = GAM_year$Year_basis



# Data list ---------------------------------------------------------------


ncounts = nrow(data_1)
nyears = max(data_1$year_vec)

count = data_1$how_many
strat = data_1$strata_vec
year = data_1$year_vec
site = data_1$circle_vec
hours = data_1$scaled_effort




stan_data = list(#scalar indicators
  nsites = nsites,
  nstrata = nstrata,
  ncounts = ncounts,
  nyears = nyears,
  
  #basic data
  count = count,
  strat = strat,
  year = year,
  site = site,
  
  #spatial structure
  N_edges = N_edges,
  node1 = node1,
  node2 = node2,
  
  #GAM structure
  nknots_year = nknots_year,
  year_basis = year_basis,
  
  #Effort information
  hours = hours,
  
  #Ragged array information to link sites to strata
  nsites_strata = nsites_strata,
  maxnsites_strata = maxnsites_strata,
  ste_mat = ste_mat,
  #nu = nu,
  
  #weights
  nonzeroweight = nonzeroweight
)




# Fit model ---------------------------------------------------------------


output_dir <- "output/"
out_base <- paste0(species,"_CBC")
csv_files <- paste0(out_base,"-",1:3,".csv")


save(list = c("stan_data",
              "neighbours",
              "realized_strata_map",
              "data_1"),
     file = paste0("data/",species_f,"CBC_data.RData"))


print(paste("beginning",species,"with",nstrata,"strata",Sys.time()))

mod.file = "models/gamye_iCAR_CBC.stan"

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
                             sdbeta = runif(nknots_year,0.01,0.1),
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


#stanfit1 <- as_cmdstan_fit(files = paste0(output_dir,csv_files))


save(list = c("stanfit","stan_data","csv_files"),
     file = paste0(output_dir,"/",out_base,"_gamye_iCAR.RData"))






