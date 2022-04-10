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
seasons <- read.csv("data/Shorebird_seasons.csv")
regs_alt_season <- seasons[which(seasons$Species == species),"regions_w_alt_season"]
regs_alt_season <- strsplit(regs_alt_season,
                            split = " - ",
                            fixed = TRUE)[[1]]


 source("functions/GAM_basis_function_mgcv.R")

FYYYY = 1980


dts <- filter(ssData,CommonName == sp,
              YearCollected >= FYYYY)
nyrs_study <- 2019-FYYYY #length of the time-series being modeled


# selecting which sites have sufficient data to include -------------------

dts$present <- FALSE

dts[which(dts$ObservationCount > 0),"present"] <- TRUE

#number of non-zero observations by site
nobs_site <- dts %>% 
  filter(present == TRUE) %>% 
  group_by(Region,SurveyAreaIdentifier) %>% 
  summarise(nobs = n())

yspn = function(x){
  diff(range(x))
}
#number of years with non-zero observations, by site
nyrs_site <- dts %>% 
  filter(present == TRUE) %>%  
  group_by(Region,SurveyAreaIdentifier,YearCollected) %>% 
  summarise(nobs = n()) %>% 
  group_by(Region,SurveyAreaIdentifier) %>% 
  summarise(nyears = n(),
            span_years = yspn(YearCollected),
            fyear = min(YearCollected),
            lyear = max(YearCollected))


min_nyears <- 2 # species must have been observed in at least 2 years
minspan <- 10 # species must have been observed at least twice, 10-years apart
#sites with > 10 year span of non-zero, observations
sites_keep <- nyrs_site[which(nyrs_site$span_years >= minspan,
                              nyrs_site$nyears >= min_nyears),"SurveyAreaIdentifier"]

# drop sites that don't meet above criterion ------------------------------------
dts <- filter(dts,SurveyAreaIdentifier %in% sites_keep$SurveyAreaIdentifier) 


# selecting strata to retain ----------------------------------------------
#number of years with non-zero observations, by hexagonal strata (hex_name)
nyrs_region <- dts %>% 
  filter(present == TRUE) %>%  
  group_by(hex_name,YearCollected) %>% 
  summarise(nobs = n()) %>% 
  group_by(hex_name) %>% 
  summarise(nyears = n(),
            span_years = yspn(YearCollected),
            fyear = min(YearCollected),
            lyear = max(YearCollected))

#strats with 20 or more year span of non-zero, observations - 
p_time_series = 0.5 #strata are required to have data that span 50% of the time-series in a region.
regions_keep <- nyrs_region[which(nyrs_region$span_years >= nyrs_study*p_time_series),"hex_name"]



# drop strata that don't meet above criterion -----------------------------
dts <- filter(dts,hex_name %in% regions_keep$hex_name) 

real_grid <- poly_grid %>% 
  filter(.,hex_name %in% unlist(regions_keep$hex_name)) 


# Grouping the strata by the broad regions where most sites fall ----------
# these broad regions are somewhat subjective, but are mostly helpful
# for defining the two separate seasonal patterns of migration for 
# Red Knot (and many other species modeled in Smith et al. 2022)
dom_value <- function(x){
  tt = table(x)
  tt = sort(tt)
  return(names(tt)[1])
}
hex_by_reg <- dts %>% group_by(hex_name) %>% 
  summarise(Region = dom_value(Region))

real_grid_regs <- left_join(real_grid, hex_by_reg,by = "hex_name")

dts <- dts %>% select(-Region) %>% 
  left_join(.,hex_by_reg)

strats_dts <- data.frame(hex_name = real_grid$hex_name,
                         stratn = 1:length(real_grid$hex_name))

dts <- left_join(dts,strats_dts,by = "hex_name")
#hexagon grid map
real_grid <- inner_join(real_grid,strats_dts)
#hesagon grid map with the broad regions added
real_grid_regs <- inner_join(real_grid_regs,strats_dts)


# DOY season definition ---------------------------------------------------

# doy is an ordinal date of the year
fday = min(dts$doy)-1


dts <- dts %>% mutate(count = as.integer(ObservationCount),
                      year = as.integer(YearCollected),
                      yr = as.integer(year-(FYYYY-1)),
                      strat = stratn,
                      date = doy-fday,
                      site = as.integer(factor(SurveyAreaIdentifier))) 

nsites = max(dts$site)
nstrata = max(dts$strat)


## indexing of sites by strata for annual index calculations
sByReg = unique(dts[,c("site","stratn")])
sByReg <- arrange(sByReg,stratn,site)
# 
nsites_strat <- table(sByReg$stratn)
maxsites_strata <- max(nsites_strat)
ste_mat <- matrix(data = 0,
                  ncol = maxsites_strata,
                  nrow = nstrata)
for(j in 1:nstrata){
  ste_mat[j,1:nsites_strat[j]] <- as.integer(unlist(sByReg[which(sByReg$stratn == j),"site"]))
}

# Identifying the strata and region combinations --------------------------
# necessary for two separate seasons
strat_regions <- unique(dts[,c("Region","stratn")])
strat_regions <- strat_regions[order(strat_regions$stratn),]

strat_regions$seas_strat <- 1
strat_regions$seas_strat[which(strat_regions$Region %in% regs_alt_season)] <- 2


dts <- left_join(dts,strat_regions,by = c("Region","stratn"))

seasons <- as.matrix(strat_regions[,c("stratn","seas_strat")])

#hesagon grid map with the broad regions added
real_grid_regs_season <- inner_join(real_grid_regs,strat_regions)


# generate neighbourhoods -------------------------------------------------


neighbours = neighbours_define(real_grid_regs_season,
                               species = species,
                               alt_strat = "stratn",
                               plot_dir = "maps/",
                               plot_file = "_strata_map",
                               voronoi = TRUE,
                               strat_link_fill = 1000000)#extend the voronoi polygons 1000 km from centres to account for large spacing of some strata



N_edges = neighbours$N_edges
node1 = neighbours$node1
node2 = neighbours$node2






ncounts = nrow(dts)

nyears = max(dts$yr)

# GAM seasonal basis function ---------------------------------------------

nknots_seasont = 10
season_df <- data.frame(date = min(dts$date):max(dts$date))

GAM_season <- gam_basis(season_df$date,
                          nknots = nknots_seasont,
                          sm_name = "season")

ndays <- max(season_df$date)
season_basis <- GAM_season$season_basis
nknots_season <- GAM_season$nknots_season # 

# GAM year basis function ---------------------------------------------
years_df <- data.frame(Year = min(dts$yr):max(dts$yr))

nknots = floor(nrow(years_df)/4)

GAM_year <- gam_basis(years_df$Year,
                      nknots = nknots,
                      sm_name = "Year")

nknots_year = GAM_year$nknots_Year
year_basis = GAM_year$Year_basis

# 

  
  stan_data <- list(
    
    nyears = nyears,
    nstrata = nstrata,
    nsites = nsites,
    ncounts = ncounts,
    ndays = ndays,
    
    count = as.integer(unlist(dts$count)),
    year = as.integer(unlist(dts$yr)),
    site = as.integer(unlist(dts$site)),
    strat = as.integer(unlist(dts$stratn)),
    N_edges = neighbours$N_edges,
    node1 = neighbours$node1,
    node2 = neighbours$node2,
    nsites_strata = as.integer(nsites_strat),
    maxnsites_strata = maxsites_strata,
    ste_mat = ste_mat,
    nknots_year = nknots_year,
    year_basis = year_basis,
    nknots_season = nknots_season,
    season_basis = season_basis,
    date = as.integer(unlist(dts$date)),
    seas_strat = as.integer(unlist(dts$seas_strat)),
    seasons = seasons)
  
  mod.file = "models/GAMYE_iCAR_shorebird_two_season.stan"
 

  init_def <- function(){ list(noise_raw = rnorm(stan_data$ncounts,0,0.1),
                               ste_raw = rnorm(stan_data$nsites,0,0.1),
                               STRATA = 0,
                               yeareffect_raw = rnorm(stan_data$nyears,0,0.1),
                               beta_raw_season_1 = rnorm(stan_data$ndays,0,0.1),
                               beta_raw_season_2 = rnorm(stan_data$ndays,0,0.1),
                               sdnoise = 0.2,
                               sdste = 0.1,
                               sdBETA = 1,
                               sdbeta = runif(stan_data$nknots_year,0.1,0.2),
                               sdseason = c(0.1,0.1),
                               sdyear = 0.1,
                               BETA_raw = rnorm(stan_data$nknots_year,0,0.1),
                               beta_raw = matrix(rnorm(stan_data$nknots_year*stan_data$nstrata,0,0.01),
                                              nrow = stan_data$nstrata,ncol = stan_data$nknots_year))}
  
  
  

  # Fit model ---------------------------------------------------------------
  
  
  output_dir <- "output/"
  out_base <- paste0(species_f,"_Shorebird")
  csv_files <- paste0(out_base,"-",1:3,".csv")
  
  data_1 = dts
  realized_strata_map = real_grid_regs_season 
  save(list = c("stan_data",
                "neighbours",
                "realized_strata_map",
                "data_1"),
       file = paste0("data/",species_f,"shorebird_data.RData"))
  
  
  
  print(paste("beginning",species,"with",nstrata,"strata",Sys.time()))
  
  
  ## compile model
  model <- cmdstan_model(mod.file)
  
  
  stanfit <- model$sample(
    data=stan_data,
    refresh=200,
    chains=3, iter_sampling=4000,
    iter_warmup=2000,
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
  
  


