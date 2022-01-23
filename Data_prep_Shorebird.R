# SPECIES MCMC data-prep -------------------------------------------------------
library(tidyverse)
library(cmdstanr)
setwd("C:/Users/adam_/OneDrivedelete/Documents/GitHub/GAMYE_Elaboration")


source("functions/neighbours_define_alt.R")
species = "Red Knot"
species_f = "Red_Knot"
sp = species


load("data/shorebird_hexagon_grid.RData")
load("Data/shorebird_full_observation_dataset.Rdata")
load("Data/shorebird_site_map.RData") 
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

real_grid <- poly_grid %>% filter(hex_name %in% regions_keep$hex_name) 


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
sByReg = unique(dts[,c("site","strat")])
sByReg <- arrange(sByReg,strat,site)
# 
nsites_strat <- table(sByReg$strat)
maxsites_strata <- max(nsites_strat)
ste_mat <- matrix(data = 0,
                  ncol = maxsites_strata,
                  nrow = nstrata)
for(j in 1:nstrata){
  ste_mat[j,1:nsites_strat[j]] <- as.integer(unlist(sByReg[which(sByReg$strat == j),"site"]))
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


neighbours = neighbours_define(real_grid_regs,
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

nknots_season = 10
season_df <- data.frame(date = min(dts$date):max(dts$date))

GAM_season <- gam_basis(season_df$date,
                          nknots = nknots_season,
                          sm_name = "season")

ndays <- max(season_df$date)
season_basis <- GAM_season$season_basis
# 

# GAM year basis function ---------------------------------------------
years_df <- data.frame(Year = min(dts$yr):max(dts$yr))

nknots = floor(nrow(years_df)/4)

GAM_year <- gam_basis(years_df$Year,
                      nknots = nknots,
                      sm_name = "Year")

nknots_year = GAM_year$nknots_Year
year_basis = GAM_year$Year_basis

# 



  
  stan_data <- list(count = as.integer(unlist(dts$count)),
                    year = as.integer(unlist(dts$yr-midyear)),
                    year_raw = as.integer(unlist(dts$yr)),
                    site = as.integer(unlist(dts$site)),
                    strat = as.integer(unlist(dts$strat)),
                    date = as.integer(unlist(dts$date)),
                    seas_strat = as.integer(unlist(dts$seas_strat)),
                    
                    nyears = nyears,
                    nstrata = nstrata,
                    nsites = nsites,
                    ncounts = ncounts,
                    ndays = ndays,
                    
                    nsites_strat = as.integer(nsites_strat),
                    maxsites_strata = maxsites_strata,
                    sites = sites,
                    seasons = seasons,
                    
                    #site_size = sizes_by_site$size_cent,
                    
                    # season_basis = basis_season$season_basis,
                    season_basispred = basis_season$season_basispred,
                    nknots_season = basis_season$nknots_season,
                    
                    year_basispred = basis_year$year_basispred,
                    nknots_year = basis_year$nknots_year,
                    
                    #midyear = midyear,
                    
                    N_edges = car_stan_dat$N_edges,
                    node1 = car_stan_dat$node1,
                    node2 = car_stan_dat$node2)
  
  mod.file1 = "models/GAMYE_strata_two_season_gammaprior_beta_normal.stan"
  prior = "gamma"
  noise_dist1 = "normal"
  
  mod.file2 = "models/GAMYE_strata_two_season_gammaprior_beta.stan"
  prior = "gamma"
  noise_dist2 = "t"
  

  init_def <- function(){ list(noise_raw = rnorm(stan_data$ncounts,0,0.1),
                               alpha_raw = rnorm(stan_data$nsites,0,0.1),
                               ALPHA1 = 0,
                               year_effect_raw = rnorm(stan_data$nyears,0,0.1),
                               B_season_raw1 = rnorm(stan_data$ndays,0,0.1),
                               B_season_raw2 = rnorm(stan_data$ndays,0,0.1),
                               sdnoise = 0.2,
                               sdalpha = 0.1,
                               sdyear_gam = 1,
                               sdyear_gam_strat = runif(stan_data$nknots_year,0.1,0.2),
                               sdseason = c(0.1,0.1),
                               sdyear = 0.1,
                               B_raw = rnorm(stan_data$nknots_year,0,0.1),
                               b_raw = matrix(rnorm(stan_data$nknots_year*stan_data$nstrata,0,0.01),
                                              nrow = stan_data$nstrata,ncol = stan_data$nknots_year))}
  
  
  

  parms = c("sdnoise",
            #"nu", #
            "sdalpha",
            "b",
            "B",
            "alpha",
            "ALPHA1",
            "sdyear",
            "sdyear_gam_strat",
            "sdyear_gam",
            "year_effect",
            "sdseason",
            "B_season_raw1",
            "B_season_raw2",
            "season_pred",
            "n",
            "nsmooth",
            "N",
            "NSmooth",
            "log_lik")
  
}else{
  stan_data <- list(count = as.integer(unlist(dts$count)),
                    year = as.integer(unlist(dts$yr-midyear)),
                    year_raw = as.integer(unlist(dts$yr)),
                    site = as.integer(unlist(dts$site)),
                    strat = as.integer(unlist(dts$strat)),
                    date = as.integer(unlist(dts$date)),
                    
                    nyears = nyears,
                    nstrata = nstrata,
                    nsites = nsites,
                    ncounts = ncounts,
                    ndays = ndays,
                    
                    nsites_strat = as.integer(nsites_strat),
                    maxsites_strata = maxsites_strata,
                    sites = sites,
                    
                    #site_size = sizes_by_site$size_cent,
                    
                    # season_basis = basis_season$season_basis,
                    season_basispred = basis_season$season_basispred,
                    nknots_season = basis_season$nknots_season,
                    
                    year_basispred = basis_year$year_basispred,
                    nknots_year = basis_year$nknots_year,
                    
                    #midyear = midyear,
                    
                    N_edges = car_stan_dat$N_edges,
                    node1 = car_stan_dat$node1,
                    node2 = car_stan_dat$node2)
  
  mod.file1 = "models/GAMYE_strata_gammaprior_beta_normal.stan"
  prior = "gamma"
  noise_dist1 = "normal"
  
  mod.file2 = "models/GAMYE_strata_gammaprior_beta.stan"
  prior = "gamma"
  noise_dist2 = "t"
  
  init_def <- function(){ list(noise_raw = rnorm(stan_data$ncounts,0,0.1),
                               alpha_raw = rnorm(stan_data$nsites,0,0.1),
                               ALPHA1 = 0,
                               year_effect_raw = rnorm(stan_data$nyears,0,0.1),
                               B_season_raw = rnorm(stan_data$ndays,0,0.1),
                               sdnoise = 0.2,
                               sdalpha = 0.1,
                               sdyear_gam = 1,
                               sdyear_gam_strat = runif(stan_data$nknots_year,0.1,0.2),
                               sdseason = 0.1,
                               sdyear = 0.1,
                               B_raw = rnorm(stan_data$nknots_year,0,0.1),
                               b_raw = matrix(rnorm(stan_data$nknots_year*stan_data$nstrata,0,0.01),
                                              nrow = stan_data$nstrata,ncol = stan_data$nknots_year))}
  
  
  
  parms = c("sdnoise",
            #"nu", #
            "sdalpha",
            "b",
            "B",
            "alpha",
            "ALPHA1",
            "sdyear",
            "sdyear_gam_strat",
            "sdyear_gam",
            "year_effect",
            "sdseason",
            "B_season_raw",
            "season_pred",
            "n",
            "nsmooth",
            "N",
            "NSmooth",
            "log_lik")
}

# Explore site-level trajectories of observed means ------------------------

# site_means <- dts %>% group_by(year,SurveyAreaIdentifier) %>%
#   summarise(means = log(mean(count,na.rm = T)+1))
# nrow(site_means)/nyears
# smp = ggplot(data = site_means,aes(x = year,y = means,colour = SurveyAreaIdentifier))+
#   geom_point(alpha = 0.05)+
#   geom_smooth(method = "lm",se = FALSE)+
#   #scale_y_continuous(trans = "log10")+
#   theme(legend.position = "none")
# print(smp)#   

# prepare stan data -------------------------------------------------------


save(list = c("stan_data",
              "dts",
              "real_grid",
              "real_grid_regs",
              "strats_dts",
              "strat_regions",
              "mod.file1",
              "prior",
              "noise_dist1",
              "mod.file2",
              "noise_dist2",
              "parms",
              "init_def",
              "vintj",
              "nb_db",
              "cc"),
     file = paste0("data/data",sp,"_cmdstanr_data",p_time_series,"_",minspan,"_",min_nyears,".RData"))




} ### end species loop


 
 
 
 #print graphs
 
 
 pdf(paste0("Figures/","All_seasonal_counts ",grid_spacing/1000,".pdf"),
     width = 11,height = 8.5)
 for(sp in sps){
   print(mean_counbts_doy_out[[sp]]+
           labs(title = sp))
   
 }
 dev.off()
 
 pdf(file = paste0("Figures/","ALL_Strata_",grid_spacing/1000,".pdf"))
 
 for(sp in sps){
   print(ggp_out[[sp]]+
           labs(title = sp))
   
 }
 
 dev.off()
 
 pdf(paste0("Figures/","All_annual_counts ",grid_spacing/1000,".pdf"),
     width = 11,height = 8.5)
 for(sp in sps){
   print(mean_counbts_year_out[[sp]]+
           labs(title = sp))
   
 }
 dev.off()
 

