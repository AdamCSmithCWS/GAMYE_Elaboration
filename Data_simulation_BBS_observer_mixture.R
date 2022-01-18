### simulating fake BBS time-series data

library(tidyverse)
library(bbsBayes)
library(mgcv) #has functions to simulate data from a GAM

# select real BBS data for CEWA -------------------------------------------

BBS_data <- stratify("bbs_usgs")


species = "Pine Warbler"  
species = "Yellow-headed Blackbird"  
species_f <- gsub(species,pattern = " ",replacement = "_")
  
  
  for(tp in c("non_linear","linear")){
source("Functions/prepare-jags-data-alt.R")

real_data <- prepare_jags_data_alternate(strat_data = BBS_data,
                          species_to_run = species,
                          model = "gamye",
                          min_year = NULL) #n-year time-series ???
min_year <- min(real_data$r_year)

to_save <- c("real_data","min_year")

# Generate balanced dataset -----------------------------------------------

# create full data frame 
real_df <- data.frame(Year = real_data$r_year,
                      Year_Factored = real_data$year,
                      Stratum = real_data$strat_name,
                      Stratum_Factored = real_data$strat,
                      Observer = real_data$ObsN,
                      Route = real_data$route,
                      First_Year = real_data$firstyr) %>% 
  mutate(Route_Factored = as.integer(factor(Route)))

to_save <- c(to_save,"real_df")


# dataframe of routes and strata
routes_df <- real_df %>% 
  select(Route,Route_Factored,Stratum,Stratum_Factored) %>% 
  distinct() 
to_save <- c(to_save,"routes_df")

# dataframe of first year for each observer on a route
observer_route_df <- real_df %>% 
  select(Route,Route_Factored,Observer,Year,First_Year) %>%
  filter(First_Year == 1) %>% 
  distinct() %>% 
  arrange(Route,Year)
to_save <- c(to_save,"observer_route_df")

# Generate balanced dataset with each year added in the realised first year it was surveyed
balanced <- NULL
for(rr in unique(observer_route_df$Route)){
  tmp <- observer_route_df %>% filter(Route == rr) %>% 
    arrange(Year)
  if(nrow(tmp) > 1){
  for(i in 1:(nrow(tmp)-1)){
    tmp1 <- data.frame(Route = rr,
                       Route_Factored = tmp[i,"Route_Factored"],
                       Observer = tmp[i,"Observer"],
                       Year = tmp[i,"Year"]:(tmp[i+1,"Year"]-1))
    balanced <- bind_rows(balanced,tmp1)
  }
  }else{
    i = 0
  }
  tmp1 <- data.frame(Route = rr,
                     Route_Factored = tmp[i+1,"Route_Factored"],
                     Observer = tmp[i+1,"Observer"],
                     Year = tmp[i+1,"Year"]:2019)
  balanced <- bind_rows(balanced,tmp1)
  
  
}


# Fill in strata info on balanced
balanced <- balanced %>% 
  left_join(.,routes_df,by = c("Route","Route_Factored"))


# Generate mean smooth ----------------------------------------

source("Functions/GAM_basis_function_mgcv.R")
years_df <- data.frame(Year = min(balanced$Year):max(balanced$Year))

nknots = floor(nrow(years_df)/4)

GAM_year <- gam_basis(years_df$Year,
                      nknots = nknots,
                      sm_name = "Year")

to_save <- c(to_save,"GAM_year")

#select random initial BETA values to seed the spatial variation in trajetories
set.seed(2017)
BETA_mid <- rnorm(GAM_year$nknots_Year,0,1.5)

# if linear pattern is desired, set all but final BETA == 0
if(tp == "linear"){ #for linear trend sets all BETAs except 13 to 0
BETA_mid[1:(GAM_year$nknots_Year-1)] <- 0
}
to_save <- c(to_save,"BETA_mid")

mean_log_smooth <-  GAM_year$Year_basis %*% BETA_mid

#plot(exp(mean_log_smooth),type = "l",ylim = c(0,max(exp(mean_log_smooth))))

## strata neighbourhoods Generate ---------------------------------
nstrata <- real_data$nstrata


strata_df <- balanced %>% 
  select(Stratum,Stratum_Factored) %>% 
  distinct() %>% 
  arrange(Stratum_Factored)

to_save <- c(to_save,"strata_df")


strata_map <- bbsBayes::load_map(stratify_by = "bbs_usgs") %>% 
  rename(Stratum = ST_12) %>% 
  right_join(.,strata_df,by = "Stratum") %>% 
  arrange(Stratum_Factored) ### this arranging is critical to the correct neighbourhoods

to_save <- c(to_save,"strata_map")

st_coord <- suppressWarnings(sf::st_centroid(strata_map)) %>% 
  sf::st_coordinates()%>% 
  as.data.frame() %>% 
  mutate(Stratum_Factored = 1:nstrata)# this 1:nstrata re-definition works because of the 
                    ### arrange statement in the strata-map block 10-lines up.



strata_df <- strata_df %>% 
  left_join(.,st_coord,by = "Stratum_Factored")

source("Functions/neighbours_define.R")

neighbours <- neighbours_define(real_strata_map = strata_map,
                                plot_dir = "maps/",
                                species = paste0("Simulated",species_f),
                                alt_strat = "Stratum")
to_save <- c(to_save,"neighbours")

neighbours2 <- neighbours_define(real_strata_map = strata_map,
                                 plot_dir = "maps/",
                                 species = paste0("Simulated_voronoi",species_f),
                                 alt_strat = "Stratum",
                                 voronoi = TRUE)
to_save <- c(to_save,"neighbours2")

## Generate stratum smooths and intercepts ----------------------------------------

### use simple, smooth, x-y coordinate variation in betas and stratas

strata_df <- strata_df %>% 
  mutate(yscale = scale(Y,scale = TRUE),
         xscale = scale(X,scale = TRUE),
         sumxy = yscale+xscale)


# strat_tempplot <- ggplot(data = strata_df,aes(x = X,y = Y))+
#   geom_point(aes(colour = sumxy),size = 3)+
#   scale_colour_viridis_c()
# 
# print(strat_tempplot)




## correlation matrix
neigh_mat <- neighbours$adj_matrix

### strata intercepts
nstrata <- nrow(strata_df)
STRATA_True <- 3 #peak abundance ~ 20 birds/route
strata_df <- strata_df %>% 
  mutate(strata_True = STRATA_True - abs(yscale)) # abundance peaks at middle latitudes
strata_True <- as.numeric(strata_df$strata_True)


### strata betas
Beta_True <- matrix(NA,nrow = nknots,
                    ncol = nstrata)
sum_xy <- as.numeric(strata_df$sumxy)
yscale <- as.numeric(strata_df$yscale)

for(k in 1:nknots){
  
Beta_True[k,] <- yscale*0.75 + BETA_mid[k]
if(tp == "linear" & k < nknots){ #for linear trend sets all BETAs except 13 to 0
  Beta_True[k,] <- rep(0,nstrata)
}

  }

  to_save <- c(to_save,"Beta_True")
  to_save <- c(to_save,"strata_True")
  
  BETA_True <- rowMeans(Beta_True)
  to_save <- c(to_save,"BETA_True")
  
 
## stratum-level smooths ----------------------------------------

  strat_log_smooths <- GAM_year$Year_basis %*% Beta_True 

  true_log_smooths <- as.data.frame(strat_log_smooths)
  true_log_smooths[,"Year"] <- min_year : 2019
  
  log_smooth_plot <- true_log_smooths %>% 
    pivot_longer(cols = starts_with("V"),
                 names_to = "Stratum_Factored",
                 names_prefix = "V",
                 values_to = "Smooth")%>% 
    mutate(index = exp(Smooth),
           Stratum_Factored = as.integer(Stratum_Factored)) %>% 
    left_join(.,strata_df,by = c("Stratum_Factored")) %>% 
    arrange(Y,Year) 

  pfs <- ggplot(data = log_smooth_plot,aes(x = Year,y = index,colour = Y))+
    geom_line(size = 2)+
    scale_color_viridis_c()+
    scale_y_continuous(limits = c(0,NA))+
    facet_wrap(~Stratum,scales = "free_y",
               nrow = ceiling(sqrt(nstrata)),
               ncol = ceiling(sqrt(nstrata)))
  
  #print(pfs)
  
  ## Add random annual fluctuations ----------------------------------------

  sdyeareffect <- 0.1 # ~10% mean annual fluctuation
  
  ye_funct <- function(x,sd = sdyeareffect){
    ye = rnorm(length(x),0,sd)
  }
  
  
  log_true_traj <- log_smooth_plot %>% 
    group_by(Stratum) %>% 
    mutate(YearEffect = ye_funct(Smooth),
           True_log_traj = Smooth + YearEffect,
           True_traj = exp(True_log_traj))
  

  pf <- ggplot(data = log_true_traj,aes(x = Year,y = True_traj,colour = Y))+
    geom_line(size = 2)+
    scale_color_viridis_c()+
    scale_y_continuous(limits = c(0,NA))+
    facet_wrap(~Stratum,scales = "free_y",
               nrow = ceiling(sqrt(nstrata)),
               ncol = ceiling(sqrt(nstrata)))
  
  
  pdf(paste0("Figures/",species_f,"_",tp,"_True_smooth.pdf"),
      width = 11,
      height = 8.5)
  print(pfs)
  print(pf)
  
 dev.off()
 

# Add observer, route intercepts ----------------------------

  balanced <- balanced %>% left_join(.,log_true_traj,
                                     by = c("Stratum","Year","Stratum_Factored")) %>% 
    select(-c(index,True_traj))
  
  ## OBserver Effects
  sdobs <- 0.2
  psi <- 0.6 #proportion of the observers that usually count the species
  psi2 = 1-psi #proportion of the observers that rarely count the species
  mu_obs = -4 #log-scale adjustment to the mean observer effect for observers that rarely count the species
  
  nobservers <- length(unique(balanced$Observer))
  True_observer_effects <- rnorm(nobservers,0,sdobs) #True observer effects by 
  rare_obs <- sample(1:nobservers,size = nobservers*(psi2))
  True_observer_effects[rare_obs] <- True_observer_effects[rare_obs]+mu_obs  
  observer_df <- data.frame(True_observer_effects = True_observer_effects,
                            Observer = unique(balanced$Observer)) %>% 
    mutate(Observer_Factored = as.integer(factor(Observer)))

  to_save <- c(to_save,"observer_df","psi","psi2",
               "mu_obs")
  
  balanced <- balanced %>% 
    left_join(observer_df,by = "Observer")
  
  ## Route Effects
  sdroute <- 0.2
  nroutes <- max(routes_df$Route_Factored)
  
  routes_df <- routes_df %>% 
    mutate(True_route_effects = rnorm(nroutes,0,sdroute))
  
  to_save <- c(to_save,"routes_df")
  
  balanced <- balanced %>% 
    left_join(.,routes_df,
              by = c("Route",
                     "Route_Factored",
                     "Stratum",
                     "Stratum_Factored"))
    
  
  ## Stratum Effects
  
  strata_df <- strata_df %>% 
    mutate(True_strata_effects = strata_True)
  
  to_save <- c(to_save,"strata_df")
  
  balanced <- balanced %>% 
    left_join(.,strata_df,
              by = c("X","Y",
                     "Stratum",
                     "Stratum_Factored",
                     "yscale",
                     "xscale",
                     "sumxy",
                     "strata_True"))
  
  
  
  

# Add simulated counts ----------------------------------------------------

  sdnoise = 0.1
  over_disp <- rnorm(nrow(balanced),0,sdnoise)
balanced <- balanced %>% 
    mutate(log_expected = over_disp + True_log_traj + True_observer_effects + True_route_effects + True_strata_effects,
           expected_count = exp(log_expected),
           count = rpois(nrow(balanced),expected_count))%>% 
  mutate(Year_Index = Year-(min(Year)-1))

  
to_save <- c(to_save,"balanced")


realised <- real_df %>% select(-c("Observer")) %>% 
  left_join(.,balanced,
            by = c("Year",
                   "Stratum",
                   "Stratum_Factored",
                   "Route",
                   "Route_Factored")) %>% 
  mutate(Year_Index = Year-(min(Year)-1))


to_save <- c(to_save,"realised")

if(tp == "non_linear"){
  # original observation data with updated observer adn route info ----------
original_data_df <- data.frame(Year = real_data$r_year,
                        Stratum = real_data$strat_name,
                        Stratum_Factored = real_data$strat,
                        Observer = real_data$ObsN,
                        Route = real_data$route,
                        First_Year = real_data$firstyr,
                        count = real_data$count) %>% 
  mutate(Route_Factored = as.integer(factor(Route)),
         Observer_Factored = as.integer(real_data$ObsN))

original_data_df <- original_data_df %>% 
  mutate(Year_Index = Year-(min(Year)-1))

# save(list = c(to_save,"original_data_df"),
#      file = paste0("Data/Real_data_",species_f,"_BBS.RData"))
save(list = c(to_save,"original_data_df"),
     file = paste0("Data/Real_data_",species_f,"_obs_mix_BBS.RData"))

}

# remove all but two years from each observer and route in some st --------

routes_mask <- routes_df %>% 
  filter(grepl(Stratum,pattern = "CA-",fixed = TRUE) )

# minimum number of observation events to retain so that all routes and observers are included

event_mask_retain <- realised %>% 
  filter(Route %in% routes_mask$Route) %>% 
  group_by(Route,Observer,First_Year) %>% 
  sample_n(size = 1) %>% 
  as.data.frame()

realised_mask <- realised %>% 
  filter(!Route %in% routes_mask$Route) %>% 
  bind_rows(.,event_mask_retain)
  

to_save <- c(to_save,"realised_mask","routes_mask","event_mask_retain")



### realised is a final dataset to demonstrate a realistic BBS dataset
### with a known pattern of population change

### realised_mask is an equivalent dataset that has data mostly missing
### from one portion of the range 
### mostly missing means only 1-years observations are available from each route ad observer
### so that the route and observer structure is retained, but there is no
### information available on change over time, after the observer and route intercepts



save(list = to_save,
     file = paste0("Data/Simulated_data_",species_f,"_",tp,"_obs_mix_BBS.RData"))

}

  slp = function(x,y){
    m = lm(x~y)
    sl = coefficients(m)[[2]]
    return(sl)
  }
  tmp = realised %>% 
    group_by(Stratum,Year,Route) %>% 
    summarise(mean = mean(log_expected),
              n = n()) %>% 
    group_by(Stratum) %>% 
    summarise(slp = slp(mean,Year),
              n = sum(n))
  
  tmpmap <- left_join(strata_map,tmp,by = "Stratum")
 
  pl = ggplot(data = tmpmap)+
    geom_sf(aes(fill = slp))+
    scale_fill_viridis_c()
print(pl)  
