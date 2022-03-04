### simulating fake BBS time-series data
# library(rsoi)
# pdo <- download_pdo()
# save(list = c("pdo"),file = "data/pdo.RData")

library(tidyverse)
library(bbsBayes)
library(mgcv) #has functions to simulate data from a GAM

# select real BBS data for CEWA -------------------------------------------

BBS_data <- stratify("bbs_usgs")


#species = "Pine Warbler"  
species = "Yellow-headed Blackbird"  
species_f <- gsub(species,pattern = " ",replacement = "_")
  
  tp = "breakpoint_cycle"
    
source("Functions/prepare-jags-data-alt.R")

real_data <- prepare_jags_data_alternate(strat_data = BBS_data,
                          species_to_run = species,
                          model = "gamye",
                          min_year = NULL) #n-year time-series ???
min_year <- min(real_data$r_year)

to_save1 <- c("real_data","min_year")

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

to_save1 <- c(to_save1,"real_df")


# dataframe of routes and strata
routes_df <- real_df %>% 
  select(Route,Route_Factored,Stratum,Stratum_Factored) %>% 
  distinct() 
to_save1 <- c(to_save1,"routes_df")

# dataframe of first year for each observer on a route
observer_route_df <- real_df %>% 
  select(Route,Route_Factored,Observer,Year,First_Year) %>%
  filter(First_Year == 1) %>% 
  distinct() %>% 
  arrange(Route,Year)
to_save1 <- c(to_save1,"observer_route_df")

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
to_save1 <- c(to_save1,"balanced")

# Generate mean smooth ----------------------------------------
source("Functions/GAM_basis_function_mgcv.R")
years_df <- data.frame(Year = min(balanced$Year):max(balanced$Year))

nknots = floor(nrow(years_df)/4)

GAM_year <- gam_basis(years_df$Year,
                      nknots = nknots,
                      sm_name = "Year")


to_save1 <- c(to_save1,"GAM_year")

# add a non-linear trend - broken stick with a single breakpoint and a
# latitudinal variation in the first and last slopes
Y_change = 1995
SLOPE_1 = 0
SLOPE_2 = -2

# a cyclical pattern similar to something climatic (~approximately 5 year cycle)
# shift the cycle intensity from east to west
# 
### using the a loess smooth of the realised time-series of the pacific decadal oscillation
### span of the smooth set to 10 years
load("data/pdo.RData")
pdo_ann <- pdo %>% 
  group_by(Year) %>% 
  summarise(mean_pdo = mean(PDO)) %>% 
  filter(Year > 1965)
sp = 10/nrow(pdo_ann)
sm = loess(mean_pdo~Year,
           data = pdo_ann,
           span = sp)
CYCLE <- predict(sm)## CYCLE is now a mean smoothed time-series to add to the main breakpoint trends
to_save1 <- c(to_save1,"CYCLE","pdo_ann",
              "Y_change","SLOPE_1","SLOPE_2")

years_df <- data.frame(Year = min(balanced$Year):max(balanced$Year))

## strata neighbourhoods Generate ---------------------------------
nstrata <- real_data$nstrata


strata_df <- balanced %>% 
  select(Stratum,Stratum_Factored) %>% 
  distinct() %>% 
  arrange(Stratum_Factored)

to_save1 <- c(to_save1,"strata_df","nstrata")


strata_map <- bbsBayes::load_map(stratify_by = "bbs_usgs") %>% 
  rename(Stratum = ST_12) %>% 
  right_join(.,strata_df,by = "Stratum") %>% 
  arrange(Stratum_Factored) ### this arranging is critical to the correct neighbourhoods

to_save1 <- c(to_save1,"strata_map")

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
to_save1 <- c(to_save1,"neighbours")

# neighbours2 <- neighbours_define(real_strata_map = strata_map,
#                                  plot_dir = "maps/",
#                                  species = paste0("Simulated_voronoi",species_f),
#                                  alt_strat = "Stratum",
#                                  voronoi = TRUE)
# to_save1 <- c(to_save1,"neighbours2")

## Generate stratum smooths and intercepts ----------------------------------------

### use simple, smooth, x-y coordinate variation in betas and stratas

strata_df <- strata_df %>% 
  mutate(yscale = as.numeric(scale(Y,scale = TRUE)),
         xscale = as.numeric(-1*scale(X,scale = TRUE)),
         sumxy = yscale+xscale,
         xy = yscale*xscale,
         y_change = round(Y_change+yscale),
         slope_1 = log(((SLOPE_1+yscale*-0.5)/100)+1),
         slope_2 = log(((SLOPE_2+SLOPE_1+yscale*-0.5+xscale*-0.5)/100)+1),
         cycle_str = (0.02*exp(0.6*(max(yscale)-(yscale)))))



# strat_tempplot <- ggplot(data = strata_df,aes(x = X,y = Y))+
#   geom_point(aes(colour = sumxy),size = 3)+
#   scale_colour_viridis_c()
# 
# print(strat_tempplot)




## correlation matrix
neigh_mat <- neighbours$adj_matrix


save(list = c("to_save1",to_save1),
     file = paste0("data/","simulation_",tp,"_basic_data.RData"))

#MAs <- round(log(c(1,5,10,20,50)),2)
MAs <- round(log(c(0.1,0.5,1,5,10,20,50)),2)# true mean abundances for different simulations

# Loop Mean Abundance -----------------------------------------------------


for(ma in MAs){
  
  load(paste0("data/","simulation_",tp,"_basic_data.RData"))
  
  to_save <- c(to_save1,"strata_True",
               "observer_df","routes_df","strata_df","realised",
               "mask_map","strata_mask","realised_mask",
               "routes_mask","event_mask_retain",
               "log_true_traj")
### strata intercepts
nstrata <- nrow(strata_df)
STRATA_True <- ma #mean abundance
strata_df <- strata_df %>% 
  mutate(strata_True = STRATA_True - abs(yscale)) # abundance peaks at middle latitudes
strata_True <- as.numeric(strata_df$strata_True)

smooth_func <- function(my = 1992,
                        cy = 1995,
                        y = 1966,
                        s1 = 0,
                        s2 = 0.1,
                        cyc = 0.1,
                        cyc_r = 0.05){
  
  sm = ifelse(y <= cy,
              s1*(y-cy)+(cyc_r*cyc*2),
              s2*(y-cy)+(cyc_r*cyc*2))
  
  return(sm)
}

### strata trajectories
log_smooth_plot <- expand_grid(Year = 1966:2019,
                               Stratum_Factored = strata_df$Stratum_Factored) %>% 
  left_join(.,strata_df,by = "Stratum_Factored") %>% 
  mutate(cycle = CYCLE[Year-1965],
         Smooth = smooth_func(my = 1992,
                              cy = y_change,
                              y = Year,
                              s1 = slope_1,
                              s2 = slope_2,
                              cyc = cycle,
                              cyc_r = cycle_str),
         index = exp(Smooth))



 
## stratum-level smooths ----------------------------------------

ch_years <- log_smooth_plot %>% 
  filter(Year == y_change)
  pfs <- ggplot(data = log_smooth_plot,
                aes(x = Year,y = index,colour = Y))+
    geom_point(data = ch_years,colour = "red",size = 0.2)+
    geom_line(size = 1)+
    scale_color_viridis_c()+
    geom_point(data = ch_years,colour = "red")+
    #scale_y_continuous(limits = c(0,NA))+
    facet_wrap(~Stratum,scales = "fixed",
               nrow = ceiling(sqrt(nstrata)),
               ncol = ceiling(sqrt(nstrata)))
  
  print(pfs)
  
  ## Add random annual fluctuations ----------------------------------------

  sdyeareffect <- 0.1 # ~10% mean annual fluctuation
  
  ye_funct <- function(x,sd = sdyeareffect){
    ye = rnorm(length(x),0,sd)
  }
  
  
  log_true_traj <- log_smooth_plot %>% 
    group_by(Stratum) %>% 
    mutate(YearEffect = ye_funct(Smooth),
           True_log_traj = Smooth + YearEffect,
           True_traj = exp(True_log_traj),
           True_scaled_smooth = exp(Smooth + strata_True))
  

  pf <- ggplot(data = log_true_traj,aes(x = Year,y = True_traj,colour = Y))+
    geom_line(size = 2)+
    scale_color_viridis_c()+
    scale_y_continuous(limits = c(0,NA))+
    facet_wrap(~Stratum,scales = "free_y",
               nrow = ceiling(sqrt(nstrata)),
               ncol = ceiling(sqrt(nstrata)))
  
  
  strat_grid <- geofacet::grid_auto(strata_map,
                                    codes = "Stratum_Factored",
                                    names = "Stratum",
                                    seed = 2019)
  
  g_inds <- ggplot(data = log_true_traj,aes(x = Year,
                                            y = True_scaled_smooth,
                                            colour = Y))+
    geom_line(size = 1)+
    scale_color_viridis_c()+
    scale_y_continuous(limits = c(0,NA))+
    geofacet::facet_geo(~Stratum,grid = strat_grid,
                      scales = "free")+
    theme(strip.text = element_text(size = 6),
          strip.background = element_blank(),
          panel.spacing = unit(0.1,"mm"),
          axis.text.x = element_text(size = 5))

  
  
  


# Add observer, route intercepts ----------------------------

  balanced <- balanced %>% left_join(.,log_true_traj,
                                     by = c("Stratum","Year","Stratum_Factored")) %>% 
    select(-c(index,True_traj))
  
  ## OBserver Effects
  sdobs <- 0.2
  nobservers <- length(unique(balanced$Observer))
  True_observer_effects <- rnorm(nobservers,0,sdobs) #True observer effects by 
  
  observer_df <- data.frame(True_observer_effects = True_observer_effects,
                            Observer = unique(balanced$Observer)) %>% 
    mutate(Observer_Factored = as.integer(factor(Observer)))

  
  balanced <- balanced %>% 
    left_join(observer_df,by = "Observer")
  
  ## Route Effects
  sdroute <- 0.2
  nroutes <- max(routes_df$Route_Factored)
  
  routes_df <- routes_df %>% 
    mutate(True_route_effects = rnorm(nroutes,0,sdroute))
  
  balanced <- balanced %>% 
    left_join(.,routes_df,
              by = c("Route",
                     "Route_Factored",
                     "Stratum",
                     "Stratum_Factored"))
    
  
  ## Stratum Effects
  
  strata_df <- strata_df %>% 
    mutate(True_strata_effects = strata_True)
  
  
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

mean_exp_count <- balanced %>% 
  group_by(Stratum_Factored,Year) %>% 
  summarise(mean_exp = mean(expected_count))

mean_ec_plot <- ggplot(data = mean_exp_count,
                       aes(x = Year,y = mean_exp))+
  geom_line()+
  scale_y_continuous(limits = c(0,NA))+
  facet_wrap(vars(Stratum_Factored),
             scales = "free_y")



pdf(paste0("Figures/Simulated_",ma,"_",tp,"_True_smooth.pdf"),
    width = 11,
    height = 8.5)
print(pfs)
print(pf)
print(mean_ec_plot)
print(g_inds)

dev.off()


realised <- real_df %>% select(-c("Observer")) %>% 
  left_join(.,balanced,
            by = c("Year",
                   "Stratum",
                   "Stratum_Factored",
                   "Route",
                   "Route_Factored")) %>% 
  mutate(Year_Index = Year-(min(Year)-1))



if(ma == MAs[1]){
  # original data with updated observer adn route info ----------
original_data_df <- data.frame(Year = real_data$r_year,
                        Stratum = real_data$strat_name,
                        Stratum_Factored = real_data$strat,
                        Observer = real_data$ObsN,
                        Route = real_data$route,
                        First_Year = real_data$firstyr,
                        count = real_data$count) %>% 
  mutate(Route_Factored = as.integer(factor(Route)),
         Observer_Factored = as.integer(factor(real_data$ObsN)))

original_data_df <- original_data_df %>% 
  mutate(Year_Index = Year-(min(Year)-1))


}


# Mask data in some strata ------------------------------------------------
set.seed(1)

masked_strata <- sample(1:nstrata,size = floor(ceiling(nstrata/4)),replace = FALSE)
strata_mask <- strata_df %>% 
  mutate(masked = ifelse(Stratum_Factored %in% masked_strata,TRUE,FALSE))
  
mask_map <- strata_map %>% 
  mutate(masked = ifelse(Stratum_Factored %in% masked_strata,TRUE,FALSE))


# testplot = ggplot()+
#   geom_sf(data = mask_map,aes(fill = masked))
# print(testplot)

# remove all but two years from each observer and route in some st --------

routes_mask <- routes_df %>% 
  filter(Stratum_Factored %in% masked_strata)

# minimum number of observation events to retain so that all routes and observers are included

event_mask_retain <- realised %>% 
  filter(Route %in% routes_mask$Route) %>% 
  group_by(Route,Observer,First_Year) %>% 
  sample_n(size = 1) %>% 
  as.data.frame()

realised_mask <- realised %>% 
  filter(!Route %in% routes_mask$Route) %>% 
  bind_rows(.,event_mask_retain)
  


# tst1 = realised_mask %>% 
#   group_by(Stratum_Factored) %>% 
#   summarise(n_m = n())
# tst2 = realised %>% 
#   group_by(Stratum_Factored) %>% 
#   summarise(n_r = n())
# tst = inner_join(tst1,tst2)
# 
# mask_map <- mask_map %>% 
#   left_join(.,tst)
# 

# 
# testplot1 = ggplot()+
#    geom_sf(data = mask_map,aes(fill = masked))
# testplot2 = ggplot()+
#   geom_sf(data = mask_map,aes(fill = n_m/n_r))
# # testplot3 = ggplot()+
# #   geom_sf(data = mask_map,aes(fill = n_m))
# #library(patchwork)
#    print(testplot1 + testplot2)
#   
### realised is a final dataset to demonstrate a realistic BBS dataset
### with a known pattern of population change

### realised_mask is an equivalent dataset that has data mostly missing
### from one portion of the range 
### mostly missing means only 1-years observations are available from each route ad observer
### so that the route and observer structure is retained, but there is no
### information available on change over time, after the observer and route intercepts


to_save <- unique(to_save)

save(list = to_save,
     file = paste0("Data/Simulated_data_",ma,"_",tp,"_BBS.RData"))

if(ma == MAs[1]){
  save(list = c(to_save,"original_data_df"),
       file = paste0("Data/Real_data_",species_f,"_BBS.RData"))
  
}

rm(list = to_save)


}#end ma loop
