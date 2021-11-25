### simulating fake BBS time-series data

library(tidyverse)
library(bbsBayes)
library(mgcv) #has functions to simulate data from a GAM


# select real BBS data for PAWR -------------------------------------------

BBS_data <- stratify("bbs_usgs")


for(species in c("Pacific Wren","Bobolink")){
  
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

if(species == "Pacific Wren"){ss <- 2021}
if(species == "Bobolink"){ss <- 2017}
set.seed(ss)
BETA_True <- rnorm(GAM_year$nknots_Year,0,1.5)

to_save <- c(to_save,"BETA_True")

mean_log_smooth <-  GAM_year$Year_basis %*% BETA_True

#plot(exp(mean_log_smooth),type = "l",ylim = c(0,max(exp(mean_log_smooth))))

## strata neighbourhoods Generate ---------------------------------
source("Functions/neighbours_define.R")
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
  mutate(Stratum_Factored = 1:nstrata)


strata_df <- strata_df %>% 
  left_join(.,st_coord,by = "Stratum_Factored")


neighbours <- neighbours_define(real_strata_map = strata_map,
                                plot_dir = "maps/",
                                species = "Simulated",
                                alt_strat = "Stratum")
to_save <- c(to_save,"neighbours")

## Generate stratum smooths ----------------------------------------

sd_spat_beta <- 0.3
to_save <- c(to_save,"sd_spat_beta")

## correlation matrix
neigh_mat <- neighbours$adj_matrix

## centre stratum
strat_mid <- which.max(colSums((neigh_mat)))
strata_df[strat_mid,]
to_save <- c(to_save,"strat_mid")


nstrata <- nrow(strata_df)
beta_True <- matrix(NA,nrow = nknots,
                    ncol = nstrata)
beta_True[,strat_mid] <- BETA_True

  wn = which(neigh_mat[strat_mid,] == 1)
  for(sj in wn){
    if(any(is.na(beta_True[,sj]))){
      for(b in 1:nknots){
        bm = mean(beta_True[b,strat_mid],na.rm = TRUE)
        beta_True[b,sj] <- rnorm(1,bm,sd_spat_beta)
      }
    }
  }
  
  while(any(is.na(beta_True))){
    wna <- which(is.na(beta_True[1,])) #strata that have no betas yet
    if(length(wna) > 1){
      ww <- which(!is.na(beta_True[1,])) # strata that do have betas
      wnxa <- (neigh_mat[ww,wna]) #matrix of columns of missing
      si <- wna[which.max(apply(wnxa,2,sum))]
    }else{
      si <- wna
    }
    wnx <- which(neigh_mat[,si] == 1)

      if(any(is.na(beta_True[,si]))){
        for(b in 1:nknots){
          bm = mean(beta_True[b,wnx],na.rm = TRUE)
          beta_True[b,si] <- rnorm(1,bm,sd_spat_beta)
        }
      }

    
  }
    
  to_save <- c(to_save,"beta_True")
  
  

## stratum-level smooths ----------------------------------------

  strat_log_smooths <- GAM_year$Year_basis %*% beta_True 

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

  pf <- ggplot(data = log_smooth_plot,aes(x = Year,y = index,colour = Y))+
    geom_line(size = 2)+
    scale_color_viridis_c()+
    scale_y_continuous(limits = c(0,NA))+
    facet_wrap(~Stratum_Factored,scales = "free_y",
               nrow = ceiling(sqrt(nstrata)),
               ncol = floor(sqrt(nstrata)))
  
  print(pf)
  
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
    facet_wrap(~Stratum_Factored,scales = "free_y",
               nrow = ceiling(sqrt(nstrata)),
               ncol = floor(sqrt(nstrata)))
  
  print(pf)
  
  

# Add observer, route, and strata intercepts ----------------------------

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

  to_save <- c(to_save,"observer_df")
  
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
              by = c("Route","Route_Factored","Stratum","Stratum_Factored"))
    
  
  ## Stratum Effects
  sdstrata <- 0.5
  
  strata_df <- strata_df %>% 
    mutate(True_strata_effects = rnorm(nstrata,0,sdstrata))
  
  to_save <- c(to_save,"strata_df")
  
  balanced <- balanced %>% 
    left_join(.,strata_df,
              by = c("X","Y","Stratum","Stratum_Factored"))
  
  
  
  

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
  left_join(.,balanced) %>% 
  mutate(Year_Index = Year-(min(Year)-1))


to_save <- c(to_save,"realised")



save(list = to_save,
     file = paste0("Simulated_data_",gsub(species,pattern = " ",replacement = "_"),"_BBS.RData"))


}
