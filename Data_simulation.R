### simulating fake BBS time-series data

library(tidyverse)
library(bbsBayes)
library(mgcv) #has functions to simulate data from a GAM


# select real BBS data for PAWR -------------------------------------------

BBS_data <- stratify("bbs_usgs")


need to generate True observer information (not the observer-route combinations)

real_data <- prepare_data(strat_data = BBS_data,
                          species_to_run = "Pacific Wren",
                          model = "gamye",
                          min_year = NULL) #n-year time-series ???
min_year <- min(real_data$r_year)

# Generate balanced dataset -----------------------------------------------

# create full data frame 
real_df <- get_prepared_data(real_data)#generates the dataframe of counts etc.

# dataframe of routes and strata
routes_df <- real_df %>% 
  select(Route,Stratum,Stratum_Factored) %>% 
  distinct() 

# dataframe of first year for each observer on a route
observer_route_df <- real_df %>% 
  select(Route,Observer_Factored,Year,First_Year) %>%
  filter(First_Year == 1) %>% 
  distinct() %>% 
  arrange(Route,Year)

# Generate balanced dataset with each year added in the realised first year it was surveyed
balanced <- NULL
for(rr in unique(observer_route_df$Route)){
  tmp <- observer_route_df %>% filter(Route == rr) %>% 
    arrange(Year)
  if(nrow(tmp) > 1){
  for(i in 1:(nrow(tmp)-1)){
    tmp1 <- data.frame(Route = rr,
                       Observer_Factored = tmp[i,"Observer_Factored"],
                       Year = tmp[i,"Year"]:(tmp[i+1,"Year"]-1))
    balanced <- bind_rows(balanced,tmp1)
  }
  }else{
    i = 1
  }
  tmp1 <- data.frame(Route = rr,
                     Observer_Factored = tmp[i,"Observer_Factored"],
                     Year = tmp[i,"Year"]:2019)
  balanced <- bind_rows(balanced,tmp1)
  
  
}


# Fill in strata info on balanced
balanced <- balanced %>% 
  left_join(.,routes_df,by = c("Route"))


# Generate mean smooth ----------------------------------------

source("Functions/GAM_basis_function_mgcv.R")
years_df <- data.frame(Year = min(balanced$Year):max(balanced$Year))

nknots = floor(nrow(years_df)/4)

GAM_year <- gam_basis(years_df$Year,
                      nknots = nknots,
                      sm_name = "Year")

set.seed(2021)
BETA_True <- rnorm(GAM_year$nknots_Year,0,1.5)

mean_log_smooth <-  GAM_year$Year_basis %*% BETA_True

plot(exp(mean_log_smooth),type = "l",ylim = c(0,max(exp(mean_log_smooth))))

## strata neighbourhoods Generate ---------------------------------
source("Functions/neighbours_define.R")

strata_df <- balanced %>% 
  select(Stratum,Stratum_Factored) %>% 
  distinct() %>% 
  arrange(Stratum_Factored)


strata_map <- bbsBayes::load_map(stratify_by = "bbs_usgs") %>% 
  rename(Stratum = ST_12) %>% 
  right_join(.,strata_df,by = "Stratum") %>% 
  arrange(Stratum_Factored) ### this arranging is critical to the correct neighbourhoods

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

## Generate stratum smooths ----------------------------------------

sd_spat_beta <- 0.3
## correlation matrix
neigh_mat <- neighbours$adj_matrix

## centre stratum
strat_mid <- which.max(colSums((neigh_mat)))
strata_df[strat_mid,]


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
                                     by = c("Stratum","Year"))
  
  sdobs <- 0.2
  nobservers <- length(unique(balanced$Observer_Factored))
  True_observers <- rnorm()

