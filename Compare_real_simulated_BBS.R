### Fitting model to simulated BBS data
library(bbsBayes)
library(tidyverse)
library(cmdstanr)
library(patchwork)

setwd("C:/Users/adam_/OneDrivedelete/Documents/GitHub/GAMYE_Elaboration")


species <- "Pine Warbler"
  species_f <- gsub(species,pattern = " ",replacement = "_")
 
   
          
          load(paste0("Data/Real_data_",species_f,"_BBS.RData"))
 

# compare count distributions means and sd by -----------------------------

# shared_names <- names(original_data_df)[which(names(original_data_df) != "count")]
#  original_data_df$true_count <- original_data_df$count         
realised$version <- "simulated"
original_data_df$version <- "true"          
          both <- realised %>% 
            bind_rows(original_data_df) 
            
## by route ----------------------------------------------------------------
rt_sum <- both %>% 
  group_by(Route,Stratum,version) %>% 
  summarise(mean = mean(count),
            sd = sd(count),
            min = min(count),
            max = max(count))

    rt_plots <- vector(mode = "list",length = 4)
    names(rt_plots) <- c("mean","sd","min","max")
 for(i in names(rt_plots)){
   rn <- function(x){
     gsub(pattern = i,
          replacement = "tc",
          x = x)
   }
   tmp <- rt_sum %>% 
     rename_with(.fn = rn)
   
   
rt_plots[[i]] = ggplot(data = tmp)+
  geom_histogram(aes(x = tc))+
  labs(title = paste(i,"count by Route"))+
  facet_wrap(~version,nrow = 2,ncol = 1)

 }
 
    print(rt_plots[[1]]+ rt_plots[[2]]+ rt_plots[[3]]+ rt_plots[[4]])
             

## by strata ---------------------------------------------------------------

    st_sum <- both %>% 
      group_by(Stratum,version) %>% 
      summarise(mean = mean(count),
                sd = sd(count),
                min = min(count),
                max = max(count)) %>% 
      left_join(.,strata_df,by = "Stratum")
    
    st_plots <- vector(mode = "list",length = 4)
    names(st_plots) <- c("mean","sd","min","max")
    
    stlat_plots <- st_plots
    
    for(i in names(st_plots)){
      rn <- function(x){
        gsub(pattern = i,
             replacement = "tc",
             x = x)
      }
      tmp <- st_sum %>% 
        rename_with(.fn = rn)
      
      
      st_plots[[i]] = ggplot(data = tmp)+
        geom_histogram(aes(x = tc))+
        labs(title = paste(i,"count by Stratum"))+
        facet_wrap(~version,nrow = 2,ncol = 1)
      
      stlat_plots[[i]] = ggplot(data = tmp)+
        geom_point(aes(x = X,y = Y,colour = tc))+
        scale_colour_viridis_c()+
        labs(title = paste(i,"count by Stratum"))+
        facet_wrap(~version,nrow = 2,ncol = 1)
      
    }
    
    print(st_plots[[1]]+ st_plots[[2]]+ st_plots[[3]]+ st_plots[[4]])
    print(stlat_plots[[1]]+ stlat_plots[[2]]+ stlat_plots[[3]]+ stlat_plots[[4]])
    



# by observer -------------------------------------------------------------

    obs_sum <- both %>% 
      group_by(Observer,version) %>% 
      summarise(mean = mean(count),
                sd = sd(count),
                min = min(count),
                max = max(count))
    
    obs_plots <- vector(mode = "list",length = 4)
    names(obs_plots) <- c("mean","sd","min","max")
    for(i in names(obs_plots)){
      rn <- function(x){
        gsub(pattern = i,
             replacement = "tc",
             x = x)
      }
      tmp <- obs_sum %>% 
        rename_with(.fn = rn)
      
      
      obs_plots[[i]] = ggplot(data = tmp)+
        geom_histogram(aes(x = tc))+
        labs(title = paste(i,"count by Observer"))+
        facet_wrap(~version,nrow = 2,ncol = 1)
      
    }
    
    print(obs_plots[[1]]+ obs_plots[[2]]+ obs_plots[[3]]+ obs_plots[[4]])
    

          
                    
          
          
          
   
