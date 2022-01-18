### Fitting model to simulated BBS data
library(bbsBayes)
library(tidyverse)
# library(cmdstanr)
library(patchwork)

setwd("C:/Users/adam_/OneDrivedelete/Documents/GitHub/GAMYE_Elaboration")


#species <- "Pine Warbler"
species = "Yellow-headed Blackbird"  
species_f <- gsub(species,pattern = " ",replacement = "_")
 
   
          
          load(paste0("Data/Real_data_",species_f,"_BBS.RData"))
          #load(paste0("Data/Real_data_",species_f,"_obs_mix_BBS.RData"))
          

# compare count distributions means and sd by -----------------------------

# shared_names <- names(original_data_df)[which(names(original_data_df) != "count")]
#  original_data_df$true_count <- original_data_df$count         
realised$version <- "simulated"
original_data_df$version <- "true"          
          both <- realised %>% 
            bind_rows(original_data_df) %>% 
            mutate(count = log(count+1))
            
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
 
    #print(rt_plots[[1]]+ rt_plots[[2]]+ rt_plots[[3]]+ rt_plots[[4]])
             

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
    
    #print(st_plots[[1]]+ st_plots[[2]]+ st_plots[[3]]+ st_plots[[4]])
    #print(stlat_plots[[1]]+ stlat_plots[[2]]+ stlat_plots[[3]]+ stlat_plots[[4]])
    



# by observer -------------------------------------------------------------
rtmeans <- rt_sum %>% 
      select(Route,mean,version) %>% 
      rename(route_mean = mean)
    
    obs_sum <- both %>% 
      group_by(Observer,Route,version) %>% 
      summarise(mean = mean(count),
                sd = sd(count),
                min = min(count),
                max = max(count),
                n_count = n()) %>% 
      left_join(.,rtmeans) %>% 
      mutate(cent_mean = mean-route_mean)
    
    obs_plots <- vector(mode = "list",length = 8)
    names(obs_plots) <- c("cent_mean","mean","sd","min","max",
                          "mean_by_sd","sd_by_n","mean_by_n")
    for(i in names(obs_plots)[1:5]){
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

    obs_plots[[6]] = ggplot(data = tmp)+
      geom_point(aes(x = sd,y = mean),
                 #position = position_jitter(width = 0.1),
                 alpha = 0.03)+
      theme_bw()+
      labs(title = paste("sd by mean count by Observer"))+
      facet_wrap(~version,nrow = 2,ncol = 1)
    
    
    #print(obs_plots[[6]])
    
    obs_plots[[7]] = ggplot(data = tmp)+
      geom_point(aes(y = sd,x = n_count),
                 #position = position_jitter(width = 0.1),
                 alpha = 0.03)+
      theme_bw()+
      labs(title = paste("sd by n count by Observer"))+
      facet_wrap(~version,nrow = 2,ncol = 1)
    
    
    #print(obs_plots[[7]])
    
    obs_plots[[8]] = ggplot(data = tmp)+
      geom_point(aes(y = mean,x = n_count),
                 #position = position_jitter(width = 0.1),
                 alpha = 0.03)+
      theme_bw()+
      labs(title = paste("sd by n count by Observer"))+
      facet_wrap(~version,nrow = 2,ncol = 1)
    
    
    #print(obs_plots[[8]])
    
    # "mean_by_sd","sd_by_n","mean_by_n"
    # print(obs_plots[[1]]+ obs_plots[[2]]+ obs_plots[[3]]+ obs_plots[[5]])
    # 

    # 
    pdf(paste0("Figures/",species_f,"_sim_vs_real.pdf"),
        width = 11,
        height = 8.5)
    # pdf(paste0("Figures/",species_f,"_obs_mix_sim_vs_real.pdf"),
    #     width = 11,
    #     height = 8.5)
    print(rt_plots[[1]]+ rt_plots[[2]]+ rt_plots[[3]]+ rt_plots[[4]])
    print(st_plots[[1]]+ st_plots[[2]]+ st_plots[[3]]+ st_plots[[4]])
    print(stlat_plots[[1]]+ stlat_plots[[2]]+ stlat_plots[[3]]+ stlat_plots[[4]])
    
    print(obs_plots[[2]]+ obs_plots[[3]]+ obs_plots[[4]]+ obs_plots[[5]])
    print(obs_plots[[1]]+ obs_plots[[6]]+ obs_plots[[7]]+ obs_plots[[8]])
    dev.off()
    
          
                    
          
          
          
   
