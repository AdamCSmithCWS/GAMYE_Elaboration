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
  labs(title = paste(i,"By Route"))+
  facet_wrap(~version,nrow = 2,ncol = 1)

 }
 
    print(rt_plots[[1]]+ rt_plots[[2]]+ rt_plots[[3]]+ rt_plots[[4]])
             

## by strata ---------------------------------------------------------------


## by lat and long -------------------------------------------------------


## by year -----------------------------------------------------------------



          
                    
          
          
          
   
