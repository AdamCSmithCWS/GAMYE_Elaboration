### trends indices graphs and maps for real data

library(tidyverse)
library(cmdstanr)
library(posterior)
library(sf)
source("functions/indices_cmdstan.R")
source("functions/posterior_summary_functions.R")


fls <- data.frame(species_f = c("Yellow-headed_Blackbird",
                              "Zonotrichia_albicollis",
                              "Red_Knot"),
                  species = c("Yellow-headed Blackbird",
                              "White-throated Sparrow",
                              "Red Knot"),
                  data = c("BBS",
                           "CBC",
                           "Shorebird"),
                  out_base = c(paste0("Yellow-headed_Blackbird","_real_","BBS"),
                               paste0("Zonotrichia_albicollis","_CBC"),
                               paste0("Red Knot","_Shorebird")),
                  y1 = c(1966,
                         1966,
                         1980),
                  strat_map_name = c("Stratum_Factored",
                                     "",
                                     ""))




# Species loop ------------------------------------------------------------
output_dir <- "output/"

for(i in 1:nrow(fls)){

  species = fls[i,"species"]
  species_f <- fls[i,"species_f"]
  dd <- fls[i,"data"]
  out_base <- fls[i,"out_base"]
  year_1 = fls[i,"y1"]
  
  st_n = fls[i,"strat_map_name"]
  
  csv_files <- paste0(out_base,"-",1:3,".csv")
  
  load(paste0("Data/",species_f,dd,"_data.RData"))
  load(paste0(output_dir,"/",out_base,"_gamye_iCAR.RData"))
  #stan_data
  #neighbours
  #realized_strata_map
  #data_1
  
# Strata Annual indices file -----------------------------------------------------
  strat_df <- as.data.frame(realized_strata_map)
  
 strat_inds <- index_function(fit = stanfit,
                              parameter = "n",
                              year_1 = year_1,
                              strat = st_n)
  indices <- strat_inds$indices %>% 
    inner_join(.,strat_df)%>% 
    mutate(version = "full")

  strat_inds_smooth <- index_function(fit = stanfit,
                               parameter = "nsmooth",
                               year_1 = year_1,
                               strat = st_n)
  indices_smooth <- strat_inds_smooth$indices %>% 
    inner_join(.,strat_df) %>% 
    mutate(version = "smooth")
  
  indices_all <- bind_rows(indices,
                           indices_smooth) %>% 
    rename_with(.,~gsub(st_n,replacement = "strat_plot",
                        .x)) %>% 
    mutate(original_strat_name = strat_plot) %>% 
    rename_with(.,~gsub(replacement = st_n,pattern = "strat_plot",
                        .x))
  
  pd = ceiling(sqrt(length(unique(indices_all$original_strat_name))))  
  
  
  
  pl_inds <- ggplot(data = indices_all,aes(x = Year,y = median))+
    geom_ribbon(aes(ymin = lci,ymax = uci,fill = version),alpha = 0.1)+
    geom_line(aes(colour = version))+
    labs(title = species)+
    facet_wrap(vars(original_strat_name),
               nrow = pd,
               ncol = pd,
               scales = "free_y")
    
print(pl_inds)
  
  # Overall Annual indices









  
  
}#end species loop






