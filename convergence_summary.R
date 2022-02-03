### Fitting model to simulated BBS data
library(bbsBayes)
library(tidyverse)
library(cmdstanr)
library(posterior)
source("Functions/posterior_summary_functions.R")

#setwd("C:/Users/adam_/OneDrivedelete/Documents/GitHub/GAMYE_Elaboration")

#species <- "Pine Warbler"

species = "Yellow-headed Blackbird"  

species_f <- gsub(species,pattern = " ",replacement = "_")


# real data first ---------------------------------------------------------

output_dir <- "Output/"
out_base <- paste0(species_f,"_real_","BBS")
mod.file = "models/gamye_iCAR_bbs.stan"

csv_files <- paste0(out_base,"-",1:3,".csv")

load(paste0("Data/Real_data_",species_f,"_BBS.RData"))
load(paste0(output_dir,"/",out_base,"_gamye_iCAR.RData"))

csv_files <- paste0(output_dir,out_base,"-",1:3,".csv")
stanfit <- as_cmdstan_fit(files = csv_files)

# diag <- stanfit$cmdstan_diagnose()


stanf_df <- stanfit$draws(format = "df")

#removes the count-level parameters
drws <- stanf_df %>% 
  select(!contains(c("E[","noise_raw"))) 

summ <- summarise_draws(drws)
summ <- summ %>% 
  mutate(species = species,
         data = "Real_BBS")




# Simulated data ----------------------------------------------------------


for(tp in c("non_linear","linear")){

  load(paste0("Data/Simulated_data_",species_f,"_",tp,"_BBS.RData"))
  
  for(sns in c("","nonSpatial_"))
    for(mk in c("","mask_")){
      
      if(sns == "nonSpatial_" & mk == "mask_"){next}
      
      
      output_dir <- "output/"
      #out_base <- paste0(species_f,"_sim_",tp,"_BBS")
      out_base <- paste0(species_f,"_sim_",sns,mk,tp,"_BBS")
      #csv_files <- paste0(out_base,"-",1:3,".csv")
      out_base_sim <- paste0("_sim_",sns,mk,tp,"_BBS")
      
      load(paste0(output_dir,"/",out_base,"_gamye_iCAR.RData"))
      
      csv_files <- paste0(output_dir,out_base,"-",1:3,".csv")
      stanfit <- as_cmdstan_fit(files = csv_files)
      

    
    
    
    stanf_df <- stanfit$draws(format = "df")
    
    #removes the count-level parameters
    drws <- stanf_df %>% 
      select(!contains(c("E[","noise_raw"))) 
    
    summt <- summarise_draws(drws)
    
    # sdsumm <- summt %>% 
    #   filter(grepl("sd",x = variable))
    # 
    summt <- summt %>% 
      mutate(species = species,
             data = out_base_sim)
    
    summ <- bind_rows(summ,summt)
    
    print(out_base_sim)
} 
}
    failed_rhat <- summ %>% 
      filter(rhat > 1.05)
    
    failed_ess_bulk <- summ %>% 
      filter(ess_bulk < 100)
save(list = c("summ"),
     file = paste0("output/convergence_summary_",species_f,".RData"))    
    