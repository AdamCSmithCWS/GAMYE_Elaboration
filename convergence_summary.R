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

out_base <- paste0(species_f,"_real_","BBS")
mod.file = "models/gamye_iCAR_bbs.stan"

csv_files <- paste0(out_base,"-",1:3,".csv")

load(paste0("Data/Real_data_",species_f,"_BBS.RData"))
load(paste0(output_dir,"/",out_base,"_gamye_iCAR.RData"))

# diag <- stanfit$cmdstan_diagnose()


stanf_df <- stanfit$draws(format = "df")

#removes the count-level parameters
drws <- stanf_df %>% 
  select(!contains(c("E[","noise_raw"))) 

summ <- summarise_draws(drws)


# Simulated data ----------------------------------------------------------



for(tp in c("non_linear","linear")){
  
  #STRATA_True <- log(2)
  output_dir <- "output/"
  out_base <- paste0(species_f,"_sim_",tp,"_BBS")
  csv_files <- paste0(out_base,"-",1:3,".csv")
  
    load(paste0("Data/Simulated_data_",species_f,"_",tp,"_BBS.RData"))
    load(paste0(output_dir,"/",out_base,"_gamye_iCAR.RData"))
