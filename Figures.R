#### Figures for publication
library(tidyverse)
library(cmdstanr)
library(posterior)
library(sf)
source("functions/indices_cmdstan.R")
source("functions/posterior_summary_functions.R")
source("Functions/palettes.R")


# 1 map of example strata connections ---------------------------------
species = "Yellow-headed Blackbird"  
species_f <- gsub(species,pattern = " ",replacement = "_")
  
load(paste0("Simulated",species_f,"_route_maps.RData"))

# 2 BETA accuracy -----------------------------------------------------

tp <- "non_linear"

load(paste0("Data/Simulated_data_",species_f,"_",tp,"_BBS.RData"))

sns <- ""
mk <- ""
output_dir <- "output/"
out_base <- paste0(species_f,"_sim_",sns,mk,tp,"_BBS")
out_base_sim <- paste0("_sim_",sns,mk,tp,"_BBS")

load(paste0("data/",out_base_sim,"_accuracy_comp.RData"))


# 3 beta accuracy -----------------------------------------------------------



# 4 trajectory accuracy ---------------------------------------------------



# 5 masked strata - comparison of spatial and non-spatial -----------------



# 6 real data overall trajectories for 3 species --------------------------




# 7 long-term and short-term (3-gen) trend maps ---------------------------
# six panel, paired maps









