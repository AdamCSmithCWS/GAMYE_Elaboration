
library(tidyverse)
library(bbsBayes)
source("functions/neighbours_define_alt.R")
species = "Zonotrichia_albicollis"

data_1 = read.csv(paste0("data/",species,"_modeled_records.csv"))

data_1 <- data_1 %>% 
  mutate(strat = paste(country,state,bcr,sep = "-"))



strata_map <- bbsBayes::load_map(stratify_by = "bbs_usgs") %>% 
  rename(strat = ST_12) 

strat_df = data_1 %>% 
  select(contains("strat"),
         bcr,state,country) %>% 
  distinct()

realized_strata_map <- strata_map %>% 
  right_join(.,strat_df,by = "strat")


neighbours = neighbours_define(realized_strata_map,
                  species = species,
                  alt_strat = "strata_vec",
                  plot_dir = "maps/",
                  plot_file = "_strata_map")








