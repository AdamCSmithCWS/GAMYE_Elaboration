## trend map annimation

library(tidyverse)
library(gganimate)
library(sf)
library(cmdstanr)
library(gifski)

source("functions/indices_cmdstan.R")
source("functions/posterior_summary_functions.R")
source("Functions/palettes.R")


fls <- data.frame(species_f = c("Yellow-headed_Blackbird",
                                "Cinclus_mexicanus",
                                "Red_Knot"),
                  species = c("Yellow-headed Blackbird",
                              "Cinclus_mexicanus",
                              "Red Knot"),
                  data = c("BBS",
                           "CBC",
                           "Shorebird"),
                  out_base = c(paste0("Yellow-headed_Blackbird","_real_","BBS"),
                               paste0("Cinclus_mexicanus","_CBC_B"),
                               paste0("Red Knot","_Shorebird")),
                  y1 = c(1966,
                         1966,
                         1980),
                  strat_map_name = c("Stratum_Factored",
                                     "strata_vec",
                                     "stratn"))



output_dir <- "output/"


for(i in c(1:3)){#1:nrow(fls)){
  
  species = fls[i,"species"]
  species_f <- fls[i,"species_f"]
  dd <- fls[i,"data"]
  out_base <- fls[i,"out_base"]
  year_1 = fls[i,"y1"]
  
  st_n = fls[i,"strat_map_name"]
  
  
  
  load(paste0("Data/",species_f,dd,"_data.RData"))
  load(paste0(output_dir,out_base,"_gamye_iCAR.RData"))
  csv_files <- paste0(output_dir,out_base,"-",1:3,".csv")
  
  #stan_data
  #neighbours
  #realized_strata_map
  #data_1
  
  ## 
  
  
  stanfit <- as_cmdstan_fit(files = csv_files)
  strat_df <- as.data.frame(realized_strata_map)
  


strat_inds_smooth <- index_function(fit = stanfit,
                                    parameter = "nsmooth",
                                    year_1 = year_1,
                                    strat = st_n)

indices_smooth <- strat_inds_smooth$indices %>% 
  inner_join(.,strat_df) %>% 
  mutate(version = "smooth")


tyrs = unique(c(2018:year_1))

trends_out <- NULL

for(j in 1:length(tyrs)){
  yy = tyrs[j]
  yy2 = tyrs[j]+1
  
  tt <- trends_function(ind_list = strat_inds_smooth,
                        start_year = yy,
                        end_year = yy2) %>% 
    mutate(species = species,
           first_year = yy,
           last_year = yy2,
           region_type = "Stratum")
  
  ttmd <- realized_strata_map %>% 
    left_join(.,tt,by = st_n) %>% 
    mutate(trend_plot = trend_plot_cats(trend),
           species = species)
  
  trends_out = bind_rows(trends_out,ttmd)
  
  
}



tmp = trends_out %>% 
  filter(species == species)

prov_state <- bbsBayes::load_map(stratify_by = "state")

xb = range(st_coordinates(realized_strata_map)[,"X"])
yb = range(st_coordinates(realized_strata_map)[,"Y"])


angif <- ggplot(data = tmp)+
  geom_sf(data = prov_state,alpha = 0,
          colour = grey(0.8),inherit.aes = FALSE)+
  geom_sf(aes(fill = trend_plot))+
  coord_sf(xlim = xb,
           ylim = yb)+
  theme_bw()+
  xlab("")+
  ylab("")+
  scale_colour_manual(values = map_palette_s, 
                      aesthetics = c("fill"),
                      guide = guide_legend(reverse=TRUE),
                      name = paste0("Trend"))+
  transition_manual(frames = first_year)+
  ggtitle('Annual trend: {frame + 1966}',
          subtitle = paste(dd,species_f))

# save_animation(angif,file = paste0("Figures/Animated_trends_",species_f,".gif",))
animate(angif,nframes = length(tyrs), fps = 1,
        renderer = gifski_renderer(file = paste0("Figures/Animated_trends_",species_f,".gif")))
#anim_save(animation = angif,file = paste0("Figures/Animated_trends_",species_f,".gif"))




}




tt_map_data[[paste0("TY",yy,"-",yy2)]] <- ttmd
  
  
  xb = range(st_coordinates(realized_strata_map)[,"X"])
  yb = range(st_coordinates(realized_strata_map)[,"Y"])
  
  prov_state <- bbsBayes::load_map(stratify_by = "state")
  
  ttm  <- ggplot(data = ttmd)+
    geom_sf(data = prov_state,alpha = 0,
            colour = grey(0.8),inherit.aes = FALSE)+
    geom_sf(aes(fill = trend_plot))+
    coord_sf(xlim = xb,
             ylim = yb)+
    theme_bw()+
    scale_colour_manual(values = map_palette_s, 
                        aesthetics = c("fill"),
                        guide = guide_legend(reverse=TRUE),
                        name = paste0("Trend"))
  
  
  # print(ttm)
  
  
  tt_map[[paste0("TY",yy,"-",yy2)]] <- ttm
  
  stratum_trends <- bind_rows(stratum_trends,tt)
}
