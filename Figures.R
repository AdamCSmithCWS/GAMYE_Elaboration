#### Figures for publication
library(tidyverse)
library(cmdstanr)
library(posterior)
library(sf)
library(patchwork)
source("functions/indices_cmdstan.R")
source("functions/posterior_summary_functions.R")
source("Functions/palettes.R")


# 1 map of example strata connections ---------------------------------
species = "Yellow-headed Blackbird"  
species_f <- gsub(species,pattern = " ",replacement = "_")
  
load(paste0("maps/Simulated",species_f,"_route_maps_data.RData"))


xb = range(st_coordinates(real_strata_map)[,"X"])
yb = range(st_coordinates(real_strata_map)[,"Y"])

set.seed(1)
real_strata_map <- real_strata_map %>% 
  mutate(rand_strat = sample(1:nrow(real_strata_map)))

ggp <- ggplot()+ 
  geom_sf(data = real_strata_map,
          aes(fill = rand_strat),
          alpha = 1,
          colour = grey(0.8))+
  geom_segment(data=DA,
               aes(x = long, y = lat,
                   xend=long_to,yend=lat_to),
               inherit.aes = FALSE,
               size=0.5,alpha=0.2) +
  geom_sf(data = centres, alpha = 0.9,colour = "white",
          size = 0.5) + 
  xlab("")+
  ylab("")+
  scale_fill_viridis_c()+
  theme_minimal() +
    coord_sf(xlim = xb,ylim = yb)+
    theme(legend.position = "none")

#print(ggp)


pdf(file = paste0("Figures/Figure_1.pdf"),
    width = 3.5,
    height = 4)
print(ggp)
dev.off()




# 2 BETA accuracy -----------------------------------------------------

tp <- "non_linear"

load(paste0("Data/Simulated_data_",species_f,"_",tp,"_BBS.RData"))

sns <- ""
mk <- ""
output_dir <- "output/"
out_base <- paste0(species_f,"_sim_",sns,mk,tp,"_BBS")
out_base_sim <- paste0("_sim_",sns,mk,tp,"_BBS")

load(paste0("data/",out_base_sim,"_accuracy_comp.RData"))

fig2 <- ggplot(data = BETA_comp)+
  geom_point(aes(x = k,y = True_BETA),
             size = 1)+
  geom_errorbar(aes(x = k,y = mean,ymin = lci,ymax = uci),
                width = 0, 
                colour = grey(0.8),
                alpha = 0.5)+
  geom_point(aes(x = k,y = mean),
             colour = grey(0.8),
             alpha = 0.75)+
  xlab("Knot position")+
  ylab("Hyperparameters")+
  theme_classic()

  
pdf(file = paste0("Figures/Figure_2.pdf"),
    width = 3.5,
    height = 3)
print(fig2)  
dev.off()



# 3 beta accuracy -----------------------------------------------------------
nstrata = length(unique(beta_comp$Stratum_Factored))

fig3 <- ggplot(data = beta_comp)+
  geom_point(aes(x = k,y = beta_True),
             size = 1)+
  geom_errorbar(aes(x = k,y = mean,ymin = lci,ymax = uci),
                width = 0,
                colour = grey(0.8),
                alpha = 0.5)+
  geom_point(aes(x = k,y = mean),
             colour = grey(0.8),
             alpha = 0.75)+
  xlab("Knot position")+
  ylab("Parameters")+
  theme_classic()+
  facet_wrap(vars(Stratum_Factored),ncol = ceiling(sqrt(nstrata)),
             nrow = ceiling(sqrt(nstrata)),
             scales = "free")
  


pdf(file = paste0("Figures/Figure_3.pdf"),
    width = 7,
    height = 8)

print(fig3)  
dev.off()



# 4 trajectory accuracy ---------------------------------------------------
nsmooth_comp2 <- nsmooth_comp2 %>% 
  left_join(.,obs_means, by = c("Stratum_Factored","Year"))


fig4 = ggplot(data = nsmooth_comp2,aes(y = True_nsmooth,
                                                x = Year))+
  geom_point(data = nsmooth_comp2,aes(x = Year,y = mean_count),
             alpha = 0.1,
             size = 0.2,
             inherit.aes = FALSE)+
  geom_ribbon(aes(ymin = Q_025,ymax = Q_975,fill = version),alpha = 0.2)+
  geom_line(aes(colour = version))+
  scale_colour_viridis_d(aesthetics = c("colour","fill"))+
  scale_y_continuous(limits = c(0,NA))+
  facet_wrap(~Stratum_Factored,nrow = ceiling(sqrt(nstrata)),
             ncol = ceiling(sqrt(nstrata)),
             scales = "free_y")+
  theme_classic() +
  theme(legend.position = "none")
print(fig4)


# 6 Trend comparisons ---------------------------------------------------

sw_trends <- NULL
strat_trends <- NULL
tp <- "non_linear"
tp <- "linear"
mk <- "mask_"
output_dir <- "output/"
load(paste0("Data/Simulated_data_",species_f,"_",tp,"_BBS.RData"))

strat_df <- as.data.frame(strata_mask) %>% 
  select(Stratum_Factored,masked)

for(sns in c("","nonSpatial_alt_")){#,"nonSpatial_"))
  for(mk in c("","mask_")){
    out_base_sim <- paste0("_sim_",sns,mk,tp,"_BBS")
    
    lbl <- "Spatial"
    if(sns == "nonSpatial_alt_"){
      lbl <- "NonSpatial"
    }
    if(mk != ""){
    lbl <- paste0(lbl," Masked")
      }
   
    
    load(paste0("data/",out_base_sim,"_accuracy_comp.RData"))
    sw_t <- all_trends %>% 
      filter(Region_type == "Survey_Wide_Mean",
             last_year == 2019) %>% 
      mutate(version = lbl)
    sw_trends <- bind_rows(sw_trends,sw_t)
    
    strat_t <- all_trends %>% 
      filter(Region_type == "Stratum_Factored",
             last_year == 2019) %>% 
      mutate(version = lbl)
    strat_trends <- bind_rows(strat_trends,strat_t)
  }
}

trends_plot <- ggplot(data = sw_trends,aes(x = true_trend,y = trend))+
  geom_errorbar(aes(ymin = lci,ymax = uci,colour = first_year),width = 0,
                alpha = 0.4)+
  geom_point(aes(colour = first_year))+
  geom_abline(slope = 1,intercept = 0)+
  facet_wrap(vars(version),
             nrow = 2,
             ncol = 2)

print(trends_plot)

strat_trends <- strat_trends %>% 
  left_join(.,strat_df,by = "Stratum_Factored")

trends2_plot <- ggplot(data = strat_trends,aes(x = true_trend,y = trend))+
  geom_errorbar(aes(ymin = lci,ymax = uci,colour = first_year),width = 0,
                alpha = 0.4)+
  geom_point(aes(colour = first_year))+
  geom_abline(slope = 1,intercept = 0)+
  facet_wrap(vars(version),
             nrow = 2,
             ncol = 2)

print(trends2_plot)

strat_trends <- strat_trends %>% 
  mutate(t_dif = true_trend - trend,
         t_abs_dif = abs(t_dif))

trends3_plot <- ggplot(data = strat_trends,aes(x = first_year,y = t_abs_dif,group = version))+
  # geom_errorbar(aes(ymin = lci,ymax = uci,colour = first_year),width = 0,
  #               alpha = 0.4)+
  geom_point(aes(colour = version),
             position = position_dodge(width = 3))
print(trends3_plot)

strat_trends_dif <- strat_trends %>% 
  group_by(first_year,version,masked) %>% 
  summarise(mean_abs_dif = mean(t_abs_dif,na.rm = T),
            median_abs_dif = median(t_abs_dif,na.rm = T),
            sd_abs_dif = sd(t_abs_dif,na.rm = T))

strat_trends_dif_nm <- strat_trends_dif %>% 
  filter(masked == FALSE)
trends4_plot <- ggplot(data = strat_trends_dif_nm,
                       aes(x = first_year,y = mean_abs_dif,group = version))+
  # geom_errorbar(aes(ymin = lci,ymax = uci,colour = first_year),width = 0,
  #               alpha = 0.4)+
  geom_point(aes(colour = version),
             position = position_dodge(width = 3))
print(trends4_plot)


strat_m_trends_dif <- strat_trends_dif %>% 
  filter(masked == TRUE)
trends5_plot <- ggplot(data = strat_m_trends_dif,
                       aes(x = first_year,y = mean_abs_dif,group = version))+
  # geom_errorbar(aes(ymin = lci,ymax = uci,colour = first_year),width = 0,
  #               alpha = 0.4)+
  geom_point(aes(colour = version),
             position = position_dodge(width = 3))
print(trends5_plot)


# 5 masked strata - comparison of spatial and non-spatial -----------------



# 6 real data overall trajectories for 3 species --------------------------




# 7 long-term and short-term (3-gen) trend maps ---------------------------
# six panel, paired maps

load("output/real_data_summaries.RData")


fls <- data.frame(species_f = c("Yellow-headed_Blackbird",
                                "Cinclus_mexicanus",
                                "Red_Knot"),
                  species = c("Yellow-headed Blackbird",
                              "American Dipper",
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



t1 <- tt_map_list[[1]][["TY1966-2019"]]  +
  labs(subtitle = "BBS",
       title = "Long-term")
t2a <- tt_map_list[[1]][["TY1970-1980"]] +
  labs(title = "First ten years")
t2b <- tt_map_list[[1]][["TY2009-2019"]] +
  labs(title = "Last ten years")
t3 <- tt_map_list[[2]][["TY1966-2019"]] +
  labs(subtitle = "CBC")
t4a <- tt_map_list[[2]][["TY1970-1980"]]
t4b <- tt_map_list[[2]][["TY2009-2019"]]
t5 <- tt_map_list[[3]][["TY1980-2019"]] +
  labs(subtitle = "Shorebird")
t6a <- tt_map_list[[3]][["TY1980-1990"]]
t6b <- tt_map_list[[3]][["TY2009-2019"]]

tcomb = t1 + t2a +t2b + t3 + t4a +t4b + t5 + t6a +t6b +
  plot_layout(ncol = 3,byrow = TRUE,
              guides = "collect")

pdf(file = paste0("Figures/Figure_7.pdf"),
    width = 7,
    height = 8)


print(tcomb)
dev.off()









