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
  labs(title = "BBS",
       subtitle = "Long-term")
t2a <- tt_map_list[[1]][["TY1970-1980"]] +
  labs(subtitle = "First ten years")
t2b <- tt_map_list[[1]][["TY2009-2019"]] +
  labs(subtitle = "Last ten years")
t3 <- tt_map_list[[2]][["TY1966-2019"]] +
  labs(title = "CBC")
t4a <- tt_map_list[[2]][["TY1970-1980"]]
t4b <- tt_map_list[[2]][["TY2009-2019"]]
t5 <- tt_map_list[[3]][["TY1980-2019"]] +
  labs(title = "Shorebird")
t6a <- tt_map_list[[3]][["TY1980-1990"]]
t6b <- tt_map_list[[3]][["TY2009-2019"]]

tcomb = t1 + t2a +t2b + t3 + t4a +t4b + t5 + t6a +t6b +
  plot_layout(ncol = 3,byrow = TRUE,
              guides = "collect")
print(tcomb)




for(j in 1:nrow(fls)){
  species = fls[j,"species"]
  species_f = fls[j,"species_f"]
  y1 = fls[j,"y1"]
  
}






