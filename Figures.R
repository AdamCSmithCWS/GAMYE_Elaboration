#### Figures for publication
library(tidyverse)
library(cmdstanr)
library(posterior)
library(sf)
library(patchwork)
library(geofacet)
library(brms)
library(ggrepel)
source("functions/indices_cmdstan.R")
source("functions/posterior_summary_functions.R")
source("Functions/palettes.R")
species = "Yellow-headed Blackbird"  
species_f <- gsub(species,pattern = " ",replacement = "_")

load(paste0("maps/Simulated",species_f,"_route_maps_data.RData"))


# 1 map of example strata connections ---------------------------------


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


# 1alt spline basis plot --------------------------------------------------

# plot the basis functions ------------------------------------------------
tp <- "non_linear"

load(paste0("Data/Simulated_data_",species_f,"_",tp,"_BBS.RData"))

basis_plot <- data.frame(GAM_year$Year_basis) %>% 
  mutate(n = rep(min_year:2019)) %>% 
  pivot_longer(.,
               cols = starts_with("X"),
               names_to = "knot",
               values_to = "basis") %>% 
  mutate(K = as.integer(str_remove(knot,pattern = "X"))) %>% 
  filter(K < 9)
lbl <- basis_plot %>% 
  filter(n == 2019)
lbl2 <- basis_plot %>% 
  filter(n == 1966)

bpl <- ggplot(data = basis_plot,aes(x = n,y = basis,group = K))+
  geom_line(aes(colour = K))+
  scale_x_continuous(limits = c(1960,2030))+
  xlab("")+
  ylab("Basis function")+
  geom_label_repel(data = lbl,
                   aes(x = n,y = basis,label = K,colour = K),
                   min.segment.length = 0,
                   nudge_x = 4,
                   seed = 2019)+
  geom_label_repel(data = lbl2,
                   aes(x = n,y = basis,label = K,colour = K),
                   min.segment.length = 0,
                   nudge_x = -4,
                   seed = 2019)+
  scale_color_viridis_c(option = "D")+
  theme(legend.position = "none")+
  theme_classic()

print(bpl)


dat = data.frame(y = rep(min_year:2019),
                 x = 1)

Md = mgcv::smoothCon(s(y,k = GAM_year$nknots_Year-2, bs = "tp"),data = dat,
                    absorb.cons=TRUE,#this drops the constant, so that identifibility constraints are abosrbed into the basis
                    diagonal.penalty=TRUE) ## If TRUE then the smooth is reparameterized to turn the penalty into an identity matrix, with the final diagonal elements zeroed (corresponding to the penalty nullspace). 

identical(Md[[1]]$X,GAM_year$Year_basis)

basis_plot2 <- data.frame(Md[[1]]$X) %>% 
  mutate(n = rep(min_year:2019)) %>% 
  pivot_longer(.,
               cols = starts_with("X"),
               names_to = "knot",
               values_to = "basis") %>% 
  mutate(K = as.integer(str_remove(knot,pattern = "X")))


lbl2 <- basis_plot2 %>% 
  filter(n == 2019)
bpl2 <- ggplot(data = basis_plot2,aes(x = n,y = basis,group = K))+
  geom_line(aes(colour = K))+
  scale_x_continuous(limits = c(1965,2030))+
  xlab("")+
  ylab("Basis function")+
  geom_label_repel(data = lbl2,
             aes(x = n,y = basis,label = K,colour = K),
             min.segment.length = 0,
             nudge_x = 4,
             seed = 2019)+
  scale_color_viridis_c()+
  theme(legend.position = "none")+
  theme_classic()

print(bpl2)


# 2 geo true trajectories and counts --------------------------------------

## show all of the true trajectories and mean counts
## in a geo facet so the reader can see how the counts and trajectories
## vary - geofacet of the mean expected counts from balanced
## and observed counts from realised.


# 3 BETA accuracy -----------------------------------------------------

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





# 3 trajectory geofacet ---------------------------------------------------


output_dir <- "output/"
#tp = "non_linear"
load(paste0("Data/Simulated_data_",species_f,"_",tp,"_BBS.RData"))
load(paste0("Data/",species_f,"BBS","_data.RData"))


  strat_grid <- geofacet::grid_auto(realized_strata_map,
                                    codes = "Stratum_Factored",
                                    seed = 2)

  
  fig4_geo = ggplot(data = nsmooth_comp2,aes(y = True_nsmooth,
                                         x = Year))+
    geom_ribbon(aes(ymin = lci,ymax = uci,
                    fill = version),alpha = 0.3)+
    geom_point(data = nsmooth_comp2,aes(x = Year,y = mean_count),
               alpha = 0.1,
               size = 0.2,
               inherit.aes = FALSE)+
    geom_line(aes(colour = version))+
    scale_colour_viridis_d(aesthetics = c("colour","fill"),
                           begin = 0,
                           end = 0.5,
                           direction = -1)+
    scale_y_continuous(limits = c(0,NA))+
    geofacet::facet_geo(~Stratum_Factored,grid = strat_grid,
                        scales = "free")+
    xlab("")+
    ylab("Mean annual smooth trajectory")+
    theme_bw() +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.text.x = element_text(size = 5))
  
  pdf(file = paste0("Figures/Figure_3.pdf"),
      width = 7,
      height = 10)
  print(fig4_geo)
  dev.off()


# 4 trajectory Masked data --------------------------------------------------

  #tp <- "non_linear"
  
  load(paste0("Data/Simulated_data_",species_f,"_",tp,"_BBS.RData"))
  

  sns <- ""
  mk <- ""
  output_dir <- "output/"
  out_base <- paste0(species_f,"_sim_",sns,mk,tp,"_BBS")
  out_base_sim <- paste0("_sim_",sns,mk,tp,"_BBS")
  
  load(paste0("data/",out_base_sim,"_accuracy_comp.RData"))
  
  smooth_spatial <- nsmooth_comp2 %>% 
    filter(version == "Estimated") %>% 
    mutate(version = "Spatial Full",
           mean_count = NA,
           lci = NA,
           uci = NA)
  
  sns <- ""
  mk <- "mask_"
  output_dir <- "output/"
  out_base <- paste0(species_f,"_sim_",sns,mk,tp,"_BBS")
  out_base_sim <- paste0("_sim_",sns,mk,tp,"_BBS")
  
  load(paste0("data/",out_base_sim,"_accuracy_comp.RData"))
  


  # strat_grid <- geofacet::grid_auto(realized_strata_map,
  #                                   codes = "Stratum_Factored",
  #                                   seed = 3)
  # 
  
  smooth_spatial_mask <- nsmooth_comp2 %>% 
    filter(masked == TRUE,
           version == "Estimated Masked") %>% 
    mutate(version = "Spatial")
  
  smooth_spatial <- smooth_spatial %>% 
    filter(Stratum_Factored %in% unique(smooth_spatial_mask$Stratum_Factored))
  
  
  smooth_t_mask <- nsmooth_comp2 %>% 
    filter(masked == TRUE,
           version == "True") %>% 
    mutate(version = "True",
           mean_count = NA)
  
  
  sns <- "nonSpatial_alt_"
  mk <- "mask_"
  output_dir <- "output/"
  out_base <- paste0(species_f,"_sim_",sns,mk,tp,"_BBS")
  out_base_sim <- paste0("_sim_",sns,mk,tp,"_BBS")
  
  load(paste0("data/",out_base_sim,"_accuracy_comp.RData"))
  
  
  smooth_nonSp_mask <- nsmooth_comp2 %>% 
    filter(masked == TRUE,
           version == "Estimated Masked") %>% 
    mutate(version = "NonSpatial",
           mean_count = NA)
  
  #smooth_spatial,
  
  nsmooth_comp_mask <- bind_rows(smooth_nonSp_mask,
                                 smooth_spatial_mask,
                                 smooth_t_mask)
  
  fig4_comb = ggplot(data = nsmooth_comp_mask,aes(y = True_nsmooth,
                                             x = Year))+
    geom_ribbon(aes(ymin = lci,ymax = uci,
                    fill = version),alpha = 0.3)+
    geom_point(data = nsmooth_comp_mask,aes(x = Year,y = mean_count),
               alpha = 0.1,
               size = 0.2,
               inherit.aes = FALSE)+
    geom_line(aes(colour = version))+
    scale_colour_viridis_d(aesthetics = c("colour","fill"),
                           begin = 0,
                           end = 1,
                           direction = -1)+
    scale_y_continuous(limits = c(0,NA))+
    facet_wrap(~Stratum_Factored,
                        scales = "free")+
    xlab("")+
    ylab("Mean annual smooth trajectory")+
    theme_bw() +
    theme(legend.position = "bottom",
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.text.x = element_text(size = 5))
  
  
 
  pdf(file = paste0("Figures/Figure_4.pdf"),
      width = 7,
      height = 10)
  print(fig4_comb)
    dev.off()

# 5 Trend comparisons ---------------------------------------------------
    MAs <- round(log(c(0.1,0.5,1,50)),2)# true mean abundances for different simulations
    
output_dir <- "output/"
tp = "non_linear"
load(paste0("Data/Simulated_data_",species_f,"_",tp,"_BBS.RData"))

strat_df <- as.data.frame(strata_mask) %>% 
  select(Stratum_Factored,masked)
sw_trends <- NULL
strat_trends <- NULL

for(sns in c("","nonSpatial_alt_")){#,"nonSpatial_"))
  for(ma in MAs){
    
    out_base_sim <- paste0("sim_",sns,tp,"_",ma,"_BBS")
    
    lbl <- "Spatial"
    if(sns == "nonSpatial_alt_"){
      lbl <- "NonSpatial"
    }
    lblm <- paste0(signif(exp(ma),1))
      
   
    
    load(paste0("data/",out_base_sim,"_accuracy_comp.RData"))
    sw_t <- all_trends %>% 
      filter(Region_type == "Survey_Wide_Mean",
             #last_year-first_year == 2019,
             last_year-first_year == 5) %>% 
      mutate(version = lbl,
             true_mean = lblm)
    sw_trends <- bind_rows(sw_trends,sw_t)
    
    strat_t <- all_trends %>% 
      filter(Region_type == "Stratum_Factored",
             #last_year-first_year == 2019,
             last_year-first_year == 5) %>% 
      mutate(version = lbl,
             true_mean = lblm)
    strat_trends <- bind_rows(strat_trends,strat_t)
  }
}



strat_trends_nm <- strat_trends %>% 
  left_join(.,strat_df,by = "Stratum_Factored") %>%
  #filter(first_year %in% c(1970:2009)) %>% 
  mutate(t_dif = true_trend - trend,
         t_abs_dif = abs(t_dif),
         t_se = ((uci-lci)/(1.96*2)),
         t_prec = 1/t_se^2)%>% 
  mutate(trend_time = factor(paste0(first_year,"-",last_year)))



mean_difs <- strat_trends_nm %>% 
  group_by(true_mean,
           version,
           trend_time) %>% 
  summarise(mean_abs_dif = mean(t_abs_dif,na.rm = T),
            lci = quantile(t_abs_dif,0.05,na.rm = T),
            uci = quantile(t_abs_dif,0.95,na.rm = T))

sw_trends_ab <- sw_trends %>% 
  mutate(t_dif = true_trend - trend,
         t_abs_dif = abs(t_dif),
         t_se = ((uci-lci)/(1.96*2)),
         t_prec = 1/t_se^2)%>% 
  mutate(trend_time = factor(paste0(first_year,"-",last_year)))

mean_difs_sw <- sw_trends_ab %>% 
  group_by(true_mean,
           version,
           trend_time) %>% 
  summarise(mean_abs_dif = mean(t_abs_dif,na.rm = T),
            lci = quantile(t_abs_dif,0.05,na.rm = T),
            uci = quantile(t_abs_dif,0.95,na.rm = T))

trends4means_plot <- ggplot(data = mean_difs,
                              aes(x = trend_time,
                                  y = mean_abs_dif,
                                  colour = version),
                              position = position_dodge(width = 0.5))+
  geom_errorbar(aes(colour = version,ymin = lci,
                    ymax = uci),
                alpha = 0.3,
                width = 0,
                position = position_dodge(width = 0.5))+
  geom_point(position = position_dodge(width = 0.5))+
  scale_y_continuous(limits = c(0,NA))+
  ylab("Absolute difference in trend (Estimated - True)")+
  xlab("Timespan of Trend")+
  theme_bw()+
  facet_wrap(vars(true_mean),
             nrow = 2,
             scales = "fixed")
print(trends4means_plot)

m1 = brm(t_abs_dif  ~ version*true_mean+trend_time + (1|Stratum_Factored),
        data = strat_trends_nm)
summary(m1)

m2 = brm(t_abs_dif  ~ version*true_mean+trend_time,
         data = sw_trends_ab)
summary(m2)


# trends4a_plot <- ggplot(data = strat_trends_nm,
#                        aes(x = first_year,y = t_abs_dif,fill = version))+
#   geom_boxplot(aes(fill = version),
#              position = position_dodge(width = 1))+
#   scale_y_continuous(limits = c(0,NA))+
#   ylab("Absolute difference in trend (Estimated - True)")+
#   xlab("Timespan of Trend")+
#   theme_bw()+
#   facet_wrap(vars(true_mean),
#              nrow = 2,
#              scales = "free")
# 
#   print(trends4a_plot)
# 
# 
#   
  # strat_trends_m <- strat_trends %>%
  #   filter(first_year %in% c(1970:2009),
  #          masked == TRUE,
  #          version %in% c("NonSpatial Masked","Spatial Masked")) %>%
  #   mutate(first_year = factor(paste0(first_year,"-2019")))

#   m2 = lm(t_abs_dif~version*true_mean+first_year,
#           data = strat_trends_m)
#   summary(m2)
#   
#   
#   trends4b_plot <- ggplot(data = strat_trends_m,
#                           aes(x = first_year,y = t_abs_dif,fill = version))+
#     geom_boxplot(aes(fill = version),
#                  position = position_dodge(width = 1),
#                  coef = 3)+
#     scale_y_continuous(limits = c(0,NA))+
#     ylab("Absolute difference in trend (Estimated - True)")+
#     xlab("Timespan of Trend")+
#     theme_bw()+
#     facet_wrap(vars(true_mean),
#                nrow = 2,
#                scales = "free")
#   
#   print(trends4b_plot)
 
  mean_difs_m <- strat_trends_m %>% 
    group_by(true_mean,
             version,
             first_year) %>% 
    summarise(mean_abs_dif = mean(t_abs_dif,na.rm = T),
              lci = quantile(t_abs_dif,0.05,na.rm = T),
              uci = quantile(t_abs_dif,0.95,na.rm = T))
  
  trends4means_plot_m <- ggplot(data = mean_difs_m,
                          aes(x = first_year,
                              y = mean_abs_dif,
                              colour = version),
                          position = position_dodge(width = 0.5))+
    geom_errorbar(aes(colour = version,ymin = lci,
                      ymax = uci),
                  alpha = 0.3,
                  width = 0,
                  position = position_dodge(width = 0.5))+
    geom_point(position = position_dodge(width = 0.5))+
    scale_y_continuous(limits = c(0,NA))+
    ylab("Absolute difference in trend (Estimated - True)")+
    xlab("Timespan of Trend")+
    theme_bw()+
    facet_wrap(vars(true_mean),
               nrow = 2,
               scales = "free")
  
  pdf(file = paste0("Figures/Figure_5.pdf"),
      width = 7,
      height = 4)
  print(trends4means_plot_m)
  dev.off()
  
  

# 5a Smooth accuracry -----------------------------------------------------

  
  
  MAs <- round(log(c(0.1,0.5,1,5,10,50)),2)# true mean abundances for different simulations
  
  output_dir <- "output/"
  tp = "non_linear"
  load(paste0("Data/Simulated_data_",MAs[1],"_",tp,"_BBS.RData"))
  strat_df <- as.data.frame(strata_mask) %>% 
    select(Stratum_Factored,masked)
  strat_smooths <- NULL
  
  for(sns in c("","nonSpatial_alt_")){#,"nonSpatial_"))
    for(ma in MAs){
      
      out_base_sim <- paste0("sim_",sns,tp,"_",ma,"_BBS")
      
      lbl <- "Spatial"
      if(sns == "nonSpatial_alt_"){
        lbl <- "NonSpatial"
      }
      lblm <- paste0(signif(exp(ma),1))
      load(paste0("data/",out_base_sim,"_accuracy_comp.RData"))
     
      smooth_t <- smooth_comp %>% 
        mutate(version = lbl,
               true_mean = lblm)
      strat_smooths <- bind_rows(strat_smooths,smooth_t)
    }
  }
  
  mean_smooth_comp <- strat_smooths %>% 
    group_by(version,Stratum_Factored,true_mean) %>% 
    summarise(mean_err = mean(smooth_error),
              mean_abs_err = mean(smooth_abs_error),
              lci_abs_err = quantile(smooth_abs_error,0.25),
              uci_abs_err = quantile(smooth_abs_error,0.75)) %>% 
    mutate(stratf = factor(Stratum_Factored))

  
  # m1 = brm(mean_abs_err  ~ version*stratf+ (1|Year),
  #          data = mean_smooth_comp)
  # summary(m1)

  plot_smooth_comp <- ggplot(data = mean_smooth_comp,
                             aes(x = true_mean,y = mean_abs_err,
                                 groups = version,colour = version))+
    geom_point(position = position_dodge(width = 0.2))+
    geom_errorbar(aes(ymin = lci_abs_err,ymax = uci_abs_err),
                  width = 0,
                  alpha = 0.2,position = position_dodge(width = 0.2))+
    facet_wrap(vars(Stratum_Factored),
               scales = "fixed")
 
print(plot_smooth_comp)
  
  # 6 real data overall trajectories for 3 species --------------------------

  
  load("output/real_data_summaries.RData")
  
  
  fls <- data.frame(species_f = c("Yellow-headed_Blackbird",
                                  "Cinclus_mexicanus",
                                  "Red_Knot"),
                    species = c("Yellow-headed Blackbird",
                                "Cinclus_mexicanus",
                                "Red Knot"),
                    species_l = c("BBS - Yellow-headed Blackbird",
                                  "CBC - American Dipper",
                                  "Shorebird - Red Knot"))
  
  Indices_all_out <- Indices_all_out %>% 
    left_join(.,fls,by = "species")
Iplot <- ggplot(data = Indices_all_out,
                aes(x = true_year,y = median))+
  geom_ribbon(aes(ymin = lci,ymax = uci,fill = version),
              alpha = 0.2)+
  geom_line(aes(colour = version))+
  scale_y_continuous(limits = c(0,NA))+
  scale_colour_viridis_d(aesthetics = c("colour","fill"),
                         begin = 0,
                         end = 0.6,
                         direction = 1)+
  ylab("Survey-wide mean annual predictions")+
  xlab("")+
  theme_bw()+
  theme(legend.position = "none")+
  facet_wrap(vars(species_l),
             nrow = 3,
             ncol = 1,
             scales = "free_y")

pdf(file = "Figures/Figure_6.pdf",
    width = 3.5,
    height = 5)
print(Iplot)
dev.off()



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




# 8 Spaghetti plot Red Knot -------------------------------------------------
load("output/real_data_summaries.RData")


fls <- data.frame(species_f = c("Yellow-headed_Blackbird",
                                "Cinclus_mexicanus",
                                "Red_Knot"),
                  species = c("Yellow-headed Blackbird",
                              "Cinclus_mexicanus",
                              "Red Knot"),
                  species_l = c("BBS - Yellow-headed Blackbird",
                                "CBC - American Dipper",
                                "Shorebird - Red Knot"))

Indices_all_out <- Indices_all_out %>% 
  left_join(.,fls,by = "species")

I_RK <- Indices_all_out %>% 
  filter(species == "Red Knot")
Iplot <- ggplot(data = I_RK,
                aes(x = true_year,y = median))+
  geom_ribbon(aes(ymin = lci,ymax = uci,fill = version),
              alpha = 0.2)+
  geom_line(aes(colour = version))+
  scale_y_continuous(limits = c(0,NA))+
  scale_colour_viridis_d(aesthetics = c("colour","fill"),
                         begin = 0,
                         end = 0.6,
                         direction = 1)+
  ylab("Survey-wide mean annual predictions")+
  xlab("")+
  theme_classic()+
  theme(legend.position = "none")


load("output/Red_Knot_SW_indices.RData")
set.seed(4)
sel <- sample(1:max(sw_smooth$samples$.draw),200)
rk_inds <- sw_smooth$samples %>% 
  filter(.draw %in% sel)

start_val = rk_inds %>% 
  filter(true_year == 1980)%>% 
  mutate(start_value = log10(.value)) %>% 
  select(.draw,start_value) 
rk_inds <- rk_inds %>% 
  left_join(.,start_val,by = ".draw")

spagl <- ggplot(data = rk_inds,
                aes(x = true_year,
                    y = .value,
                    group = .draw,
                    colour = start_value))+
  geom_line(alpha = 0.3)+
  xlab("")+
  ylab("log(Survey-wide smooth)")+
  theme_classic()+
  scale_colour_viridis_c(end = 0.9)+
  scale_y_continuous(trans = "log10")+
  theme(legend.position = "none")

pdf(file = paste0("Figures/Figure_8.pdf"),
    width = 7,
    height = 4)

print(Iplot + spagl) 

dev.off()





