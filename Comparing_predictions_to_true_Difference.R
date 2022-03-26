setwd("C:/GitHub/GAMYE_Elaboration")

library(tidyverse)
library(cmdstanr)
library(posterior)
library(sf)
library(geofacet)
source("functions/indices_cmdstan.R")
source("functions/posterior_summary_functions.R")
source("Functions/palettes.R")


MAs <- round(log(c(0.1,0.5,1,5,10,50)),2)# true mean abundances for different simulations


fls <- data.frame(species_f = c("Yellow-headed_Blackbird",
                                "Cinclus_mexicanus",
                                "Red_Knot",
                                "Dickcissel",
                                rep("Yellow-headed_Blackbird",length(MAs)*3),
                                "Yellow-headed_Blackbird",
                                "Cinclus_mexicanus"),
                  species = c("Yellow-headed Blackbird",
                              "Cinclus_mexicanus",
                              "Red Knot",
                              "Dickcissel",
                              rep("Yellow-headed Blackbird",length(MAs)*3),
                              "Yellow-headed Blackbird",
                              "Cinclus_mexicanus"),
                  data = c("BBS",
                           "CBC",
                           "Shorebird",
                           "BBS",
                           rep("BBS",length(MAs)*3),
                           "BBS",
                           "CBC"),
                  out_base = c(paste0("Yellow-headed_Blackbird","_real_","BBS"),
                               paste0("Cinclus_mexicanus","_CBC_B"),
                               paste0("Red_Knot","_Shorebird"),
                               paste0("Dickcissel","_real_","BBS"),
                               paste0("sim_Spatial_Differencebreakpoint_cycle_",MAs,"_BBS"),
                               paste0("sim_Non_Spatial_Differencebreakpoint_cycle_",MAs,"_BBS"),
                               paste0("sim_Non_Hierarchical_Differencebreakpoint_cycle_",MAs,"_BBS"),
                               paste0("Yellow-headed_Blackbird","_real_Non_spatial","BBS"),
                               paste0("Cinclus_mexicanus","_CBC_Non_spatial")),
                  y1 = c(1966,
                         1966,
                         1980,
                         1966,
                         rep(1966,length(MAs)*3),
                         1966,
                         1966),
                  strat_map_name = c("Stratum_Factored",
                                     "strata_vec",
                                     "stratn",
                                     "Stratum_Factored",
                                     rep("Stratum_Factored",length(MAs)*3),
                                     "Stratum_Factored",
                                     "strata_vec"),
                  real = c(TRUE,TRUE,TRUE,TRUE,
                           rep(FALSE,length(MAs)*3),
                           TRUE,TRUE))



output_dir <- "output/"
load("output/Difference_real_data_summaries.RData")
# 
# 
# 
# # Compare spatial and nonspatial for BBS and CBC --------------------------
# 
# for(species in c("Yellow-headed Blackbird",
#                  "Cinclus_mexicanus")){
#   dd <- ifelse(species == "Yellow-headed Blackbird","BBS","CBC")
#   
#   species_f <- gsub(pattern = " ",replacement = "_",
#                     species)
#   
#   load(paste0("Data/",species_f,dd,"_data.RData"))
#   ma <- -5
#   mean_ab <- signif(exp(ma),2)
#   
#   
#   if(dd == "BBS"){
#     ind_sel <- indices_all_out %>% 
#       filter(species == species,
#              data == dd,
#              simulated_data == "real",
#              version == "smooth")
#     
#     ind_sel_f <- indices_all_out %>% 
#       filter(species == species,
#              data == dd,
#              simulated_data == "real",
#              version == "full")
#     
#     
#     mean_obs <- ind_sel_f %>% 
#       filter(mean_abundance == mean_ab,
#              version == "full") %>% 
#       select(true_year,mean_obs,zero,
#              n_surveys,Stratum,Stratum_Factored,
#              Year) %>% 
#       distinct()
#     
# 
#   }else{
#     realized_strata_map <- realized_strata_map %>% 
#       mutate(Stratum = strat)
#     
#     ind_sel <- indices_all_out %>% 
#       filter(species == species,
#              data == dd,
#              simulated_data == "real",
#              version == "smooth")%>% 
#       mutate(Stratum = strat)
#     
#     ind_sel_f <- indices_all_out %>% 
#       filter(species == species,
#              data == dd,
#              simulated_data == "real",
#              version == "full")%>% 
#       mutate(Stratum = strat)
#     
#     
#     mean_obs <- ind_sel_f %>% 
#       filter(mean_abundance == mean_ab,
#              version == "full") %>% 
#       select(true_year,mean_obs,zero,
#              n_surveys,Stratum,Stratum_Factored,
#              Year) %>% 
#       distinct()
#     
#   }
#   
#   strat_grid <- geofacet::grid_auto(realized_strata_map,
#                                     #codes = "Stratum_Factored",
#                                     names = "Stratum",
#                                     seed = 2019)
#   
#  
#   
#   g_inds <- suppressMessages(ggplot(data = ind_sel,aes(x = true_year,y = median))+
#                                geom_point(data = mean_obs,
#                                           aes(x = true_year,y = mean_obs*zero,
#                                               alpha = n_surveys),
#                                           size = 0.3,
#                                           inherit.aes = FALSE)+
#                                geom_ribbon(aes(ymin = lci,ymax = uci,fill = model),
#                                            alpha = 0.3)+
#                                geom_line(aes(colour = model))+
#                                scale_colour_viridis_d(end = 0.85,begin = 0.2,
#                                                       aesthetics = c("fill","colour"))+
#                                #labs(title = paste("Simulated data mean Abundance",mean_ab))+
#                                scale_y_continuous(limits = c(0,NA))+
#                                geofacet::facet_geo(~Stratum,grid = strat_grid,
#                                                    scales = "free_y")+
#                                theme_bw()+
#                                xlab("Mean annual abundance")+
#                                ylab("")+
#                                theme(strip.text = element_text(size = 6),
#                                      strip.background = element_blank(),
#                                      panel.spacing = unit(0.1,"mm"),
#                                      axis.text.x = element_text(size = 5))) 
#   
#   #print(g_inds)
#   
#   pdf(paste0("Figures/Geofacet_",species_f,"_spatial_vs_non.pdf"),
#       width = 8.5,
#       height = 11)
#   print(g_inds)
#   dev.off()
#   
#   
#   
#   
#   
#   
# }
# 


# Explore predicted vs true trajectories and trends for simulations -----------------------

for(ma in MAs){
mean_ab <- signif(exp(ma),2)
load(paste0("Data/Simulated_data_",ma,"_breakpoint_cycle_BBS.RData"))
realized_strata_map = strata_map

strat_grid <- geofacet::grid_auto(realized_strata_map,
                                  codes = "Stratum_Factored",
                                  names = "Stratum",
                                  seed = 2019)

ind_sel <- indices_all_out %>% 
  filter(!is.na(True_scaled_smooth),
         mean_abundance == mean_ab,
         version == "smooth")
ind_sel$model <- c(rep("Spatial",nrow(ind_sel)/2),
                   rep("Non-Spatial",nrow(ind_sel)/2))

mean_obs <- indices_all_out %>% 
  filter(!is.na(True_scaled_smooth),
         mean_abundance == mean_ab,
         version == "full") %>% 
  select(true_year,mean_obs,zero,
         n_surveys,Stratum,Stratum_Factored,
         Year) %>% 
  distinct()


g_inds <- suppressMessages(ggplot(data = ind_sel,aes(x = true_year,y = median))+
                             geom_line(aes(x = true_year,
                                           y = True_scaled_smooth),
                                       colour = "black",
                                       alpha = 0.9,
                                       size = 1,
                                       inherit.aes = FALSE)+
                             geom_point(data = mean_obs,
                                        aes(x = true_year,y = mean_obs*zero,
                                            alpha = n_surveys),
                                        size = 0.3,
                                        inherit.aes = FALSE)+
                             geom_ribbon(aes(ymin = lci,ymax = uci,fill = model),
                                         alpha = 0.3)+
                             geom_line(aes(colour = model))+
                             scale_colour_viridis_d(end = 0.85,begin = 0.2,
                                                    aesthetics = c("fill","colour"))+
                             #labs(title = paste("Simulated data mean Abundance",mean_ab))+
                             scale_y_continuous(limits = c(0,NA))+
                             geofacet::facet_geo(~Stratum,grid = strat_grid,
                                                 scales = "free_y")+
                             theme_bw()+
                             xlab("Mean annual abundance")+
                             ylab("")+
                             theme(strip.text = element_text(size = 6),
                                   strip.background = element_blank(),
                                   panel.spacing = unit(0.1,"mm"),
                                   axis.text.x = element_text(size = 5))) 
  
  #print(g_inds)

pdf(paste0("Figures/Geofacet",mean_ab,"spatial_vs_non.pdf"),
    width = 8.5,
    height = 11)
print(g_inds)
dev.off()

}


# Compare trend estimates -------------------------------------------------

trend_comparison <- NULL

for(ma in MAs){
  mean_ab <- signif(exp(ma),2)
  load(paste0("Data/Simulated_data_",ma,"_breakpoint_cycle_BBS.RData"))
  realized_strata_map = strata_map
  

  t_sel <- stratum_trends %>% 
    filter(mean_abundance == mean_ab,
           data == "BBS") %>% 
    mutate(nyears = last_year - first_year)

  
  t_sel <- t_sel %>% 
    select(Stratum_Factored,first_year,last_year,
           nyears, model,
           trend,mean_abundance) %>%
    distinct() %>% 
    pivot_wider(.,names_from = model,
                values_from = trend)
    
  true_smooth <- log_true_traj %>% 
    select(Stratum,Stratum_Factored,
           Year, True_scaled_smooth)
  tr_yrs <- t_sel %>% 
    select(first_year,last_year) %>% 
    distinct()
  
  true_trends <- NULL
  for(j in 1:nrow(tr_yrs)){
    y1 <- as.numeric(tr_yrs[j,1])
    y2 <- as.numeric(tr_yrs[j,2])
    ny = y2-y1
    
    tmpt <- true_smooth %>% 
      filter(Year %in% c(y1,y2)) %>% 
      pivot_wider(names_from = Year,
                  values_from = True_scaled_smooth,
                  names_prefix = "Y") %>% 
      rename_with(., ~gsub(replacement = "start",
                           pattern = paste0("Y",y1),.x,
                           fixed = TRUE))%>% 
      rename_with(., ~gsub(replacement = "end",
                           pattern = paste0("Y",y2),.x,
                           fixed = TRUE)) %>% 
      mutate(true_trend = texp(end/start,ny),
             first_year = y1,
             last_year = y2) %>% 
      select(Stratum,Stratum_Factored,first_year,last_year,true_trend)
    
    true_trends <- bind_rows(true_trends,tmpt)
    
  }
  
  t_sel <- left_join(t_sel,true_trends,
                     by = c("Stratum_Factored",
                            "first_year",
                            "last_year"))
  
  t_sel <- t_sel %>% 
    mutate(Spatial_error = abs(Spatial - true_trend),
           NonSpatial_error = abs(`Non-Spatial` - true_trend),
           NonHierarchical_error = abs(`Non-Hierarchical` - true_trend),
           dif_error = Spatial_error - NonSpatial_error,
           dif_error_S_H = Spatial_error - NonHierarchical_error)
  
  trend_comparison <- bind_rows(trend_comparison,t_sel)
}


  mean_dif <- trend_comparison %>% 
    group_by(nyears,mean_abundance) %>% 
    summarise(mean_dif = mean(dif_error_S_H),
              sd_dif = sd(dif_error_S_H),
              lq = quantile(dif_error_S_H,0.05),
              uq = quantile(dif_error_S_H,0.95),
              lci = mean_dif-(1.96*sd_dif),
              uci = mean_dif+(1.96*sd_dif))
  
  mns <- ggplot(data = mean_dif,aes(x = nyears,y = mean_dif))+
    geom_errorbar(aes(ymin = lq,ymax = uq),
                  alpha = 0.3,
                  width = 0)+
    geom_point()+
    facet_wrap(vars(mean_abundance),
               nrow = 2,
               ncol = 3,
               scales = "free_y")
  
  print(mns)
  trend_comparison$nyearsF <- factor(trend_comparison$nyears)
  box <- ggplot(data = trend_comparison,
                aes(y = dif_error,x = nyearsF))+
    geom_boxplot(alpha = 0.5)+
    facet_wrap(vars(mean_abundance),
               nrow = 2,
               ncol = 3,
               scales = "free_y")+
    ylab("Difference absolute error in trends Spatial - Non-spatial")+
    xlab("Length of trend (number of years)")+
    geom_abline(slope = 0,intercept = 0,colour = "blue")+
    theme_bw()
  pdf("Figures/Difference_trend_error_Spatial_NonHier.pdf",
      width = 7,
      height = 8)
  print(box)
  dev.off()
  
  
  # Explore Convergence -------------------------------------------------------


failed_rhat <- conv_summaries %>% 
  filter(rhat > 1.05)

failed_ess_bulk <- conv_summaries %>% 
  filter(ess_bulk < 100)


conv_summaries <- conv_summaries %>% 
  filter(species != "Dickcissel") %>% 
  mutate(Simulated_data = ifelse(grepl(x = model, pattern = "^sim_"),
                                 "Simulated","Real"),
         Spatial_model = ifelse(grepl(x = model, pattern = "patial"),
                                "Non-spatial","Spatial"))

x_numbs <- function(x){
  as.numeric(gsub("_","",gsub("[[:alpha:]]", "", x)))
  
}
conv_summaries_sim <- conv_summaries %>% 
  filter(Simulated_data == "Simulated") %>% 
  mutate(Mean_abundance = signif(exp(x_numbs(model)),2))


conv_summaries_real <- conv_summaries %>% 
  filter(Simulated_data == "Real")

sdbetas_sim <- conv_summaries_sim %>% 
  filter(grepl("sdbeta",x = variable)) %>% 
  arrange(ess_bulk) 

sdb_ess_plot = ggplot(data = sdbetas_sim,aes(x = ess_bulk))+
  geom_histogram(bins = 30)+
  facet_wrap(vars(Spatial_model,Mean_abundance),
             nrow = 2,
             ncol = 6,
             scales = "free_y")+
  coord_cartesian(xlim = c(0,6000))+
  geom_vline(xintercept = 500)+
  theme_bw()

pdf("Figures/sd_beta_bulk_ess_simulated.pdf",
    width = 8,
    height = 5)
print(sdb_ess_plot)
dev.off()




# simulated BETA ----------------------------------------------------------



BETAs <- conv_summaries_sim %>% 
  filter(grepl("^BETA\\[",x = variable)) %>% 
  arrange(ess_bulk)

BETAs_plot = ggplot(data = BETAs,aes(x = ess_bulk))+
  geom_histogram(bins = 30)+
  facet_wrap(vars(Spatial_model,Mean_abundance),
             nrow = 2,
             ncol = 6,
             scales = "free_y")+
  coord_cartesian(xlim = c(0,NA))+
  geom_vline(xintercept = 500)+
  xlab("Effective Sample Sizes BETA for simulated data by model and mean abundance")+
  theme_bw()

pdf("Figures/BETA_bulk_ess_simulated.pdf",
    width = 8,
    height = 5)
print(BETAs_plot)
dev.off()


# Real BETA ---------------------------------------------------------------



BETAsr <- conv_summaries_real %>% 
  filter(grepl("^BETA\\[",x = variable)) %>% 
  arrange(ess_bulk)

BETAs_plot_r = ggplot(data = BETAsr,aes(x = ess_bulk))+
  geom_histogram(bins = 30)+
  facet_wrap(vars(data,Spatial_model),
             nrow = 3,
             ncol = 2,
             scales = "free_y")+
  coord_cartesian(xlim = c(0,NA))+
  geom_vline(xintercept = 500)+
  xlab("Effective Sample Sizes BETA for real data by model and dataset")+
  theme_bw()

pdf("Figures/BETA_bulk_ess_real.pdf",
    width = 8,
    height = 5)
print(BETAs_plot_r)
dev.off()



# Simulated betas ---------------------------------------------------------


betas <- conv_summaries_sim %>% 
  filter(grepl("^beta\\[",x = variable)) %>% 
  arrange(ess_bulk)

betas_plot = ggplot(data = betas,aes(x = ess_bulk))+
  geom_histogram(bins = 30)+
  facet_wrap(vars(Spatial_model,Mean_abundance),
             nrow = 2,
             ncol = 6,
             scales = "free_y")+
  coord_cartesian(xlim = c(0,NA))+
  geom_vline(xintercept = 500,
             colour = "darkgreen")+
  xlab("Effective Sample Sizes, beta for simulated data by model and mean abundance")+
  theme_bw()

pdf("Figures/betas_bulk_ess_simulated.pdf",
    width = 8,
    height = 5)
print(betas_plot)
dev.off()


# real data betas ---------------------------------------------------------


betas_r <- conv_summaries_real %>% 
  filter(grepl("^beta\\[",x = variable)) %>% 
  arrange(ess_bulk)

betas_plot_r = ggplot(data = betas_r,aes(x = ess_bulk))+
  geom_histogram(bins = 30)+
  facet_wrap(vars(data,Spatial_model),
             nrow = 3,
             ncol = 2,
             scales = "free_y")+
  coord_cartesian(xlim = c(0,NA))+
  geom_vline(xintercept = 200,
             colour = "blue",
             alpha = 0.2)+
  geom_vline(xintercept = 0)+
  xlab("Effective Sample Sizes, beta for real data by model and data")+
  theme_bw()

pdf("Figures/betas_bulk_ess_real.pdf",
    width = 8,
    height = 5)
print(betas_plot_r)
dev.off()


# real smooths ------------------------------------------------------------



smooths_r <- conv_summaries_real %>% 
  filter(grepl("^nsmooth\\[",x = variable)) %>% 
  arrange(ess_bulk)

smooths_plot_r = ggplot(data = smooths_r,aes(x = ess_bulk))+
  geom_histogram(bins = 30)+
  facet_wrap(vars(data,Spatial_model),
             nrow = 3,
             ncol = 2,
             scales = "free_y")+
  coord_cartesian(xlim = c(0,NA))+
  geom_vline(xintercept = 200,
             colour = "blue",
             alpha = 0.2)+
  geom_vline(xintercept = 0)+
  xlab("Effective Sample Sizes, predicted smooth trajectories for real data by model and data")+
  theme_bw()

pdf("Figures/smooth_trajectories_bulk_ess_real.pdf",
    width = 8,
    height = 5)
print(smooths_plot_r)
dev.off()






# simulated smooths ------------------------------------------------------------



smooths <- conv_summaries_sim %>% 
  filter(grepl("^nsmooth\\[",x = variable)) %>% 
  arrange(ess_bulk)

smooths_plot = ggplot(data = smooths,aes(x = ess_bulk))+
  geom_histogram(bins = 30)+
  facet_wrap(vars(Spatial_model,Mean_abundance),
             nrow = 2,
             ncol = 6,
             scales = "free_y")+
  coord_cartesian(xlim = c(0,NA))+
  geom_vline(xintercept = 200,
             colour = "blue",
             alpha = 0.2)+
  geom_vline(xintercept = 0)+
  xlab("Effective Sample Sizes, predicted smooth trajectories for simulated
       data by model and mean abundance")+
  theme_bw()

pdf("Figures/smooth_trajectories_bulk_ess_simulated.pdf",
    width = 8,
    height = 5)
print(smooths_plot)
dev.off()





trajs <- conv_summaries %>% 
  filter(grepl("^n\\[",x = variable))%>% 
  arrange(ess_bulk)

trajs_ess_plot = ggplot(data = trajs,aes(x = ess_bulk))+
  geom_histogram(bins = 20)+
  facet_wrap(vars(model),
             nrow = 4,
             ncol = 5)+
  geom_vline(xintercept = 500)
print(trajs_ess_plot)
