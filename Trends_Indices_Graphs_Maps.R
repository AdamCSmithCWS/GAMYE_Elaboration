### trends indices graphs and maps for real data
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
                              rep("Yellow-headed_Blackbird",length(MAs)*2),
                              "Yellow-headed_Blackbird",
                              "Cinclus_mexicanus"),
                  species = c("Yellow-headed Blackbird",
                              "Cinclus_mexicanus",
                              "Red Knot",
                              "Dickcissel",
                              rep("Yellow-headed Blackbird",length(MAs)*2),
                              "Yellow-headed Blackbird",
                              "Cinclus_mexicanus"),
                  data = c("BBS",
                           "CBC",
                           "Shorebird",
                           "BBS",
                           rep("BBS",length(MAs)*2),
                           "BBS",
                           "CBC"),
                  out_base = c(paste0("Yellow-headed_Blackbird","_real_","BBS"),
                               paste0("Cinclus_mexicanus","_CBC_B"),
                               paste0("Red_Knot","_Shorebird"),
                               paste0("Dickcissel","_real_","BBS"),
                               paste0("sim_breakpoint_cycle_",MAs,"_BBS"),
                               paste0("sim_nonSpatial_alt_breakpoint_cycle_",MAs,"_BBS"),
                               paste0("Yellow-headed_Blackbird","_real_Non_spatial","BBS"),
                               paste0("Cinclus_mexicanus","_CBC_Non_spatial")),
                  y1 = c(1966,
                         1966,
                         1980,
                         1966,
                         rep(1966,length(MAs)*2),
                         1966,
                         1966),
                  strat_map_name = c("Stratum_Factored",
                                     "strata_vec",
                                     "stratn",
                                     "Stratum_Factored",
                                     rep("Stratum_Factored",length(MAs)*2),
                                     "Stratum_Factored",
                                     "strata_vec"),
                  real = c(TRUE,TRUE,TRUE,TRUE,
                           rep(FALSE,length(MAs)*2),
                           TRUE,TRUE))



# Species loop ------------------------------------------------------------
output_dir <- "output/"
stratum_trends <- NULL
all_trends <- NULL
indices_all_out <- NULL
Indices_all_out <- NULL

tt_map_list = vector(mode = "list",length = nrow(fls))
names(tt_map_list) <- fls$species_f
tt_map_data_list = tt_map_list
ind_plots_list = tt_map_list
Ind_plots_list = tt_map_list

conv_summaries <- NULL

for(i in c(1,2,17,18)){#c(1:nrow(fls))){

  species = fls[i,"species"]
  species_f <- fls[i,"species_f"]
  dd <- fls[i,"data"]
  out_base <- fls[i,"out_base"]
  year_1 = fls[i,"y1"]
  reald = fls[i,"real"]
  st_n = fls[i,"strat_map_name"]
  if(reald){
    load(paste0("Data/",species_f,dd,"_data.RData"))
    ma <- -5
    modl <- ifelse(i < 3,"Spatial","Non-Spatial")
  }else{
    if(i < 11){
    ma <- MAs[i-4]
  load(paste0("Data/Simulated_data_",ma,"_breakpoint_cycle_BBS.RData"))
  realized_strata_map = strata_map
  modl <- "Spatial"
    }else{
      
      ma <- MAs[i-(4+length(MAs))]
      load(paste0("Data/Simulated_data_",ma,"_breakpoint_cycle_BBS.RData"))
      realized_strata_map = strata_map 
      modl <- "Non-Spatial"
    }
  }
  
  load(paste0(output_dir,"/",out_base,"_gamye_iCAR.RData"))
  csv_files <- paste0(output_dir,out_base,"-",1:3,".csv")
  

  
  
  stanfit <- as_cmdstan_fit(files = csv_files)
 strat_df <- as.data.frame(realized_strata_map)
  

# convergence summary -----------------------------------------------------

  
  stanf_df <- stanfit$draws(format = "df")
  
  #optional removes the count-level parameters
   drws <- stanf_df #%>% 
  #   select(!contains(c("E[","noise_raw"))) 
  
  conv_summ <- summarise_draws(drws) %>% 
    mutate(species = species,
           model = out_base,
           data = dd)
  
  conv_summaries <- bind_rows(conv_summaries,
                              conv_summ)
  
  
# exploring the spatial abundance dist ------------------------------------
# 
#   strat_samples <- posterior_samples(fit = stanfit,
#                                      parm = "strata_raw",
#                                      dims = "stratnumber")
#   STRAT_samples <- posterior_samples(fit = stanfit,
#                                      parm = "STRATA")
#   sdstrata_samples <- posterior_samples(fit = stanfit,
#                                         parm = "sdstrata")
#   sd_tmp = sdstrata_samples %>% select(.draw,.value) %>% 
#     rename(sdstrata = .value)
#   STRAT_tmp = STRAT_samples %>% select(.draw,.value) %>% 
#     rename(STRAT = .value)
#   strata_sums <- strat_samples %>% 
#     left_join(sd_tmp) %>% 
#     left_join(STRAT_tmp) %>% 
#     mutate(strat = .value*sdstrata + STRAT) %>% 
#     group_by(stratnumber) %>% 
#     summarise(mean_est = mean(strat),
#               lci_est = quantile(strat,0.025),
#               uci_est = quantile(strat,0.975))
#   
#   
  stan_df <- data.frame(stratnumber = stan_data$strat,
                          count = stan_data$count,
                          site = stan_data$site,
                          Year = stan_data$year) 
  
  
 
  if((dd %in% c("BBS","CBC")) & reald){
    # obs_strat_means <- stan_df %>% 
    #   group_by(stratnumber) %>% 
    #   summarise(mean = mean(count),
    #             max = max(count),
    #             lmean = log(mean)) %>% 
    #   left_join(.,strata_sums)
    
  strat_z <- data.frame(stratnumber = 1:stan_data$nstrata,
                        zero = stan_data$nonzeroweight)
  }else{
    strat_z <- data.frame(stratnumber = 1:stan_data$nstrata,
                          zero = 1)
  }
  obs_means <- stan_df %>% 
    group_by(stratnumber,Year) %>% 
    summarise(n_surveys = n(),
              mean_obs = mean(count),
              max_obs = max(count)) %>% 
    left_join(.,strat_z) %>% 
    rename_with(.,~gsub(x = .x,pattern = "stratnumber",
                        replacement = st_n))

  

  
 
  
  # plot_obs_m <- realized_strata_map %>% 
  #   rename_with(~gsub(pattern = st_n,replacement = "stratnumber",x = .x)) %>% 
  #   left_join(.,obs_strat_means,by = c("stratnumber"))
  # 
  # 
  # obs_m_map <- ggplot()+
  #   geom_sf(data = plot_obs_m,aes(fill = lmean))+
  #   scale_fill_viridis_c()
  # print(obs_m_map)
  # 
  # est_m_map <- ggplot()+
  #   geom_sf(data = plot_obs_m,aes(fill = mean_est))+
  #   scale_fill_viridis_c()
  # print(est_m_map)
  # 
  # plot(obs_strat_means$lmean,obs_strat_means$mean_est)
  # abline(0,1)

# indices and trends ------------------------------------------------------

  strat_inds <- index_function(fit = stanfit,
                               parameter = "n",
                               year_1 = year_1,
                               strat = st_n)
  strat_inds_smooth <- index_function(fit = stanfit,
                                      parameter = "nsmooth",
                                      year_1 = year_1,
                                      strat = st_n)
  
  if(!reald){
    true_smooths <- log_true_traj %>% 
      ungroup() %>% 
      select(Stratum_Factored,Year,
             True_scaled_smooth) %>%
      rename(true_year = Year) 
    
    indices <- strat_inds$indices %>% 
      left_join(.,strat_df)%>% 
      mutate(version = "full") %>% 
      left_join(.,obs_means) %>% 
      left_join(.,true_smooths,by = c("true_year","Stratum_Factored"))
    
    indices_smooth <- strat_inds_smooth$indices %>% 
      left_join(.,strat_df) %>% 
      mutate(version = "smooth")%>% 
      left_join(.,true_smooths,by = c("true_year","Stratum_Factored"))
    
    
  }else{
    indices <- strat_inds$indices %>% 
      left_join(.,strat_df)%>% 
      mutate(version = "full") %>% 
      left_join(.,obs_means) 
    
    indices_smooth <- strat_inds_smooth$indices %>% 
      left_join(.,strat_df) %>% 
      mutate(version = "smooth")
    
    
  }

 if(reald){rd = "real"}else{rd = "simulated"}

  indices_all <- bind_rows(indices,
                           indices_smooth) %>% 
    rename_with(.,~gsub(st_n,replacement = "strat_plot",
                        .x)) %>% 
    mutate(original_strat_name = strat_plot,
           species = species,
           region_type = "Stratum",
           simulated_data = rd,
           mean_abundance = signif(exp(ma),2),
           data = dd,
           model = modl) %>% 
    rename_with(.,~gsub(replacement = st_n,pattern = "strat_plot",
                        .x))
  indices_all_out <- bind_rows(indices_all_out,indices_all)
  
  pd = ceiling(sqrt(length(unique(indices_all$original_strat_name))))  
  
  
  
  pl_inds <- ggplot(data = indices_all,aes(x = true_year,y = median))+
    geom_point(aes(x = true_year,y = mean_obs*zero,alpha = n_surveys),
               inherit.aes = FALSE)+
    geom_ribbon(aes(ymin = lci,ymax = uci,fill = version),alpha = 0.1)+
    geom_line(aes(colour = version))+
    labs(title = species)+
    scale_y_continuous(limits = c(0,NA))+
    facet_wrap(vars(original_strat_name),
               nrow = pd,
               ncol = pd,
               scales = "free_y")
 
  if(!reald){
    pl_inds <- pl_inds +
      geom_line(data = indices_all,
                aes(x = true_year,
                    y = True_scaled_smooth),
                colour = "black",
                alpha = 0.3,
                inherit.aes = FALSE)
  }
 # print(pl_inds)
# geofacet ----------------------------------------------------------------

  # row_col_cat <- function(x,n){
  #   cat_b = c(0,1/n*(1:n))
  #   qq = quantile(x,cat_b)
  #   qq[1] <- qq[1]-1
  #   qq[n+1] <- qq[n+1]+1
  #   cc <- as.integer(cut(x,
  #                        qq,
  #                        ordered_result = TRUE))
  #   return(cc)
  # }
  # 

  # 
  
  #   centres <- suppressWarnings(st_centroid(realized_strata_map)) %>% 
  #     rename_with(.,~gsub(pattern = st_n,
  #                         replacement = "strat_tmp",
  #                         x = .x)) %>% 
  #     arrange(strat_tmp) %>% 
  #     rename_with(.,~gsub(replacement = st_n,
  #                         pattern = "strat_tmp",
  #                         x = .x))
  #   
  # coords <- st_coordinates(centres)%>%
  #   as.data.frame() 
  # coords[,st_n] <- 1:nrow(coords)
  # coords <- coords %>% 
  #   mutate(y_cat = row_col_cat(-Y,pd)) %>% 
  #   group_by(y_cat) %>% 
  #   mutate(x_cat = rank(X))
  # 
  # inds_geo <- indices_all %>% 
  #   left_join(.,coords)
  #   
  # lbls = inds_geo %>% 
  #   ungroup() %>% 
  #   select(x_cat,y_cat,strat) %>% 
  #   distinct()
  if(dd == "Shorebird"){
  strat_grid <- geofacet::grid_auto(realized_strata_map,
                                    codes = st_n,
                                    names = "hex_name",
                                    seed = 2019)
  }
  if(dd == "BBS"){
  strat_grid <- geofacet::grid_auto(realized_strata_map,
                                    codes = st_n,
                                    names = "Stratum",
                                    seed = 2019)
  }
  if(dd == "CBC"){
    strat_grid <- geofacet::grid_auto(realized_strata_map,
                                      codes = st_n,
                                      names = "stratum",
                                      seed = 2019)
  }
  indices_geo <- indices_all %>% 
    rename_with(~gsub(x = .x,
                      pattern = st_n,
                      replacement = "strat_labs")) 
  

  g_inds <- suppressMessages(ggplot(data = indices_geo,aes(x = true_year,y = median))+
    geom_point(aes(x = true_year,y = mean_obs*zero,alpha = n_surveys),
               inherit.aes = FALSE)+
    geom_ribbon(aes(ymin = lci,ymax = uci,fill = version),alpha = 0.1)+
    geom_line(aes(colour = version))+
    labs(title = species)+
    scale_y_continuous(limits = c(0,NA))+
    geofacet::facet_geo(~strat_labs,grid = strat_grid,
                        scales = "free")+
    theme(strip.text = element_text(size = 6),
          strip.background = element_blank(),
          panel.spacing = unit(0.1,"mm"),
    axis.text.x = element_text(size = 5)))

  if(!reald){
    
    g_inds <- g_inds +
      geom_line(aes(x = true_year,
                    y = True_scaled_smooth),
                colour = "black",
                alpha = 0.3,
                inherit.aes = FALSE)
    
  }
    pdf(paste0("figures/",out_base,"_geofacet.pdf"),
      height = 14,
      width = 14)
  
  print(g_inds)
 dev.off()
 
       
ind_plots_list[[species_f]] <- pl_inds

tyrs = unique(c(2009,1999,1990,1980,1970,year_1))
tyrs = tyrs[which(tyrs >= year_1)]
tyrs2 <- rep(2019,length(tyrs))
tyrs2 <- c(tyrs2,tyrs+10)
tyrs <- c(tyrs,tyrs)


tt_map_data <- vector(mode = "list",length = length(tyrs))
names(tt_map_data) <- paste0("TY",tyrs,"-",tyrs2)
tt_map = tt_map_data

for(j in 1:length(tyrs)){
  yy = tyrs[j]
  yy2 = tyrs2[j]
  
  tt <- trends_function(ind_list = strat_inds_smooth,
                        start_year = yy,
                        end_year = yy2) %>% 
    mutate(species = species,
           first_year = yy,
           last_year = yy2,
           region_type = "Stratum",
           simulated_data = rd,
           mean_abundance = signif(exp(ma),2),
           data = dd,
           model = modl)
  
 ttmd <- realized_strata_map %>% 
    left_join(.,tt,by = st_n) %>% 
   mutate(trend_plot = trend_plot_cats(trend))
 
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
if(reald){
  pdf(file = paste0("Figures/",out_base,"trend_maps.pdf"))
}else{
  pdf(file = paste0("Figures/",out_base,"trend_maps.pdf"))
}
for(j in 1:length(tt_map)){
  print(tt_map[[j]])
}
dev.off()


  # Overall Annual indices

if(dd == "Shorebird"){
  
  sw_inds <- index_function(fit = stanfit,
                               parameter = "N",
                               year_1 = year_1,
                            strat = NULL,
                            first_dim = "y")
  Inds <- sw_inds$indices %>% 
    mutate(version = "full",
           simulated_data = rd,
           mean_abundance = signif(exp(ma),2),
           data = dd,
           model = modl)
  
  sw_smooth <- index_function(fit = stanfit,
                                      parameter = "NSmooth",
                                      year_1 = year_1,
                                      strat = NULL,
                                      first_dim = "y")
 ### one-off save of overall smooth indices for spaghetti plot
   save(list = "sw_smooth",
       file = "output/Red_Knot_SW_indices.RData")
  
  Inds_smooth <- sw_smooth$indices %>% 
    mutate(version = "smooth",
           simulated_data = rd,
           mean_abundance = signif(exp(ma),2),
           data = dd,
           model = modl)
  
  Indices_all <- bind_rows(Inds,
                           Inds_smooth)  %>% 
    mutate(species = species,
           region_type = "Survey_wide")
  
  Indices_all_out <- bind_rows(Indices_all_out,Indices_all)
  
  
  
  pl_Inds <- ggplot(data = Indices_all,aes(x = Year,y = median))+
    geom_ribbon(aes(ymin = lci,ymax = uci,fill = version),alpha = 0.1)+
    geom_line(aes(colour = version))+
    labs(title = species)
  
  ind_plots_list[[species_f]] <- pl_Inds
  
  
}else{
 
  
  sw_inds <- index_function(fit = stanfit,
                               parameter = "n",
                               year_1 = year_1,
                               strat = st_n,
                               weights_df = strat_df,
                               area = "AREA_1",#"Area",
                               summary_regions = NULL)
  Inds <- sw_inds$indices %>% 
    mutate(version = "full",
           simulated_data = rd,
           mean_abundance = signif(exp(ma),2),
           data = dd,
           model = modl)
  
  sw_smooth <- index_function(fit = stanfit,
                                      parameter = "nsmooth",
                                      year_1 = year_1,
                                      strat = st_n,
                                      weights_df = strat_df,
                                      area = "AREA_1",#"Area",
                                      summary_regions = NULL)
  Inds_smooth <- sw_smooth$indices %>% 
    mutate(version = "smooth",
           simulated_data = rd,
           mean_abundance = signif(exp(ma),2),
           data = dd,
           model = modl)
  
  Indices_all <- bind_rows(Inds,
                           Inds_smooth) %>% 
    mutate(species = species,
           region_type = "Survey_wide")
  
  Indices_all_out <- bind_rows(Indices_all_out,Indices_all)
  
  pl_Inds <- ggplot(data = Indices_all,aes(x = Year,y = median))+
    geom_ribbon(aes(ymin = lci,ymax = uci,fill = version),alpha = 0.1)+
    geom_line(aes(colour = version))+
    labs(title = species)
  
  #print(pl_Inds)
  ind_plots_list[[species_f]] <- pl_Inds
  
  
}
if(reald){
pdf(file = paste0("Figures/",species_f,"_Indices.pdf"))
}else{
pdf(file = paste0("Figures/",out_base,"_Indices.pdf"))
}
print(pl_Inds)
print(pl_inds)
dev.off()

tt_map_data_list[[species_f]] <- tt_map_data
tt_map_list[[species_f]] <- tt_map



for(j in 1:length(tyrs)){
  yy = tyrs[j]
  yy2 = tyrs2[j]
  
  TT <- trends_function(ind_list = sw_smooth,
                        start_year = yy,
                        end_year = yy2) %>% 
    mutate(species = species,
           first_year = yy,
           last_year = yy2,
           region_type = "Survey_wide")
  
all_trends <- bind_rows(all_trends,TT)

}


  
  

print(i)
}#end species loop


save(list = c("all_trends",
              "tt_map_data_list",
              "tt_map_list",
              "indices_all_out",
              "stratum_trends",
              "Indices_all_out",
              "ind_plots_list",
              "Ind_plots_list",
              "conv_summaries"),
     file = "output/real_data_summaries.RData")


