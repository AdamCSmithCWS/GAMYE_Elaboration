### trends indices graphs and maps for real data

library(tidyverse)
library(cmdstanr)
library(posterior)
library(sf)
source("functions/indices_cmdstan.R")
source("functions/posterior_summary_functions.R")
source("Functions/palettes.R")


fls <- data.frame(species_f = c("Yellow-headed_Blackbird",
                              "Zonotrichia_albicollis",
                              "Red_Knot"),
                  species = c("Yellow-headed Blackbird",
                              "White-throated Sparrow",
                              "Red Knot"),
                  data = c("BBS",
                           "CBC",
                           "Shorebird"),
                  out_base = c(paste0("Yellow-headed_Blackbird","_real_","BBS"),
                               paste0("Zonotrichia_albicollis","_CBC_B"),
                               paste0("Red Knot","_Shorebird")),
                  y1 = c(1966,
                         1966,
                         1980),
                  strat_map_name = c("Stratum_Factored",
                                     "strata_vec",
                                     "stratn"))




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


for(i in c(1,3)){#1:nrow(fls)){

  species = fls[i,"species"]
  species_f <- fls[i,"species_f"]
  dd <- fls[i,"data"]
  out_base <- fls[i,"out_base"]
  year_1 = fls[i,"y1"]
  
  st_n = fls[i,"strat_map_name"]
  
  
  
  load(paste0("Data/",species_f,dd,"_data.RData"))
  load(paste0(output_dir,"/",out_base,"_gamye_iCAR.RData"))
  csv_files <- paste0(output_dir,out_base,"-",1:3,".csv")
  
  #stan_data
  #neighbours
  #realized_strata_map
  #data_1
  
 
  
  
  stanfit <- as_cmdstan_fit(files = csv_files)
# Strata Annual indices file -----------------------------------------------------
  strat_df <- as.data.frame(realized_strata_map)
  
  

# exploring the spatial abundance dist ------------------------------------

  strat_samples <- posterior_samples(fit = stanfit,
                                     parm = "strata_raw",
                                     dims = "stratnumber")
  STRAT_samples <- posterior_samples(fit = stanfit,
                                     parm = "STRATA")
  sdstrata_samples <- posterior_samples(fit = stanfit,
                                        parm = "sdstrata")
  sd_tmp = sdstrata_samples %>% select(.draw,.value) %>% 
    rename(sdstrata = .value)
  STRAT_tmp = STRAT_samples %>% select(.draw,.value) %>% 
    rename(STRAT = .value)
  strata_sums <- strat_samples %>% 
    left_join(sd_tmp) %>% 
    left_join(STRAT_tmp) %>% 
    mutate(strat = .value*sdstrata + STRAT) %>% 
    group_by(stratnumber) %>% 
    summarise(mean_est = mean(strat),
              lci_est = quantile(strat,0.025),
              uci_est = quantile(strat,0.975))
  
  
  stan_df <- data.frame(stratnumber = stan_data$strat,
                          count = stan_data$count,
                          site = stan_data$site,
                          year = stan_data$year) 
  
  
  obs_strat_means <- stan_df %>% 
    group_by(stratnumber) %>% 
    summarise(mean = mean(count),
              max = max(count),
              lmean = log(mean)) %>% 
    left_join(.,strata_sums)
  
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
  indices <- strat_inds$indices %>% 
    inner_join(.,strat_df)%>% 
    mutate(version = "full")

  strat_inds_smooth <- index_function(fit = stanfit,
                               parameter = "nsmooth",
                               year_1 = year_1,
                               strat = st_n)
  indices_smooth <- strat_inds_smooth$indices %>% 
    inner_join(.,strat_df) %>% 
    mutate(version = "smooth")
  
  indices_all <- bind_rows(indices,
                           indices_smooth) %>% 
    rename_with(.,~gsub(st_n,replacement = "strat_plot",
                        .x)) %>% 
    mutate(original_strat_name = strat_plot,
           species = species,
           region_type = "Stratum") %>% 
    rename_with(.,~gsub(replacement = st_n,pattern = "strat_plot",
                        .x))
  indices_all_out <- bind_rows(indices_all_out,indices_all)
  
  pd = ceiling(sqrt(length(unique(indices_all$original_strat_name))))  
  
  
  
  pl_inds <- ggplot(data = indices_all,aes(x = Year,y = median))+
    geom_ribbon(aes(ymin = lci,ymax = uci,fill = version),alpha = 0.1)+
    geom_line(aes(colour = version))+
    labs(title = species)+
    facet_wrap(vars(original_strat_name),
               nrow = pd,
               ncol = pd,
               scales = "free_y")
    
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
           region_type = "Stratum")
  
 ttmd <- realized_strata_map %>% 
    left_join(.,tt,by = st_n) %>% 
   mutate(trend_plot = trend_plot_cats(trend))
 
 tt_map_data[[paste0("TY",yy,"-",yy2)]] <- ttmd
 
  
 ## import one of the bbsBays binned colour scales
 
 ttm  <- ggplot(data = ttmd)+
   geom_sf(aes(fill = trend_plot))+
   theme_bw()+
   scale_colour_manual(values = map_palette_s, 
                       aesthetics = c("fill"),
                       guide = guide_legend(reverse=TRUE),
                       name = paste0(species," Trend\n",yy,"-",yy2))
 

 # print(ttm)
  
  
 tt_map[[paste0("TY",yy,"-",yy2)]] <- ttm
 
  stratum_trends <- bind_rows(stratum_trends,tt)
}
pdf(file = paste0("Figures/",species_f,"trend_maps.pdf"))
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
    mutate(version = "full")
  
  sw_smooth <- index_function(fit = stanfit,
                                      parameter = "NSmooth",
                                      year_1 = year_1,
                                      strat = NULL,
                                      first_dim = "y")
  
  Inds_smooth <- sw_smooth$indices %>% 
    mutate(version = "smooth")
  
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
    mutate(version = "full")
  
  sw_smooth <- index_function(fit = stanfit,
                                      parameter = "nsmooth",
                                      year_1 = year_1,
                                      strat = st_n,
                                      weights_df = strat_df,
                                      area = "AREA_1",#"Area",
                                      summary_regions = NULL)
  Inds_smooth <- sw_smooth$indices %>% 
    mutate(version = "smooth")
  
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






