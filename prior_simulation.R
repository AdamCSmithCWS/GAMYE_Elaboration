### prior simulation of spatial GAM sd

library(bbsBayes)
library(tidyverse)
library(cmdstanr)

setwd("C:/GitHub/GAMYE_Elaboration")


# fit model with fixed data for all parameters except the local sm --------

species <- "Yellow-headed Blackbird"
species_f <- gsub(species,pattern = " ",replacement = "_")


bbs_trends <- read.csv("data_basic/2019All BBS trends stratum.csv")

sd_trends_long <- bbs_trends %>% 
  filter(Trend_Time == "Long-term",
         !grepl(pattern = "^unid",species)) %>% 
  group_by(species) %>% 
  summarise(sd_trends = sd(Trend),
            min_trend = min(Trend),
            max_trend = max(Trend),
            uci_trend = quantile(Trend,0.95),
            lci_trend = quantile(Trend,0.05)) %>% 
  mutate(range_trend = max_trend - min_trend,
         span_90_trend = uci_trend-lci_trend) %>% 
  arrange(-span_90_trend)

hist(sd_trends_long$range_trend)
hist(sd_trends_long$span_90_trend)


# USGS trends -------------------------------------------------------------
# download_bbs_data(sb_id = sb_items[6,2],
#                   bbs_dir = "data/")



bbs_trends_usgs <- read.csv("data/BBS_1966-2019_core_best_trend.csv")

sd_trends_long_usgs <- bbs_trends_usgs %>% 
  filter(Region.Name != "Survey-Wide",
         Region.Name != "CA1",
         Region.Name != "US1",
         !grepl(pattern = "^unid",Species.Name)) %>% 
  group_by(Species.Name) %>% 
  summarise(sd_trends = sd(Trend),
            min_trend = min(Trend),
            max_trend = max(Trend),
            uci_trend = quantile(Trend,0.95),
            lci_trend = quantile(Trend,0.05)) %>% 
  mutate(range_trend = max_trend - min_trend,
         span_90_trend = uci_trend-lci_trend) %>% 
  arrange(-span_90_trend)

hist(sd_trends_long_usgs$range_trend)
hist(sd_trends_long_usgs$span_90_trend)



for(pp in c("gamma","norm","t")){
  for(prior_scale in c(0.5,1,2,3,4)){
    
  tp = paste0(pp,"_rate_",prior_scale)

  #STRATA_True <- log(2)
  output_dir <- "output/"
  out_base <- paste0(species_f,"_sim_",tp,"_BBS")
  csv_files <- paste0(out_base,"-",1:3,".csv")
  
  if(pp == "gamma"){
    pnorm <- 0
  }
  if(pp == "norm"){
    pnorm <- 1
  }
  if(pp == "t"){
    pnorm <- 2
  }
  
  if(!file.exists(paste0(output_dir,csv_files[1]))){
    
    load(paste0("Data/Real_data_",species_f,"_BBS.RData"))
    
    tmp_data = original_data_df
    
    nstrata = max(strata_df$Stratum_Factored)
    nyears = max(tmp_data$Year_Index)
    
    
    N_edges = neighbours$N_edges
    node1 = neighbours$node1
    node2 = neighbours$node2
    
    nknots_year = GAM_year$nknots_Year
    year_basis = GAM_year$Year_basis
    
    stan_data = list(#scalar indicators
      nstrata = nstrata,
      nyears = nyears,
      
      
      #spatial structure
      N_edges = N_edges,
      node1 = node1,
      node2 = node2,
      
      #GAM structure
      nknots_year = nknots_year,
      year_basis = year_basis,
      
      prior_scale = prior_scale,
      pnorm = pnorm
      )
    
    
    
    
    # Fit model ---------------------------------------------------------------
    
    print(paste("beginning",tp,Sys.time()))
    
    mod.file = "models/gamye_iCAR_sd_prior_sim.stan"
    
    ## compile model
    model <- cmdstan_model(mod.file)
    
    
    # Initial Values ----------------------------------------------------------
    
    
    init_def <- function(){ list(sdbeta = runif(nknots_year,0.01,0.1),
                                 beta_raw = matrix(rnorm(nknots_year*nstrata,0,0.01),nrow = nstrata,ncol = nknots_year))}
    
    stanfit <- model$sample(
      data=stan_data,
      refresh=100,
      chains=2, iter_sampling=1000,
      iter_warmup=500,
      parallel_chains = 2,
      #pars = parms,
      adapt_delta = 0.8,
      max_treedepth = 14,
      seed = 123,
      init = init_def,
      output_dir = output_dir,
      output_basename = out_base)
    
    
    #stanfit1 <- as_cmdstan_fit(files = paste0(output_dir,csv_files))
    
    
    save(list = c("stanfit","stan_data","csv_files",
                  "out_base"),
         file = paste0(output_dir,"/",out_base,"_gamye_iCAR.RData"))
    
    
    
  }
  
}#end prior_scale loop
}#end pp loop



# post model summary of priors --------------------------------------------

source("Functions/posterior_summary_functions.R")

nsmooth_out <- NULL
trends_out <- NULL
summ_out <- NULL

for(pp in c("gamma","norm","t")){
  for(prior_scale in c(0.5,1,2,3,4)){
    
    tp = paste0(pp,"_rate_",prior_scale)
    
    #STRATA_True <- log(2)
    output_dir <- "output/"
    out_base <- paste0(species_f,"_sim_",tp,"_BBS")
    csv_files <- paste0(out_base,"-",1:3,".csv")
    

load(paste0(output_dir,"/",out_base,"_gamye_iCAR.RData"))

summ = stanfit$summary()

summ <- summ %>% 
  mutate(prior_scale = prior_scale,
         distribution = pp)


nsmooth_samples <- posterior_samples(stanfit,
                                 parm = "nsmooth",
                                 dims = c("Stratum_Factored","Year_Index"))

nsmooth <- nsmooth_samples %>% 
  posterior_sums(.,
                 dims = c("Stratum_Factored","Year_Index"))%>% 
  mutate(prior_scale = prior_scale,
         distribution = pp)

nyears = max(nsmooth_samples$Year_Index)
nyh <- paste0("Y",nyears)
trs <- function(y1,y2,ny){
  tt <- (((y2/y1)^(1/ny))-1)*100
}
trends <- nsmooth_samples %>% 
  filter(Year_Index %in% c(1,nyears)) %>% 
  select(.draw,.value,Stratum_Factored,Year_Index) %>% 
  pivot_wider(.,names_from = Year_Index,
              values_from = .value,
              names_prefix = "Y") %>%
  rename_with(.,~gsub(pattern = nyh,replacement = "Y2", .x)) %>% 
  group_by(.draw,Stratum_Factored) %>% 
  summarise(trend = trs(Y1,Y2,nyears))%>% 
  mutate(prior_scale = prior_scale,
         distribution = pp)


nsmooth_out <- bind_rows(nsmooth_out,nsmooth)
trends_out <- bind_rows(trends_out,trends)
summ_out <- bind_rows(summ_out,summ)

  }#prior_scale
  print(paste(pp,prior_scale))
}# pp

# save(file = "output/prior_sim_summary.RData",
#      list = c("nsmooth_out",
#               "trends_out",
#               "summ_out"))



# summarise ---------------------------------------------------------------

load("output/prior_sim_summary.RData")

trends_out <- trends_out %>% 
  mutate(abs_trend = abs(trend))


trend_h <- ggplot(data = trends_out)+
  geom_histogram(aes(x = abs_trend),
                 binwidth = 2)+
  coord_cartesian(xlim = c(0,50))+
  facet_wrap(facets = vars(distribution,prior_scale),
             nrow = 3,
             ncol = 5,
             scales = "free")
print(trend_h)


trends_sd <- trends_out %>% 
  group_by(.draw,prior_scale,distribution) %>% 
  summarise(sd_trends = sd(trend),
            min_trend = min(trend),
            max_trend = max(trend),
            lci95 = quantile(trend,0.025),
            uci95 = quantile(trend,0.975),
            span_95 = uci95 - lci95,
            range_trend = max_trend - min_trend)

trendsd_h <- ggplot(data = trends_sd)+
  geom_histogram(aes(x = range_trend),
                 binwidth = 2)+
  coord_cartesian(xlim = c(0,50))+
  facet_wrap(facets = vars(distribution,prior_scale),
             nrow = 3,
             ncol = 5,
             scales = "free")
print(trendsd_h)






# 
# trend_f <- ggplot(data = trends_out)+
#   geom_freqpoly(aes(x = trend, group = prior_scale,
#                     colour = prior_scale),
#                 bins = 1000)+
#   scale_color_viridis_c()+
#   coord_cartesian(xlim = c(-20,20),
#                   ylim = c(0,10000))+
#   #scale_x_continuous(limits = c(-30,30))+
#   facet_wrap(facets = vars(distribution),
#              nrow = 2,
#              scales = "fixed")
# print(trend_f)


# 95%ile of trends ---------------------------------------------------------
w_95 <- function(x){
  d <- diff(quantile(x,c(0.025,0.975)))
}
percentile_trends <- trends_out %>% 
  group_by(distribution,prior_scale,.draw) %>% 
  summarise(p_95 = w_95(trend)) 


perc_f <- ggplot(data = percentile_trends)+
  geom_freqpoly(aes(x = p_95, group = prior_scale,
                    colour = prior_scale),
                bins = 50)+
  scale_color_viridis_c()+
  scale_x_continuous(limits = c(0,50))+
  facet_wrap(facets = vars(distribution),
             nrow = 2,
             scales = "fixed")
print(perc_f)





