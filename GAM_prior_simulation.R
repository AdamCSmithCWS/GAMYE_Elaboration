### prior simulation of spatial GAM sd

library(bbsBayes)
library(tidyverse)
library(cmdstanr)

setwd("C:/GitHub/GAMYE_Elaboration")


# fit model with fixed data for all parameters except the local sm --------


for(pp in c("gamma","norm","t")){
  for(prior_scale in c(0.5,1,2,3,4)){
    
  tp = paste0("GAM_",pp,prior_scale,"_rate")

  #STRATA_True <- log(2)
  output_dir <- "output/"
  out_base <- tp
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
    
    
    
    nknots_year = GAM_year$nknots_Year
    year_basis = GAM_year$Year_basis
    nyears = nrow(year_basis)
    
    stan_data = list(#scalar indicators
      nyears = nyears,

      
      #GAM structure
      nknots_year = nknots_year,
      year_basis = year_basis,
      
      prior_scale = prior_scale,
      pnorm = pnorm
      )
    
    
    
    
    # Fit model ---------------------------------------------------------------
    
    print(paste("beginning",tp,Sys.time()))
    
    mod.file = "models/GAM_prior_sim.stan"
    
    ## compile model
    model <- cmdstan_model(mod.file)
    
    
    # Initial Values ----------------------------------------------------------
    
    
    init_def <- function(){ list(sdbeta = runif(1,0.01,0.1),
                                 BETA_raw = rnorm(nknots_year,0,0.01))}
    
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
    
    tp = paste0("GAM_",pp,prior_scale,"_rate")
    
    #STRATA_True <- log(2)
    output_dir <- "output/"
    out_base <- tp
    csv_files <- paste0(out_base,"-",1:3,".csv")
    

load(paste0(output_dir,"/",out_base,"_gamye_iCAR.RData"))

summ = stanfit$summary()

summ <- summ %>% 
  mutate(prior_scale = prior_scale,
         distribution = pp)

nsmooth_samples <- posterior_samples(stanfit,
                                 parm = "nsmooth",
                                 dims = c("Year_Index"))



BETA_samples <- posterior_samples(stanfit,
                                  parm = "BETA",
                                  dims = c("k"))
BETA_wide <- BETA_samples %>% 
  pivot_wider(.,id_cols = .draw,
              names_from = k,
              names_prefix = "BETA",
              values_from = .value)


nsmooth_samples <- nsmooth_samples %>% 
  left_join(., BETA_wide,by = ".draw") %>% 
  mutate(prior_scale = prior_scale,
         distribution = pp)


# 
# nsmooth_sel <- nsmooth_samples %>% 
#   filter(BETA12 < 5 , BETA12 > -5)
# 
# nsmooth_plot <- ggplot(data = nsmooth_sel,
#                        aes(x = Year_Index,
#                            y = .value,
#                            group = .draw))+
#   geom_line(alpha = 0.1)+
#   #scale_y_continuous(trans = "log10")+
#   coord_cartesian(ylim = c(0,250))
# print(nsmooth_plot)






# nsmooth <- nsmooth_samples %>% 
#   posterior_sums(.,
#                  dims = c("Year_Index"))%>% 
#   mutate(prior_scale = prior_scale,
#          distribution = pp)

nyears = max(nsmooth_samples$Year_Index)
trs <- function(y1,y2,ny){
  tt <- (((y2/y1)^(1/ny))-1)*100
}

for(y2 in c(seq(2,nyears,by = 5),nyears)){
nyh <- paste0("Y",y2)
ny = y2-1
trends <- nsmooth_samples %>% 
  filter(Year_Index %in% c(1,y2)) %>% 
  select(.draw,.value,Year_Index) %>% 
  pivot_wider(.,names_from = Year_Index,
              values_from = .value,
              names_prefix = "Y") %>%
  rename_with(.,~gsub(pattern = nyh,replacement = "Y2", .x)) %>% 
  group_by(.draw) %>% 
  summarise(trend = trs(Y1,Y2,ny))%>% 
  mutate(prior_scale = prior_scale,
         distribution = pp,
         first_year = 1,
         last_year = y2,
         nyears = ny)
trends_out <- bind_rows(trends_out,trends)
}

nsmooth_out <- bind_rows(nsmooth_out,nsmooth_samples)
summ_out <- bind_rows(summ_out,summ)

  }#prior_scale
  print(paste(pp,prior_scale))
}# pp

save(file = "output/GAM_prior_sim_summary.RData",
     list = c("nsmooth_out",
              "trends_out",
              "summ_out"))





# summarize and plot ------------------------------------------------------

load("output/GAM_prior_sim_summary.RData")


bbs_trends <- read.csv("data_basic/2019All BBS trends continent and national.csv")
bbs_short <- bbs_trends %>% 
  filter(Trend_Time == "Short-term",
         Region_type == "continental",
         !grepl(pattern = "^unid",species)) %>% # removing Cave Swallow - original abundance == 0 and rate of change is not meaningful
  select(Trend,Start_year,species) %>% 
  mutate(abs_trend = abs(Trend))

bbs_long <- bbs_trends %>% 
  filter(Trend_Time == "Long-term",
         Region_type == "continental",
         !grepl(pattern = "^unid",species),
         species != "Cave Swallow",
         species != "Northern Gannet")%>%  # removing Cave Swallow - original abundance == 0 and rate of change is not meaningful
  select(Trend,Start_year,species) %>% 
  mutate(abs_trend = abs(Trend)) %>% 
  arrange(abs_trend)



hist(bbs_long$abs_trend,breaks = 20)
hist(bbs_short$abs_trend,breaks = 20)
G_short <- max(bbs_short$abs_trend)
G_long <- max(bbs_long$abs_trend)
# USGS trends -------------------------------------------------------------

bbs_trends_usgs <- read.csv("data/BBS_1966-2019_core_best_trend.csv")



bbs_trends_usgs_long <-bbs_trends_usgs %>% 
  filter(Region == "SU1",
         !grepl(pattern = "^unid",Species.Name))%>%  
  select(Trend,Species.Name) %>% 
  mutate(abs_trend = abs(Trend)) %>% 
  arrange(-abs_trend)

G_long_usgs <- max(bbs_trends_usgs_long$abs_trend)
G_short_usgs <- 23.5 #Eurasian Collared Dove trend for short-term 1966-2019 analysis USGS
    # short-term trends not included in Science Base, but visible here:
    # https://www.mbr-pwrc.usgs.gov/bbs/reglist19v3.shtml





# exploring prior sims ----------------------------------------------------



trends_out <- trends_out %>% 
  mutate(abs_trend = abs(trend))

trends_long <- trends_out %>% 
  filter(nyears == 53 ) 

trend_h <- ggplot(data = trends_long)+
  geom_histogram(aes(x = abs_trend),binwidth = 2)+
  facet_wrap(facets = vars(distribution,prior_scale),
             nrow = 3,
             ncol = 5,
             scales = "free_y")+
  #coord_cartesian(xlim = c(-20,20))+
  coord_cartesian(xlim = c(0,30))+
  theme_bw()
print(trend_h)

quant_long_tends <- trends_long %>% 
  group_by(distribution,prior_scale) %>% 
  summarise(mean_abs_t = mean(abs(trend)),
            median_abs_t = median(abs(trend)),
            U90 = quantile(abs(trend),0.90),
            U80 = quantile(abs(trend),0.80),
            U99 = quantile(abs(trend),0.99),
            pGTmax = length(which(abs_trend > G_long_usgs))/length(abs_trend))




trends_short <- trends_out %>% 
  filter(nyears == 11 )

trend_hs <- ggplot(data = trends_short)+
  geom_histogram(aes(x = abs_trend),binwidth = 2)+
  facet_wrap(facets = vars(distribution,prior_scale),
             nrow = 3,
             ncol = 5,
             scales = "free_y")+
  coord_cartesian(xlim = c(0,20))+
  theme_bw()
print(trend_hs)

quant_short_tends <- trends_short %>% 
  group_by(distribution,prior_scale) %>% 
  summarise(mean_abs_t = mean(abs(trend)),
            median_abs_t = median(abs(trend)),
            U90 = quantile(abs(trend),0.90),
            U80 = quantile(abs(trend),0.80),
            U99 = quantile(abs(trend),0.99),
            pGTMax = length(which(abs_trend > G_short_usgs))/length(abs_trend))








# trend_f <- ggplot(data = trends_long)+
#   geom_freqpoly(aes(x = trend, group = prior_scale,
#                     colour = prior_scale),
#                 bins = 5000)+
#   scale_color_viridis_c()+
#   coord_cartesian(xlim = c(-20,20),
#                   ylim = c(0,1000))+
#   #scale_x_continuous(limits = c(-30,30))+
#   facet_wrap(facets = vars(distribution),
#              nrow = 3,
#              scales = "fixed")
# print(trend_f)








