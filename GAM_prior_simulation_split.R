### prior simulation of spatial GAM sd


setwd("C:/GitHub/GAMYE_Elaboration")

# this is all set-up
library(tidyverse)
library(cmdstanr)
library(mgcv) 



for(pp in c("norm","t")){
  for(prior_scale in c(0.5,1,2,3,4)){
    
    tp = paste0("GAM_split_",pp,prior_scale,"_rate")
    
    #STRATA_True <- log(2)
    output_dir <- "output/"
    out_base <- tp
    csv_files <- paste0(out_base,"-",1:3,".csv")
    
    
    if(pp == "norm"){
      pnorm <- 1
    }
    if(pp == "t"){
      pnorm <- 0
    }
    
    #if(!file.exists(paste0(output_dir,csv_files[1]))){
      
      nyears = 54 #to match the time-scales of BBS and CBC analyses
      dat = data.frame(year = 1:nyears)
      nknots = 13  
      nknots_realised = nknots-2 #removes the two non-penalized components (mean and linear) 
      lin_component = nknots-1 # the mgcv function below orders the basis such that
                # the linear component is hte final column in the prediction matrix
      M = mgcv::smoothCon(s(year,k = nknots, bs = "tp"),data = dat,
                          absorb.cons=TRUE,#this drops the constant
                          diagonal.penalty=TRUE) ## If TRUE then the smooth is reparameterized to turn the penalty into an identity matrix, with the final diagonal elements zeroed (corresponding to the penalty nullspace). This fits with the Bayesian interpretation of the complexity penalty as the inverse of the variance of the i.i.d. collection of parameters.
      
      
      year_basis = M[[1]]$X
      
      
      
      
      stan_data = list(#scalar indicators
        nyears = nyears,
        
        
        #GAM structure
        nknots_year = nknots_realised,
        year_basis = year_basis,
        
        prior_scale = prior_scale,
        pnorm = pnorm,
        lin_component = lin_component
      )
      
      
      
      
      # Fit model ---------------------------------------------------------------
      
      print(paste("beginning",tp,Sys.time()))
      
      mod.file = "models/GAM_split_prior_sim.stan"
      
      ## compile model
      model <- cmdstan_model(mod.file)
      
      
      # Initial Values ----------------------------------------------------------
      
      
      init_def <- function(){ list(sdbeta = runif(1,0.01,0.1),
                                   BETA_raw = rnorm(nknots_realised,0,0.01))}
      
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
    
    
    
  
  
}#end prior_scale loop
}#end pp loop


# post model summary of priors --------------------------------------------

source("Functions/posterior_summary_functions.R")

nsmooth_out <- NULL
trends_out <- NULL
summ_out <- NULL

for(pp in c("norm","t")){
  for(prior_scale in c(0.5,1,2,3,4)){
    
    tp = paste0("GAM_split_",pp,prior_scale,"_rate")
    
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



nyears = max(nsmooth_samples$Year_Index)
# function to calculate a %/year trend from a count-scale trajectory
trs <- function(y1,y2,ny){
  tt <- (((y2/y1)^(1/ny))-1)*100
}

for(tl in c(2,6,11,21)){ #estimating all possible 1-year, 10-year, and full trends
  ny = tl-1
  yrs1 <- seq(1,(nyears-ny),by = ny)
  yrs2 <- yrs1+ny
  for(j in 1:length(yrs1)){
    y2 <- yrs2[j]
    y1 <- yrs1[j]
    
nyh2 <- paste0("Y",y2)
nyh1 <- paste0("Y",y1)
trends <- nsmooth_samples %>% 
  filter(Year_Index %in% c(y1,y2)) %>% 
  select(.draw,.value,Year_Index) %>% 
  pivot_wider(.,names_from = Year_Index,
              values_from = .value,
              names_prefix = "Y") %>%
  rename_with(.,~gsub(pattern = nyh2,replacement = "YE", .x)) %>% 
  rename_with(.,~gsub(pattern = nyh1,replacement = "YS", .x)) %>% 
  group_by(.draw) %>% 
  summarise(trend = trs(YS,YE,ny))%>% 
  mutate(prior_scale = prior_scale,
         distribution = pp,
         first_year = y1,
         last_year = y2,
         nyears = ny)
trends_out <- bind_rows(trends_out,trends)
}
}
nsmooth_out <- bind_rows(nsmooth_out,nsmooth_samples)
summ_out <- bind_rows(summ_out,summ)
print(paste(pp,prior_scale))

  }#prior_scale
}# pp

save(file = "output/GAM_split_prior_sim_summary.RData",
     list = c("nsmooth_out",
              "trends_out",
              "summ_out"))





# summarize and plot ------------------------------------------------------



# bbs_trends <- read.csv("data_basic/2019All BBS trends continent and national.csv")
# bbs_short <- bbs_trends %>% 
#   filter(Trend_Time == "Short-term",
#          Region_type == "continental",
#          !grepl(pattern = "^unid",species)) %>% # removing Cave Swallow - original abundance == 0 and rate of change is not meaningful
#   select(Trend,Start_year,species) %>% 
#   mutate(abs_trend = abs(Trend))
# 
# bbs_long <- bbs_trends %>% 
#   filter(Trend_Time == "Long-term",
#          Region_type == "continental",
#          !grepl(pattern = "^unid",species),
#          species != "Cave Swallow",
#          species != "Northern Gannet")%>%  # removing Cave Swallow - original abundance == 0 and rate of change is not meaningful
#   select(Trend,Start_year,species) %>% 
#   mutate(abs_trend = abs(Trend)) %>% 
#   arrange(abs_trend)
# 
# 
# 
# hist(bbs_long$abs_trend,breaks = 20)
# hist(bbs_short$abs_trend,breaks = 20)
# G_short <- max(bbs_short$abs_trend)
# G_long <- max(bbs_long$abs_trend)
# USGS trends -------------------------------------------------------------

bbs_trends_usgs <- read.csv("data/BBS_1966-2019_core_best_trend.csv")


## selecting survey-wide trend estimates
bbs_trends_usgs_long <-bbs_trends_usgs %>% 
  filter(Region == "SU1")%>%  
  select(Trend,Species.Name) %>% 
  mutate(abs_trend = abs(Trend)) %>% #calculating absolute values of the trends
  arrange(-abs_trend)

realised_long_bbs_hist <- ggplot(data = bbs_trends_usgs_long,
                            aes(abs_trend,after_stat(density)))+
  geom_freqpoly(breaks = seq(0,13,0.5),center = 0)+
  xlab("Absolute value of long-term BBS trends USGS models (1966-2019)")+
  theme_bw()
print(realised_long_bbs_hist)

G_long_usgs <- max(bbs_trends_usgs_long$abs_trend)


G_short_usgs <- 23.5 #Eurasian Collared Dove trend for short-term 1966-2019 analysis USGS
    # short-term trends not included in Science Base, but visible here:
    # https://www.mbr-pwrc.usgs.gov/bbs/reglist19v3.shtml





# exploring prior sims ----------------------------------------------------


load("output/GAM_split_prior_sim_summary.RData")

trends_out <- trends_out %>% 
  mutate(abs_trend = abs(trend))

pdf("Figures/GAM_split_simulated_trend_distributions.pdf",
    width = 11,
    height = 8.5)
for(time in c(1,5,10,20)){
trends_long <- trends_out %>% 
  filter(nyears == time ) 

trend_h <- ggplot(data = trends_long)+
  geom_histogram(aes(x = abs_trend),binwidth = 2)+
  facet_wrap(facets = vars(distribution,prior_scale),
             nrow = 3,
             ncol = 5,
             scales = "free_y")+
  labs(title = paste0("Prior distribution of absolute value of ",time,"-year trends"))+
  coord_cartesian(xlim = c(-30,30))+
  #coord_cartesian(xlim = c(0,30))+
  theme_bw()
print(trend_h)
print(time)

}

dev.off()


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








