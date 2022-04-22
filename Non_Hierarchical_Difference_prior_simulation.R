### prior simulation of spatial GAM sd

library(bbsBayes)
library(tidyverse)
library(cmdstanr)
library(patchwork)

setwd("C:/GitHub/GAMYE_Elaboration")


# fit model with fixed data for all parameters except the local sm --------

species <- "Yellow-headed Blackbird"
species_f <- gsub(species,pattern = " ",replacement = "_")



for(pp in c("t3","t4","t10")){
  for(prior_scale in c(0.1,0.2,0.3)){
    
  
    if(pp == "t3"){
      pnorm <- 0
      df = 3
    }
    if(pp == "t4"){
      pnorm <- 0
      df = 4
    }
    if(pp == "t10"){
      pnorm <- 0
      df = 10
    }
    
    
    tp = paste0(pp,"_rate_",prior_scale)
    
    #STRATA_True <- log(2)
    output_dir <- "output/"
    out_base <- paste0("Prior_sim_Non_Hierarchical_Difference_",tp,"_BBS")
    csv_files <- paste0(out_base,"-",1:3,".csv")
    
   # if(!file.exists(paste0(output_dir,csv_files[1]))){
      
      load(paste0("Data/Real_data_",species_f,"_BBS.RData"))
      
      tmp_data = original_data_df
      
      nstrata = max(strata_df$Stratum_Factored)
      nyears = max(tmp_data$Year_Index)
      nyears_m1 = nyears-1
      midyear = floor(nyears/2)
      
      Iy1 = c((midyear-1):1)
      Iy2 = c((midyear+1):nyears)
      nIy1 = length(Iy1)
      nIy2 = length(Iy2)
      
      stan_data = list(#scalar indicators
        nstrata = nstrata,
        nyears = nyears,
        nyears_m1 = nyears_m1,
        
        
        #temporal indexing
        midyear = midyear,
        Iy1 = Iy1,
        Iy2 = Iy2,
        nIy1 = nIy1,
        nIy2 = nIy2,
        
        
        prior_scale_b = prior_scale,
        pnorm = pnorm,
        df = df,
        
        #vector of zeros to fill midyear beta values
        zero_betas = rep(0,nstrata)
      )
      
      
      
      
      # Fit model ---------------------------------------------------------------
      
      print(paste("beginning",out_base,Sys.time()))
      
      mod.file = "models/Difference_NonHierarchical_Prior_sim.stan"
      
      ## compile model
      model <- cmdstan_model(mod.file)
      
      
      # Initial Values ----------------------------------------------------------
      
      
      init_def <- function(){ list(sdbeta = runif(nyears_m1,0.01,0.1),
                                   beta_raw = matrix(rnorm(nyears_m1*nstrata,0,0.01),nrow = nstrata,ncol = nyears_m1))}
      
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
      
      
      
      save(list = c("stanfit","stan_data","csv_files",
                    "out_base"),
           file = paste0(output_dir,"/",out_base,"_difference.RData"))
      
      
      
    #}
    
  }#end prior_scale loop
}#end pp loop



# post model summary of priors --------------------------------------------

source("Functions/posterior_summary_functions.R")

n_out <- NULL
N_out <- NULL
trends_out <- NULL
summ_out <- NULL

for(pp in c("t3","t4","t10")){
  for(prior_scale in c(0.1,0.2,0.3)){
    
    
    
    if(pp == "t3"){
      pnorm <- 0
      df = 3
    }
    if(pp == "t4"){
      pnorm <- 0
      df = 4
    }
    if(pp == "t10"){
      pnorm <- 0
      df = 10
    }
    
    
    tp = paste0(pp,"_rate_",prior_scale)
    
    #STRATA_True <- log(2)
    output_dir <- "output/"
    out_base <- paste0("Prior_sim_Non_Hierarchical_Difference_",tp,"_BBS")
    csv_files <- paste0(out_base,"-",1:3,".csv")
    
    
    load(paste0(output_dir,"/",out_base,"_difference.RData"))
    
    summ = stanfit$summary()
    
    summ <- summ %>% 
      mutate(prior_scale = prior_scale,
             distribution = pp)
    
    
    n_samples <- posterior_samples(stanfit,
                                         parm = "n",
                                         dims = c("Stratum_Factored","Year_Index"))
    
    n <- n_samples %>% 
      posterior_sums(.,
                     dims = c("Stratum_Factored","Year_Index"))%>% 
      mutate(prior_scale = prior_scale,
             distribution = pp,
             param = "n")
   
    

    
    nyears = max(n_samples$Year_Index)
    # function to calculate a %/year trend from a count-scale trajectory
    trs <- function(y1,y2,ny){
      tt <- (((y2/y1)^(1/ny))-1)*100
    }
    
    for(tl in c(2,6,11,21,nyears)){ #estimating all possible 1-year, 10-year, and full trends
      ny = tl-1
      yrs1 <- seq(1,(nyears-ny),by = ny)
      yrs2 <- yrs1+ny
      for(j in 1:length(yrs1)){
        y2 <- yrs2[j]
        y1 <- yrs1[j]
        
        nyh2 <- paste0("Y",y2)
        nyh1 <- paste0("Y",y1)
        trends <- n_samples %>% 
          filter(Year_Index %in% c(y1,y2)) %>% 
          select(.draw,.value,Stratum_Factored,Year_Index) %>% 
          pivot_wider(.,names_from = Year_Index,
                      values_from = .value,
                      names_prefix = "Y") %>%
          rename_with(.,~gsub(pattern = nyh2,replacement = "YE", .x)) %>% 
          rename_with(.,~gsub(pattern = nyh1,replacement = "YS", .x)) %>% 
          group_by(.draw,Stratum_Factored) %>% 
          summarise(trend = trs(YS,YE,ny))%>% 
          mutate(prior_scale = prior_scale,
                 distribution = pp,
                 first_year = y1,
                 last_year = y2,
                 nyears = ny,
                 param = "n")
        trends_out <- bind_rows(trends_out,trends)
        
        
        # TRENDS <- N_samples %>% 
        #   filter(Year_Index %in% c(y1,y2)) %>% 
        #   select(.draw,.value,Year_Index) %>% 
        #   pivot_wider(.,names_from = Year_Index,
        #               values_from = .value,
        #               names_prefix = "Y") %>%
        #   rename_with(.,~gsub(pattern = nyh2,replacement = "YE", .x)) %>% 
        #   rename_with(.,~gsub(pattern = nyh1,replacement = "YS", .x)) %>% 
        #   group_by(.draw) %>% 
        #   summarise(trend = trs(YS,YE,ny))%>% 
        #   mutate(prior_scale = prior_scale,
        #          distribution = pp,
        #          first_year = y1,
        #          last_year = y2,
        #          nyears = ny,
        #          param = "N")
        # trends_out <- bind_rows(trends_out,TRENDS)
      }
    }
    
    
    
    n_out <- bind_rows(n_out,n)
    # N_out <- bind_rows(N_out,N)
     summ_out <- bind_rows(summ_out,summ)
    print(paste(pp,prior_scale))
    
  }#prior_scale
}# pp

save(file = "output/Non_Hierarchical_Difference_prior_sim_summary.RData",
     list = c("n_out",
              "trends_out",
              "summ_out"))



# summarise ---------------------------------------------------------------

load("output/Non_Hierarchical_Difference_prior_sim_summary.RData")


trends_sd <- trends_out %>% 
  group_by(.draw,prior_scale,distribution,param,nyears) %>% 
  summarise(sd_trends = sd(trend),
            min_trend = min(trend),
            max_trend = max(trend),
            lci95 = quantile(trend,0.025),
            uci95 = quantile(trend,0.975),
            Absuci50 = quantile(abs(trend),0.5),
            Absuci95 = quantile(abs(trend),0.975),
            Absuci75 = quantile(abs(trend),0.75),
            span_95 = uci95 - lci95,
            range_trend = max_trend - min_trend) %>% 
  mutate(scale_factor = factor(prior_scale,levels = c(0.5,1,2,3,4),ordered = TRUE))

# save(list = "trends_sd",
#      file = "data/long_term_trend_hier_Non_Spatial_sd_prior_sim.RData")


trends_sd_t <- trends_sd %>% 
  filter(nyears == 53,
         param == "n")


absplot_long <- ggplot(data = trends_sd_t,
                  aes(Absuci50,after_stat(density),
                      group = prior_scale,
                      colour = prior_scale))+
            geom_freqpoly(breaks = seq(0,700,0.5),center = 0)+
            xlab("SD of short-term BBS trends among regions from CWS models (2009-2019)")+
            theme_bw()+
            scale_y_continuous(limits = c(0,1))+
  coord_cartesian(xlim = c(0,20))
  
print(absplot_long)        



#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 


# TRENDS_abs <- trends_out %>% 
#   filter(is.na(Stratum_Factored),
#          distribution != "t4",
#          prior_scale != 0.5) %>% 
#   mutate(abs_trend = abs(trend),
#          scale_factor = factor(prior_scale,ordered = TRUE),
#          nyears_factor = factor(nyears,levels = c(1,5,10,20,53),ordered = TRUE))

trends_abs <- trends_out %>% 
  filter(!is.na(Stratum_Factored),
         distribution != "t4",
         prior_scale != 0.5) %>% 
  mutate(abs_trend = abs(trend),
         scale_factor = factor(prior_scale,ordered = TRUE),
         nyears_factor = factor(nyears,levels = c(1,5,10,20,53),ordered = TRUE))


# TRENDS_abs_l <- TRENDS_abs %>% 
#   filter(nyears == 53,
#          prior_scale != 0.5)
# TRENDS_abs_20 <- TRENDS_abs %>% 
#   filter(nyears == 20,
#          prior_scale != 0.5)
# TRENDS_abs_10 <- TRENDS_abs %>% 
#   filter(nyears == 10,
#          prior_scale != 0.5)
# TRENDS_abs_5 <- TRENDS_abs %>% 
#   filter(nyears == 5,
#          prior_scale != 0.5)
# TRENDS_abs_1 <- TRENDS_abs %>% 
#   filter(nyears == 1,
#          prior_scale != 0.5)


      
# 
# ABSplot_long <- ggplot(data = TRENDS_abs,
#                        aes(abs_trend,after_stat(density),
#                            group = scale_factor,
#                            colour = scale_factor))+
#   geom_freqpoly(breaks = seq(0,700,0.5),center = 0)+
#   xlab("Absolute value of long-term trends from prior distribution")+
#   theme_bw()+
#   #scale_y_continuous(limits = c(0,1))+
#   coord_cartesian(xlim = c(0,20))+
#   scale_colour_viridis_d()+
#   facet_wrap(vars(distribution,nyears_factor),ncol = 5,nrow = 2)
# 
# print(ABSplot_long)


absplot_long <- ggplot(data = trends_abs,
                       aes(abs_trend,after_stat(density),
                           group = scale_factor,
                           colour = scale_factor))+
  geom_freqpoly(breaks = seq(0,700,0.5),center = 0)+
  xlab("Absolute value of long-term trends from prior distribution")+
  theme_bw()+
  scale_y_continuous(limits = c(0,1))+
  coord_cartesian(xlim = c(0,50))+
  scale_colour_viridis_d()+
  facet_wrap(vars(distribution,nyears_factor),ncol = 5,nrow = 2)


print(absplot_long)


#### Generalize this plot for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 



tplot_long <- ggplot(data = trends_abs,
                     aes(trend,after_stat(density),
                           group = scale_factor,
                           colour = scale_factor))+
  geom_freqpoly(breaks = seq(-99,700,0.5),center = 0)+
  xlab("Absolute value of long-term trends from prior distribution")+
  theme_bw()+
  #scale_y_continuous(limits = c(0,1))+
  coord_cartesian(xlim = c(-50,100))+
  scale_colour_viridis_d()+
  facet_wrap(vars(distribution,nyears_factor),ncol = 5,nrow = 2,
             scales = "free")


print(tplot_long)



#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 
#### Generalize this part for all simulations 









trends_sd_t <- trends_sd %>% 
  filter(distribution == "t",
         nyears == 10)

overp_t <- realised_long_bbs_hist +
  geom_freqpoly(data = trends_abs,
                aes(sd_trends,after_stat(density),
                    colour = scale_factor),
                breaks = seq(0,100,1),center = 0,
                alpha = 0.8)+
  scale_colour_viridis_d(begin = 0.5,alpha = 0.8,
                         "Prior Scale\nSD half-t df=3")+
  coord_cartesian(xlim = c(0,20))

print(overp_t)


trends_sd_norm <- trends_sd %>% 
  filter(distribution == "norm",
         nyears == 53)


overp_norm <- realised_long_bbs_hist +
  geom_freqpoly(data = trends_sd_norm,
                aes(sd_trends,after_stat(density),
                    colour = scale_factor),
                breaks = seq(0,100,1),center = 0,
                alpha = 0.8)+
  scale_colour_viridis_d(begin = 0.5,alpha = 0.8,
                         "Prior Scale\nSD half-normal")+
  coord_cartesian(xlim = c(0,20))

print(overp_norm)



trends_sd_gamma <- trends_sd %>% 
  filter(distribution == "gamma",
         nyears == 53)


overp_gamma <- realised_long_bbs_hist +
  geom_freqpoly(data = trends_sd_gamma,
                aes(sd_trends,after_stat(density),
                    colour = scale_factor),
                breaks = seq(0,100,1),center = 0,
                alpha = 0.8)+
  scale_colour_viridis_d(begin = 0.5,alpha = 0.8,direction = -1,
                         "Prior Scale\nSD gamma shape = 2")+
  coord_cartesian(xlim = c(0,20))

print(overp_gamma)

print(overp_t/overp_gamma/overp_norm)



trends_sd_t1 <- trends_sd_t %>% 
  filter(prior_scale == 1)
trends_sd_norm1 <- trends_sd_norm %>% 
  filter(prior_scale == 1)
trends_sd_gamma1 <- trends_sd_gamma %>% 
  filter(prior_scale == 2)


overp_t1 <- realised_long_bbs_hist +
  geom_freqpoly(data = trends_sd_t1,
                aes(sd_trends,after_stat(density)),
                    colour = "darkgreen",
                breaks = seq(0,100,1),center = 0,
                alpha = 0.8)+
  coord_cartesian(xlim = c(0,20))

print(overp_t1)

overp_norm1 <- realised_long_bbs_hist +
  geom_freqpoly(data = trends_sd_norm1,
                aes(sd_trends,after_stat(density)),
                    colour = "darkgreen",
                breaks = seq(0,100,1),center = 0,
                alpha = 0.8)+
  coord_cartesian(xlim = c(0,20))

print(overp_norm1)

overp_gamma1 <- realised_long_bbs_hist +
  geom_freqpoly(data = trends_sd_gamma1,
                aes(sd_trends,after_stat(density)),
                colour = "darkgreen",
                breaks = seq(0,100,1),center = 0,
                alpha = 0.8)+
  coord_cartesian(xlim = c(0,20))

print(overp_gamma1)

# save(list = c("overp_gamma1","overp_norm1","overp_t1"),
#      file = "data/sd_GAM_spatial_saved_long-term.RData")

trends_long <- trends_sd %>% 
  filter(nyears == 53)

quant_long_tends <- trends_long %>% 
  group_by(distribution,prior_scale) %>% 
  summarise(median_sim_sd_trends = median(sd_trends),
            sim_80th = quantile(sd_trends,0.80),
            sim_90th = quantile(sd_trends,0.90),
            sim_99th = quantile(sd_trends,0.99),
            proportion_sim_GT_max = length(which(sd_trends > G_long_usgs))/length(sd_trends))


kable(quant_long_tends, booktabs = TRUE,
      digits = 3) %>%
  kable_styling(font_size = 8,
repeat_header_text = "Prior simulated distribution quantiles and proportion of the distributions greater than the realised maximum values, for three prior distributions with 5 different scales/rates")


# Short-term comparison ---------------------------------------------------



trends_sd_t_short <- trends_sd %>% 
  filter(distribution == "t",
         nyears == 10)


overp_t_short <- realised_short_bbs_hist +
  geom_freqpoly(data = sd_slope_trends_short,
                breaks = seq(0,G_short_slope_canada,0.5),center = 0,
                colour = grey(0.5))+
  geom_freqpoly(data = trends_sd_t_short,
                aes(sd_trends,after_stat(density),
                    colour = scale_factor),
                breaks = seq(0,100,1),center = 0,
                alpha = 0.8)+
  scale_colour_viridis_d(begin = 0.5,alpha = 0.8,
                         "Prior Scale\nSD half-t df=3")+
  coord_cartesian(xlim = c(0,20))

print(overp_t_short)


trends_sd_norm_short <- trends_sd %>% 
  filter(distribution == "norm",
         nyears == 10)


overp_norm_short <- realised_short_bbs_hist +
  geom_freqpoly(data = trends_sd_norm_short,
                aes(sd_trends,after_stat(density),
                    colour = scale_factor),
                breaks = seq(0,100,1),center = 0,
                alpha = 0.8)+
  geom_freqpoly(data = sd_slope_trends_short,
                breaks = seq(0,G_short_slope_canada,0.5),center = 0,
                colour = grey(0.5))+
  scale_colour_viridis_d(begin = 0.5,alpha = 0.8,
                         "Prior Scale\nSD half-normal")+
  coord_cartesian(xlim = c(0,20))

print(overp_norm_short)



trends_sd_gamma_short <- trends_sd %>% 
  filter(distribution == "gamma",
         nyears == 10)


overp_gamma_short <- realised_short_bbs_hist +
  geom_freqpoly(data = trends_sd_gamma_short,
                aes(sd_trends,after_stat(density),
                    colour = scale_factor),
                breaks = seq(0,100,1),center = 0,
                alpha = 0.8)+
  geom_freqpoly(data = sd_slope_trends_short,
                breaks = seq(0,G_short_slope_canada,0.5),center = 0,
                colour = grey(0.5))+
  scale_colour_viridis_d(begin = 0.5,alpha = 0.8,direction = -1,
                         "Prior Scale\nSD gamma shape = 2")+
  coord_cartesian(xlim = c(0,20))

print(overp_gamma_short)

print(overp_t_short/overp_gamma_short/overp_norm_short)

trends_short <- trends_sd %>% 
  filter(nyears == 10)

quant_short_tends <- trends_short %>% 
  group_by(distribution,prior_scale) %>% 
  summarise(median_sim_sd_trends = median(sd_trends),
            sim_80th = quantile(sd_trends,0.80),
            sim_90th = quantile(sd_trends,0.90),
            sim_99th = quantile(sd_trends,0.99),
            proportion_sim_GT_max = length(which(sd_trends > G_short_canada))/length(sd_trends),
            proportion_sim_GT_max_slope = length(which(sd_trends > G_short_slope_canada))/length(sd_trends))

kable(quant_short_tends, booktabs = TRUE,
      digits = 3) %>%
  kable_styling(font_size = 8,
                repeat_header_text = "Prior simulated distribution quantiles and proportion of the distributions greater than the realised maximum values, for three prior distributions with 5 different scales/rates")

trends_sd_t_short1 <- trends_sd_t_short %>% 
  filter(prior_scale == 1)
overp_t_short1 <- realised_short_bbs_hist +
  geom_freqpoly(data = sd_slope_trends_short,
                breaks = seq(0,G_short_slope_canada,0.5),center = 0,
                colour = grey(0.5))+
  geom_freqpoly(data = trends_sd_t_short1,
                aes(sd_trends,after_stat(density)),
                colour = "darkgreen",
                breaks = seq(0,100,1),center = 0,
                alpha = 0.8)+
  scale_colour_viridis_d(begin = 0.5,alpha = 0.8,
                         "Prior Scale\nSD half-t df=3")+
  coord_cartesian(xlim = c(0,20))

print(overp_t_short1)
# save(list = c("quant_short_tends",
#               "overp_t1",
#               "overp_t_short",
#               "overp_gamma_short",
#               "overp_norm_short",
#               "overp_t_short1",
#               "quant_long_tends",
#               "overp_t",
#               "overp_gamma",
#               "overp_norm"),
#      file = "data/SD_prior_sim_graphs.RData")
