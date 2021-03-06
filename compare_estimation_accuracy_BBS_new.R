## compare estimation accuracy

library(bbsBayes)
library(tidyverse)
library(cmdstanr)
library(posterior)

library(sf)
library(geofacet)
source("functions/indices_cmdstan.R")
source("functions/posterior_summary_functions.R")
source("Functions/palettes.R")

#species = "Pine Warbler"  
species = "Yellow-headed Blackbird"  
species_f <- gsub(species,pattern = " ",replacement = "_")


tp = "non_linear"
mk = ""

 
MAs <- round(log(c(0.1,0.5,1,5,10,50)),2)# true mean abundances for different simulations


for(ma in MAs){  
  load(paste0("Data/Simulated_data_",ma,"_",tp,"_BBS.RData"))
  
  
  for(sns in c("nonSpatial_alt_","")){
    
  output_dir <- "output/"
  out_base <- paste0("sim_",sns,tp,"_",ma,"_BBS")
 
  
  load(paste0(output_dir,"/",out_base,"_gamye_iCAR.RData"))

    out_base_sim <- out_base
  csv_files <- paste0(output_dir,out_base,"-",1:3,".csv")
  

stanfit <- as_cmdstan_fit(files = csv_files)


# Compare stratum level beta parameters -----------------------------------
to_save <- c("out_base","out_base_sim")


betas_est <- posterior_samples(stanfit,
                               parm = "beta",
                               dims = c("Stratum_Factored","k")) %>% 
posterior_sums(.,
               dims = c("Stratum_Factored","k")) 

beta_comp <- Beta_True %>% 
  as.data.frame() %>% 
  mutate(k = 1:13) %>% 
  pivot_longer(cols = starts_with("V"),
               names_prefix = "V",
               names_to = "Stratum_Factored",
               values_to = "beta_True") %>% 
  mutate(Stratum_Factored = as.integer(Stratum_Factored)) %>% 
  left_join(betas_est,by = c("Stratum_Factored","k"))

to_save <- c(to_save,"beta_comp")


betas_plot <- ggplot(data = beta_comp,aes(x = beta_True,y = mean))+
  geom_point(aes(colour = Stratum_Factored))+
  scale_colour_viridis_c()+
  geom_errorbar(aes(ymin = Q_025,ymax = Q_975),width = 0,alpha = 0.2)+
  geom_abline(slope = 1, intercept = 0)+ 
  facet_wrap(~k,nrow = ceiling(sqrt(stan_data$nknots_year)),
             ncol = ceiling(sqrt(stan_data$nknots_year)))


sdbetas_est <- posterior_samples(stanfit,
                               parm = "sdbeta",
                               dims = c("k")) %>% 
  posterior_sums(.,
                 dims = c("k")) 

# compare stratum level smoothed indices ----------------------------------

true_inds <- balanced %>% 
  select(Stratum_Factored,Smooth,Year,True_log_traj,True_strata_effects) %>% 
  distinct() %>% 
  mutate(True_nsmooth = exp(Smooth + True_strata_effects)) %>% 
  arrange(Stratum_Factored)

if(mk == ""){

  obs_means <- realised %>% 
  group_by(Stratum_Factored,Year) %>% 
  summarise(n_routes = n(),
            mean_count = mean(count),
            lmean_count = mean(log(count+1)),
            .groups = "drop") 
  
}else{
  obs_means <- realised_mask %>% 
    group_by(Stratum_Factored,Year) %>% 
    summarise(n_routes = n(),
              mean_count = mean(count),
              lmean_count = mean(log(count+1)),
              .groups = "drop") 
}

strat_log_smooths <- GAM_year$Year_basis %*% Beta_True 

true_log_smooths <- as.data.frame(strat_log_smooths)
true_log_smooths[,"Year"] <- min_year : 2019
log_smooth_plot <- true_log_smooths %>% 
  pivot_longer(cols = starts_with("V"),
               names_to = "Stratum_Factored",
               names_prefix = "V",
               values_to = "Smooth")%>% 
  mutate(index = exp(Smooth),
         Stratum_Factored = as.integer(Stratum_Factored)) %>% 
  left_join(.,strata_df,by = c("Stratum_Factored")) %>% 
  select(Year,Stratum_Factored,Smooth,index,Stratum) %>% 
  arrange(Year,Stratum_Factored) 


smooth_samples <- posterior_samples(stanfit,
                               parm = "smooth_pred",
                               dims = c("Year_Index","Stratum_Factored")) %>%
   mutate(Year = Year_Index+min(balanced$Year)-1) %>% 
   left_join(.,log_smooth_plot,by = c("Year","Stratum_Factored")) %>% 
   mutate(smooth_err = .value - Smooth)

 
smooth_comp <- smooth_samples %>% 
  group_by(Year,Stratum_Factored) %>% 
  summarise(smooth_est = mean(.value,na.rm = T),
            smooth_est_lci = quantile(.value,0.025,na.rm = T),
            smooth_est_uci = quantile(.value,0.975,na.rm = T),
            smooth_error = mean(smooth_err,na.rm = T),
            smooth_error_lci = quantile(smooth_err,0.025,na.rm = T),
            smooth_error_uci = quantile(smooth_err,0.975,na.rm = T),
            smooth_abs_error = mean(abs(smooth_err),na.rm = T),
            smooth_abs_error_lci = quantile(abs(smooth_err),0.025,na.rm = T),
            smooth_abs_error_uci = quantile(abs(smooth_err),0.975,na.rm = T))

to_save <- c(to_save,"smooth_comp")



strat_inds_smooth <- index_function(fit = stanfit,
                                    parameter = "nsmooth",
                                    year_1 = min(balanced$Year),
                                    strat = "Stratum_Factored")


# Trends ------------------------------------------------------------------

all_trends <- NULL

tyrs = unique(c(2014,2009,2004,1999,1994,1990,1985,1980,1975,1970,1966))
tyrs2 <- rep(2019,length(tyrs))
tyrs2 <- c(tyrs2,tyrs+5)
tyrs <- c(tyrs,tyrs)

for(j in 1:length(tyrs)){
  yy = tyrs[j]
  yy2 = tyrs2[j]
  
  TT <- trends_function(ind_list = strat_inds_smooth,
                        start_year = yy,
                        end_year = yy2) %>% 
    mutate(species = species,
           first_year = yy,
           last_year = yy2)
  

  tt_true <- true_inds %>% 
    select(Stratum_Factored,Year,True_nsmooth) %>% 
    mutate(Year = as.character(Year),
           True_nsmooth = as.numeric(True_nsmooth)) %>% 
    filter(Year %in% as.character(c(yy,yy2))) %>% 
    pivot_wider(names_from = Year,
                values_from = True_nsmooth,
                names_prefix = "Y") %>% 
    rename_with(., ~gsub(replacement = "start",
                         pattern = paste0("Y",yy),.x,
                         fixed = TRUE)) %>% 
    rename_with(., ~gsub(replacement = "end",
                         pattern = paste0("Y",yy2),.x,
                         fixed = TRUE)) %>% 
    group_by(Stratum_Factored) %>% 
    summarise(true_trend = texp(end/start,ny =(yy2-yy)))
              
  
  TT <- TT %>% 
    left_join(.,tt_true,by = "Stratum_Factored")
  
  all_trends <- bind_rows(all_trends,TT)
  
}


to_save <- c(to_save,"all_trends")


nsmooth_est <- strat_inds_smooth$indices %>% 
  mutate(version = "smooth",
         Year = true_year)


nsmooth_comp <- nsmooth_est %>% 
  left_join(true_inds,by = c("Stratum_Factored","Year"))



nsmooth_plot = ggplot(data = nsmooth_comp,aes(x = True_nsmooth,
                                              y = mean))+
  geom_point(aes(colour = Year))+
  scale_colour_viridis_c()+
  geom_errorbar(aes(ymin = lci,ymax = uci),width = 0,alpha = 0.2)+
  geom_abline(slope = 1, intercept = 0)+ 
  facet_wrap(~Stratum_Factored,nrow = ceiling(sqrt(stan_data$nstrata)),
             ncol = ceiling(sqrt(stan_data$nstrata)))




true_inds <- true_inds %>% 
  mutate(version = "True")

nsmooth_comp2 <- nsmooth_est %>% 
  select(mean,lci,uci,Stratum_Factored,Year) %>% 
  rename(True_nsmooth = mean) %>% 
  mutate(version = "Estimated") %>% 
  bind_rows(.,true_inds)%>% 
  left_join(.,obs_means,by = c("Stratum_Factored","Year"))
  

to_save <- c(to_save,"nsmooth_comp2","obs_means")

if(mk == ""){
nsmooth_plot2 = ggplot(data = nsmooth_comp2,aes(y = True_nsmooth,
                                              x = Year))+
  geom_ribbon(aes(ymin = lci,ymax = uci,fill = version),alpha = 0.2)+
  geom_line(aes(colour = version))+
  scale_colour_viridis_d(aesthetics = c("colour","fill"))+
  scale_y_continuous(limits = c(0,NA))+
  facet_wrap(~Stratum_Factored,nrow = ceiling(sqrt(stan_data$nstrata)),
             ncol = ceiling(sqrt(stan_data$nstrata)),
             scales = "free_y")
}else{
  relab <- function(x,m){
    y = vector("character",
               length(x))
    for(i in 1:length(x)){
      if(x[i] == "True"){
      y[i] = x[i]
    }else{
    if(m[i]){
      y[i] = (paste(x[i],"Masked"))
    }else{
      y[i] = x[i]
    }
    }
    }
    
    return(y)
  }
  
  nsmooth_comp2 <- nsmooth_comp2 %>% 
    left_join(.,strata_mask,by = "Stratum_Factored") %>% 
    mutate(version = relab(version,masked))
  
  
  
  nsmooth_plot2 = ggplot(data = nsmooth_comp2,aes(y = True_nsmooth,
                                                x = Year))+
  geom_ribbon(aes(ymin = lci,ymax = uci,fill = version),alpha = 0.2)+
  geom_line(aes(colour = version))+
  scale_colour_viridis_d(aesthetics = c("colour","fill"))+
  scale_y_continuous(limits = c(0,NA))+
  facet_wrap(~Stratum_Factored,nrow = ceiling(sqrt(stan_data$nstrata)),
             ncol = ceiling(sqrt(stan_data$nstrata)),
             scales = "free_y")


  print(nsmooth_plot2)
  # nsmooth_mask <- nsmooth_comp2 %>% 
#   filter(masked == TRUE)

# nsmooth_plot3 = ggplot(data = nsmooth_mask,aes(y = True_nsmooth,
#                                                 x = Year))+
#   geom_ribbon(aes(ymin = lci,ymax = uci,fill = version),alpha = 0.2)+
#   geom_line(aes(colour = version))+
#   scale_colour_viridis_d(aesthetics = c("colour","fill"))+
#   scale_y_continuous(limits = c(0,NA))+
#   facet_wrap(~Stratum_Factored,nrow = ceiling(sqrt(stan_data$nstrata)),
#              ncol = ceiling(sqrt(stan_data$nstrata)),
#              scales = "free_y")


}



# Hyperparameters ---------------------------------------------------------


true_BETA <- data.frame(True_BETA = BETA_True,
                        k = 1:stan_data$nknots_year)

BETA_est <- posterior_samples(stanfit,
                                           parm = "BETA",
                                           dims = c("k")) %>% 
  posterior_sums(.,
                 dims = c("k")) 

BETA_comp <- BETA_est %>% 
  left_join(true_BETA,by = c("k"))

to_save <- c(to_save,"BETA_comp")


BETA_plot = ggplot(data = BETA_comp,aes(x = True_BETA,
                                              y = mean))+
  geom_point(aes(colour = k))+
  scale_colour_viridis_c()+
  geom_errorbar(aes(ymin = Q_025,ymax = Q_975),width = 0,alpha = 0.2)+
  geom_abline(slope = 1, intercept = 0)



# Hyperparameters realised ------------------------------------------------


true_BETA_r <- data.frame(True_BETA = rowMeans(Beta_True),
                        k = 1:stan_data$nknots_year)



BETA_comp_r <- BETA_est %>% 
  left_join(true_BETA_r,by = c("k"))



BETA_plot_r = ggplot(data = BETA_comp_r,aes(x = True_BETA,
                                        y = mean))+
  geom_point(aes(colour = k))+
  labs(title = "Realised True BETA")+
  scale_colour_viridis_c()+
  geom_errorbar(aes(ymin = Q_025,ymax = Q_975),width = 0,alpha = 0.2)+
  geom_abline(slope = 1, intercept = 0)



# MEan Smooth comparison --------------------------------------------------

true_SMOOTH <- data.frame(True_SMOOTH = stan_data$year_basis %*% BETA_True,
                        y = 1:stan_data$nyears) %>% 
  mutate(version = "TRUE",
         esmooth = exp(True_SMOOTH),
         Year = y + min(balanced$Year)-1)



SMOOTH_est <- posterior_samples(stanfit,
                              parm = "SMOOTH_pred",
                              dims = c("y")) %>% 
  posterior_sums(.,
                 dims = c("y")) 


# Trends Hyperparamter ----------------------------------------------------


eSMOOTH_est <- posterior_samples(stanfit,
                                parm = "SMOOTH_pred",
                                dims = c("y")) %>% 
  mutate(.value = exp(.value),
         Year = y + min(balanced$Year)-1) 

for(j in 1:length(tyrs)){
  yy = tyrs[j]
  yy2 = tyrs2[j]
  
nyrs = yy2-yy
lu = 0.025
uu = 0.975
tt_tmp <- eSMOOTH_est %>%
  select(-y) %>% 
  filter(Year %in% c(yy,yy2)) %>% 
  pivot_wider(names_from = Year,
              values_from = .value,
              names_prefix = "Y") %>% 
  rename_with(., ~gsub(replacement = "start",
                       pattern = paste0("Y",yy),.x,
                       fixed = TRUE))%>% 
  rename_with(., ~gsub(replacement = "end",
                       pattern = paste0("Y",yy2),.x,
                       fixed = TRUE))%>% 
  group_by(.draw) %>% 
  summarise(t = texp(end/start,ny = nyrs),
            ch = chng(end/start),
            .groups = "drop") %>% 
  summarise(trend = mean(t),
            lci = quantile(t,lu,names = FALSE),
            uci = quantile(t,uu,names = FALSE),
            percent_change = median(ch),
            p_ch_lci = quantile(ch,lu,names = FALSE),
            p_ch_uci = quantile(ch,uu,names = FALSE),
            prob_decline = prob_dec(ch,0),
            prob_decline_GT30 = prob_dec(ch,-30),
            prob_decline_GT50 = prob_dec(ch,-50),
            prob_decline_GT70 = prob_dec(ch,-70))

TT_true <- true_SMOOTH %>% 
  select(Year,esmooth) %>% 
  mutate(Year = as.character(Year)) %>% 
  filter(Year %in% as.character(c(yy,yy2))) %>% 
  pivot_wider(names_from = Year,
              values_from = esmooth,
              names_prefix = "Y") %>% 
  rename_with(., ~gsub(replacement = "start",
                       pattern = paste0("Y",yy),.x,
                       fixed = TRUE)) %>% 
  rename_with(., ~gsub(replacement = "end",
                       pattern = paste0("Y",yy2),.x,
                       fixed = TRUE)) %>% 
  summarise(true_trend = texp(end/start,ny =(yy2-yy))) %>% 
  mutate(first_year = yy,
         last_year = yy2,
         Region_type = "Survey_Wide_Mean")


tt_tmp <- tt_tmp %>% 
  bind_cols(.,TT_true)

all_trends <- bind_rows(all_trends,
                        tt_tmp)
}

tmp = all_trends %>% 
  filter(Region_type == "Survey_Wide_Mean")

SMOOTH_comp <- SMOOTH_est %>% 
  select(mean,Q_025,Q_975,y) %>% 
  rename(True_SMOOTH = mean) %>% 
  mutate(version = "Estimated") %>% 
  bind_rows(.,true_SMOOTH)

to_save <- c(to_save,"SMOOTH_comp")

SMOOTH_plot = ggplot(data = SMOOTH_comp,aes(y = True_SMOOTH,
                                                x = y))+
  geom_ribbon(aes(ymin = Q_025,ymax = Q_975,fill = version),alpha = 0.2)+
  geom_line(aes(colour = version))+
  scale_colour_viridis_d(aesthetics = c("colour","fill"))


# Hypersmooth realised ----------------------------------------------------



true_SMOOTH_r <- data.frame(True_SMOOTH = stan_data$year_basis %*% rowMeans(Beta_True),
                          y = 1:stan_data$nyears) %>% 
  mutate(version = "TRUE")




SMOOTH_comp_r <- SMOOTH_est %>% 
  select(mean,Q_025,Q_975,y) %>% 
  rename(True_SMOOTH = mean) %>% 
  mutate(version = "Estimated") %>% 
  bind_rows(.,true_SMOOTH_r)


SMOOTH_plot_r = ggplot(data = SMOOTH_comp_r,aes(y = True_SMOOTH,
                                            x = y))+
  geom_ribbon(aes(ymin = Q_025,ymax = Q_975,fill = version),alpha = 0.2)+
  labs(title = "Realised True SMOOTH")+
  geom_line(aes(colour = version))+
  scale_colour_viridis_d(aesthetics = c("colour","fill"))



pdf(paste0("Figures/Estimated_True_Comparisons_",out_base_sim,".pdf"),
    width = 11,
    height = 8)
print(betas_plot)
print(BETA_plot)
print(BETA_plot_r)
print(nsmooth_plot)
print(nsmooth_plot2)
print(SMOOTH_plot)
print(SMOOTH_plot_r)

dev.off()

save(list = to_save,
     file = paste0("data/",out_base_sim,"_accuracy_comp.RData"))

}

  }
  


