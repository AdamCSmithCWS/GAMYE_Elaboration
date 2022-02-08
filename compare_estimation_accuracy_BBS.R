## compare estimation accuracy

library(bbsBayes)
library(tidyverse)
library(cmdstanr)
library(posterior)

source("functions/posterior_summary_functions.R")
#species = "Pine Warbler"  
species = "Yellow-headed Blackbird"  
species_f <- gsub(species,pattern = " ",replacement = "_")

for(tp in c("non_linear","linear")){
  
  load(paste0("Data/Simulated_data_",species_f,"_",tp,"_BBS.RData"))
  
  for(sns in c("","nonSpatial_","nonSpatial_alt_"))
for(mk in c("","mask_")){
  
#if(sns == "nonSpatial_" & mk == "mask_"){next}


output_dir <- "output/"
#out_base <- paste0(species_f,"_sim_",tp,"_BBS")
out_base <- paste0(species_f,"_sim_",sns,mk,tp,"_BBS")
#csv_files <- paste0(out_base,"-",1:3,".csv")
out_base_sim <- paste0("_sim_",sns,mk,tp,"_BBS")
  
load(paste0(output_dir,"/",out_base,"_gamye_iCAR.RData"))

csv_files <- paste0(output_dir,out_base,"-",1:3,".csv")
stanfit <- as_cmdstan_fit(files = csv_files)


# Compare stratum level beta parameters -----------------------------------



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


betas_plot = ggplot(data = beta_comp,aes(x = beta_True,y = mean))+
  geom_point(aes(colour = k))+
  scale_colour_viridis_c()+
  geom_errorbar(aes(ymin = Q_025,ymax = Q_975),width = 0,alpha = 0.2)+
  geom_abline(slope = 1, intercept = 0)+ 
  facet_wrap(~Stratum_Factored,nrow = ceiling(sqrt(stan_data$nstrata)),
             ncol = ceiling(sqrt(stan_data$nstrata)))





# compare stratum level smoothed indices ----------------------------------

true_inds <- balanced %>% 
  select(Stratum,Stratum_Factored,Smooth,Year,True_log_traj,True_strata_effects) %>% 
  distinct() %>% 
  mutate(True_nsmooth = exp(Smooth + True_strata_effects)) %>% 
  arrange(Stratum_Factored)



nsmooth_est <- posterior_samples(stanfit,
                               parm = "nsmooth",
                               dims = c("Stratum_Factored","Year_Index")) %>% 
  posterior_sums(.,
                 dims = c("Stratum_Factored","Year_Index")) %>% 
  mutate(Year = Year_Index+min(balanced$Year)-1)


nsmooth_comp <- nsmooth_est %>% 
  left_join(true_inds,by = c("Stratum_Factored","Year"))



nsmooth_plot = ggplot(data = nsmooth_comp,aes(x = True_nsmooth,
                                              y = mean))+
  geom_point(aes(colour = Year))+
  scale_colour_viridis_c()+
  geom_errorbar(aes(ymin = Q_025,ymax = Q_975),width = 0,alpha = 0.2)+
  geom_abline(slope = 1, intercept = 0)+ 
  facet_wrap(~Stratum_Factored,nrow = ceiling(sqrt(stan_data$nstrata)),
             ncol = ceiling(sqrt(stan_data$nstrata)))




true_inds <- true_inds %>% 
  mutate(version = "True")

nsmooth_comp2 <- nsmooth_est %>% 
  select(mean,Q_025,Q_975,Stratum_Factored,Year) %>% 
  rename(True_nsmooth = mean) %>% 
  mutate(version = "Estimated") %>% 
  bind_rows(.,true_inds)


nsmooth_plot2 = ggplot(data = nsmooth_comp2,aes(y = True_nsmooth,
                                              x = Year))+
  geom_ribbon(aes(ymin = Q_025,ymax = Q_975,fill = version),alpha = 0.2)+
  geom_line(aes(colour = version))+
  scale_colour_viridis_d(aesthetics = c("colour","fill"))+
  facet_wrap(~Stratum_Factored,nrow = ceiling(sqrt(stan_data$nstrata)),
             ncol = ceiling(sqrt(stan_data$nstrata)),
             scales = "free_y")




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
  mutate(version = "TRUE")

SMOOTH_est <- posterior_samples(stanfit,
                              parm = "SMOOTH_pred",
                              dims = c("y")) %>% 
  posterior_sums(.,
                 dims = c("y")) 



SMOOTH_comp <- SMOOTH_est %>% 
  select(mean,Q_025,Q_975,y) %>% 
  rename(True_SMOOTH = mean) %>% 
  mutate(version = "Estimated") %>% 
  bind_rows(.,true_SMOOTH)


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


}

}

