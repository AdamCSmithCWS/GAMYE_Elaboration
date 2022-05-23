## Summary of the trend estimates from prior simulations
setwd("C:/GitHub/GAMYE_Elaboration")

library(tidyverse)
library(patchwork)


bbs_indices_usgs <- read.csv("data/Index_best_1966-2019_core_best.csv",
                             colClasses = c("integer",
                                            "character",
                                            "integer",
                                            "numeric",
                                            "numeric",
                                            "numeric"))

bbs_inds <- bbs_indices_usgs %>% 
  filter(Region == "SU1" |
           !grepl(x = Region, pattern = "[[:digit:]]")) #just hte continental estimates

# function to calculate a %/year trend from a count-scale trajectory
trs <- function(y1,y2,ny){
  tt <- (((y2/y1)^(1/ny))-1)*100
}

miny = min(bbs_inds$Year)
maxy = max(bbs_inds$Year)
bbs_trends <- NULL

for(tl in c(2,6,11,21,54)){ #estimating all possible 1-year, 2-year, 5-year, 10-year, and 20-year trends, with no uncertainty, just the point estimates based on the comparison of posterior means fo annual indices
  ny = tl-1
  yrs1 <- seq(miny,(maxy-ny),by = 1)
  yrs2 <- yrs1+ny
  for(j in 1:length(yrs1)){
    y2 <- yrs2[j]
    y1 <- yrs1[j]
    
    nyh2 <- paste0("Y",y2)
    nyh1 <- paste0("Y",y1)
    
    tmp <- bbs_inds %>% 
      filter(Year %in% c(y1,y2)) %>% 
      select(AOU,Index,Year,Region) %>% 
      pivot_wider(.,names_from = Year,
                  values_from = Index,
                  names_prefix = "Y") %>%
      rename_with(.,~gsub(pattern = nyh2,replacement = "YE", .x)) %>% 
      rename_with(.,~gsub(pattern = nyh1,replacement = "YS", .x)) %>% 
      drop_na() %>% 
      group_by(AOU,Region) %>% 
      summarise(trend = trs(YS,YE,ny),
                .groups = "keep")%>% 
      mutate(first_year = y1,
             last_year = y2,
             nyears = ny,
             abs_trend = abs(trend),
             t_years = paste0(ny,"-year trends"))
    
    bbs_trends <- bind_rows(bbs_trends,tmp)
  }
}

t_quants <- bbs_trends %>% 
  group_by(t_years,Region) %>% 
  summarise(x99 = quantile(abs_trend,0.99),
            x995 = quantile(abs_trend,0.995))

bbs_trends <- bbs_trends %>% 
  mutate(t_years = factor(t_years,
                          levels = c("1-year trends",
                                     "5-year trends",
                                     "10-year trends",
                                     "20-year trends",
                                     "53-year trends"),
                          ordered = TRUE))

bbs_continental_trends <- bbs_trends %>% 
  filter(Region == "SU1")
bbs_politic_trends <- bbs_trends %>% 
  filter(Region != "SU1")


freq_brks <- c(0,seq(1,100,1))


realised_bbs_politic_freq <- ggplot(data = bbs_politic_trends,
                                 aes(abs_trend,after_stat(density)))+
  geom_freqpoly(breaks = freq_brks,center = 0)+
  xlab("Absolute value of BBS state/province trends USGS models (1966-2019)")+
  ylab("")+
  theme_bw()+
  coord_cartesian(ylim = c(0,0.7),
                  xlim = c(0,40))+
  facet_wrap(vars(t_years),
             nrow = 1,
             ncol = 5)
print(realised_bbs_politic_freq)



realised_bbs_sw_freq <- ggplot(data = bbs_continental_trends,
                                 aes(abs_trend,after_stat(density)))+
  geom_freqpoly(breaks = freq_brks,center = 0)+
  xlab("Absolute value of BBS survey-wide trends USGS models (1966-2019)")+
  ylab("")+
  theme_bw()+
  coord_cartesian(ylim = c(0,0.7),
                  xlim = c(0,40))+
  facet_wrap(vars(t_years),
             nrow = 1,
             ncol = 5)
print(realised_bbs_sw_freq)



bbs_sd_trends <- bbs_politic_trends %>% 
  group_by(AOU,t_years) %>% 
  summarise(sd_trends = sd(trend),
            min_trend = min(trend),
            max_trend = max(trend),
            q5_trend = quantile(trend,0.05),
            q95_trend = quantile(trend,0.95),
            .groups = "keep") %>% 
  filter(is.finite(sd_trends))

realised_bbs_sd <- ggplot(data = bbs_sd_trends,
                                    aes(sd_trends,after_stat(density)))+
  geom_freqpoly(breaks = freq_brks,center = 0)+
  xlab("SD (by species) of BBS state/province trends USGS models (1966-2019)")+
  ylab("")+
  theme_bw()+
  coord_cartesian(ylim = c(0,0.7),
                  xlim = c(0,40))+
  facet_wrap(vars(t_years),
             nrow = 1,
             ncol = 5)
print(realised_bbs_sd)


tb_sims <- data.frame(model = c(rep("GAMYE",4),
                                rep("Difference",3)),
                      spatial = c(TRUE,FALSE,FALSE,TRUE,
                                  TRUE,FALSE,FALSE),
                      hierarchical = c(TRUE,TRUE,TRUE,TRUE,
                                       TRUE,TRUE,FALSE),
                      fl = paste0(c("Hier_prior_sim_summary",
                                    "Hier_Non_Spatial_prior_sim_summary",
                                    "Hier_Non_Spatial_with_YE_prior_sim_summary",
                                    "Hier_Spatial_with_YE_prior_sim_summary",
                                    "Hier_Spatial_Difference_prior_sim_summary",
                                    "Hier_Non_Spatial_Difference_prior_sim_summary",
                                    "Non_Hierarchical_Difference_prior_sim_summary"),".RData"),
                      prior_sel = c(1,1,10,20,
                                    0.1,0.1,0.3))


# Full summary and plotting of all priors ---------------------------------




pdf(file = "Figures/Prior_trend_comparisons.pdf",
    width = 11,
    height = 8.5)

figso <- vector(mode = "list",length = nrow(tb_sims))
for(i in 1:nrow(tb_sims)){
  figso[[i]] <- vector(mode = "list",7)
M = tb_sims[i,"model"]
spat = tb_sims[i,"spatial"]
hier = tb_sims[i,"hierarchical"]
fl = tb_sims[i,"fl"]

load(paste0("output/",fl))

if(M == "GAMYE"){
  trends_out <- trends_out %>% 
    filter(prior_scale != 0.5)
}
modl <- paste("Non Hierarchical",M)
if(spat & hier){
  modl <- paste("Spatial Hierarchical",M)
}
if(!spat & hier){
  modl <- paste("Non Spatial Hierarchical",M)
}
trends_abs <- trends_out %>% 
  filter(!is.na(Stratum_Factored)) %>% 
  mutate(abs_trend = abs(trend),
         scale_factor = factor(prior_scale,ordered = TRUE),
         t_years = factor(paste0(nyears,"-year trends"),
                          levels = c("1-year trends",
                                     "5-year trends",
                                     "10-year trends",
                                     "20-year trends",
                                     "53-year trends"),
                          ordered = TRUE))

trends_abs1a <- trends_abs %>% 
  mutate(distribution_factor = factor(distribution,
                                      ordered = TRUE)) %>% 
  filter(distribution_factor %in% c("t3"))

if(length(unique(trends_abs1a$param)) > 1){
  
  if(spat){
    modl <- paste("Spatial Hierarchical",M,"with YE")
  }else{
    modl <- paste("Non Spatial Hierarchical",M,"with YE")
    
  }
  
  trends_abs1 <- trends_abs1a %>% 
    filter(param == "nsmooth")
  
  
  comp_plot_strat <- realised_bbs_politic_freq +
    geom_freqpoly(data = trends_abs1,
                  aes(abs_trend,after_stat(density),
                      colour = scale_factor),
                  breaks = freq_brks,
                  center = 0,
                  alpha = 0.5)+
    scale_colour_viridis_d(begin = 0.4,end = 0.9,
                           guide_legend(title = "Prior scale"))+
    labs(title = paste("Regional trends SMOOTH",modl))
  
  
  
  trends_sd <- trends_abs1 %>% 
    group_by(.draw,t_years,distribution_factor,scale_factor) %>% 
    summarise(sd_trends = sd(trend,na.rm = T),
              min_trend = min(trend,na.rm = T),
              max_trend = max(trend,na.rm = T),
              q5_trend = quantile(trend,0.05,na.rm = T),
              q95_trend = quantile(trend,0.95,na.rm = T),
              .groups = "keep") %>% 
    filter(is.finite(sd_trends))
  
  comp_plot_sd <- realised_bbs_sd +
    geom_freqpoly(data = trends_sd,
                  aes(sd_trends,after_stat(density),
                      colour = scale_factor),
                  breaks = freq_brks,
                  center = 0,
                  alpha = 0.5)+
    scale_colour_viridis_d(begin = 0.4,end = 0.9,
                           guide_legend(title = "Prior scale"))+
    labs(title = paste("SD trends SMOOTH",modl))
  
  
  trends_abs1 <- trends_abs1a %>% 
    filter(param == "n")
  
  
  comp_plot_strat2 <- realised_bbs_politic_freq +
    geom_freqpoly(data = trends_abs1,
                  aes(abs_trend,after_stat(density),
                      colour = scale_factor),
                  breaks = freq_brks,
                  center = 0,
                  alpha = 0.5)+
    scale_colour_viridis_d(begin = 0.4,end = 0.9,
                           guide_legend(title = "Prior scale"))+
    labs(title = paste("Regional trends FULL",modl))
  
  
  
  
  trends_sd <- trends_abs1 %>% 
    group_by(.draw,t_years,distribution_factor,scale_factor) %>% 
    summarise(sd_trends = sd(trend,na.rm = T),
              min_trend = min(trend,na.rm = T),
              max_trend = max(trend,na.rm = T),
              q5_trend = quantile(trend,0.05,na.rm = T),
              q95_trend = quantile(trend,0.95,na.rm = T),
              .groups = "keep") %>% 
    filter(is.finite(sd_trends))
  
  comp_plot_sd2 <- realised_bbs_sd +
    geom_freqpoly(data = trends_sd,
                  aes(sd_trends,after_stat(density),
                      colour = scale_factor),
                  breaks = freq_brks,
                  center = 0,
                  alpha = 0.5)+
    scale_colour_viridis_d(begin = 0.4,end = 0.9,
                           guide_legend(title = "Prior scale"))+
    labs(title = paste("SD trends FULL",modl))
  
  
}else{
  trends_abs1 <- trends_abs1a
  
  comp_plot_strat <- realised_bbs_politic_freq +
    geom_freqpoly(data = trends_abs1,
                  aes(abs_trend,after_stat(density),
                      colour = scale_factor),
                  breaks = freq_brks,
                  center = 0,
                  alpha = 0.5)+
    scale_colour_viridis_d(begin = 0.4,end = 0.9,
                           guide_legend(title = "Prior scale"))+
    labs(title = paste("Regional trends",modl))
  
  
  comp_plot_strat2 <- NA
  
  trends_sd <- trends_abs1 %>% 
    group_by(.draw,t_years,distribution_factor,scale_factor) %>% 
    summarise(sd_trends = sd(trend,na.rm = T),
              min_trend = min(trend,na.rm = T),
              max_trend = max(trend,na.rm = T),
              q5_trend = quantile(trend,0.05,na.rm = T),
              q95_trend = quantile(trend,0.95,na.rm = T),
              .groups = "keep") %>% 
    filter(is.finite(sd_trends))
  
  comp_plot_sd <- realised_bbs_sd +
    geom_freqpoly(data = trends_sd,
                  aes(sd_trends,after_stat(density),
                      colour = scale_factor),
                  breaks = freq_brks,
                  center = 0,
                  alpha = 0.5)+
    scale_colour_viridis_d(begin = 0.4,end = 0.9,
                           guide_legend(title = "Prior scale"))+
    labs(title = paste("SD trends",modl))
  
  comp_plot_sd2 <- NA
}





TRENDS_abs <- trends_out %>% 
  filter(is.na(Stratum_Factored),
         param %in% c("N","NSMOOTH")) %>% 
  mutate(abs_trend = abs(trend),
         scale_factor = factor(prior_scale,ordered = TRUE),
         t_years = factor(paste0(nyears,"-year trends"),
                          levels = c("1-year trends",
                                     "5-year trends",
                                     "10-year trends",
                                     "20-year trends",
                                     "53-year trends"),
                          ordered = TRUE))


  
TRENDS_abs1a <- TRENDS_abs %>% 
  mutate(distribution_factor = factor(distribution,
                                      ordered = TRUE)) %>% 
  filter(distribution_factor %in% c("t3"))

if(length(unique(TRENDS_abs1a$param)) > 1){
  
  TRENDS_abs1 <- TRENDS_abs1a %>% 
    filter(param == "NSMOOTH")
  
comp_plot_sw <- realised_bbs_sw_freq +
  geom_freqpoly(data = TRENDS_abs1,
                aes(abs_trend,after_stat(density),
                    colour = scale_factor),
                breaks = freq_brks,
                center = 0,
                alpha = 0.5)+
  scale_colour_viridis_d(begin = 0.4,end = 0.9,
                         guide_legend(title = "Prior scale"))+
  labs(title = paste("Survey-wide trends SMOOTH",modl))


TRENDS_abs1 <- TRENDS_abs1a %>% 
  filter(param == "N")

comp_plot_sw2 <- realised_bbs_sw_freq +
  geom_freqpoly(data = TRENDS_abs1,
                aes(abs_trend,after_stat(density),
                    colour = scale_factor),
                breaks = freq_brks,
                center = 0,
                alpha = 0.5)+
  scale_colour_viridis_d(begin = 0.4,end = 0.9,
                         guide_legend(title = "Prior scale"))+
  labs(title = paste("Survey-wide trends FULL",modl))


}else{
  
  TRENDS_abs1 <- TRENDS_abs1a 
  
  comp_plot_sw <- realised_bbs_sw_freq +
    geom_freqpoly(data = TRENDS_abs1,
                  aes(abs_trend,after_stat(density),
                      colour = scale_factor),
                  breaks = freq_brks,
                  center = 0,
                  alpha = 0.5)+
    scale_colour_viridis_d(begin = 0.4,end = 0.9,
                           guide_legend(title = "Prior scale"))+
    labs(title = paste("Survey-wide trends",modl))
 
  

  comp_plot_sw2 <- NA
  
}

TRENDS_abs_comp <- trends_out %>% 
  filter(is.na(Stratum_Factored),
         param %in% c("NComp","NSmoothComp")) %>% 
  mutate(abs_trend = abs(trend),
         scale_factor = factor(prior_scale,ordered = TRUE),
         t_years = factor(paste0(nyears,"-year trends"),
                          levels = c("1-year trends",
                                     "5-year trends",
                                     "10-year trends",
                                     "20-year trends",
                                     "53-year trends"),
                          ordered = TRUE))

TRENDS_abs_comp1a <- TRENDS_abs_comp %>% 
  mutate(distribution_factor = factor(distribution,
                                      ordered = TRUE)) %>% 
  filter(distribution_factor %in% c("t3"))

if(length(unique(TRENDS_abs_comp1a$param)) > 1){
  TRENDS_abs_comp1 <- TRENDS_abs_comp1a %>% 
    filter(param == "NSmoothComp")
  
  comp_plot_sw_comp <- realised_bbs_sw_freq +
    geom_freqpoly(data = TRENDS_abs_comp1,
                  aes(abs_trend,after_stat(density),
                      colour = scale_factor),
                  breaks = freq_brks,
                  center = 0,
                  alpha = 0.5)+
    scale_colour_viridis_d(begin = 0.4,end = 0.9,
                           guide_legend(title = "Prior scale"))+
    labs(title = paste("Composite Survey-wide trends SMOOTH",modl))
  
  TRENDS_abs_comp1 <- TRENDS_abs_comp1a %>% 
    filter(param == "NComp")
  
  comp_plot_sw_comp2 <- realised_bbs_sw_freq +
    geom_freqpoly(data = TRENDS_abs_comp1,
                  aes(abs_trend,after_stat(density),
                      colour = scale_factor),
                  breaks = freq_brks,
                  center = 0,
                  alpha = 0.5)+
    scale_colour_viridis_d(begin = 0.4,end = 0.9,
                           guide_legend(title = "Prior scale"))+
    labs(title = paste("Composite Survey-wide trends FULL",modl))
  
  
}else{
  
  TRENDS_abs_comp1 <- TRENDS_abs_comp1a
  
comp_plot_sw_comp <- realised_bbs_sw_freq +
  geom_freqpoly(data = TRENDS_abs_comp1,
                aes(abs_trend,after_stat(density),
                    colour = scale_factor),
                breaks = freq_brks,
                center = 0,
                alpha = 0.5)+
  scale_colour_viridis_d(begin = 0.4,end = 0.9,
                         guide_legend(title = "Prior scale"))+
  labs(title = paste("Composite Survey-wide trends",modl))

comp_plot_sw_comp2 <- NA

}

figso[[i]][[1]] <- comp_plot_strat
figso[[i]][[2]] <- comp_plot_sw
figso[[i]][[3]] <- comp_plot_sw_comp
figso[[i]][[4]] <- comp_plot_sd
figso[[i]][[5]] <- comp_plot_strat2
figso[[i]][[6]] <- comp_plot_sw_comp2
figso[[i]][[7]] <- comp_plot_sd2


fnp <- (comp_plot_strat/comp_plot_sd) +
  plot_layout(guides = "collect")

print(fnp)

if(!is.na(comp_plot_strat2[1])){
  fnp <- (comp_plot_strat2/comp_plot_sd2) +
    plot_layout(guides = "collect")
  
  print(fnp)
}
fnp <- (comp_plot_sw/comp_plot_sw_comp) +
  plot_layout(guides = "collect")

print(fnp)

if(!is.na(comp_plot_sw_comp2[1])){
  
  fnp <- (comp_plot_sw_comp2/comp_plot_sw_comp) +
    plot_layout(guides = "collect")
  
  print(fnp)
  
}

print(modl)
}
dev.off()



to_save <- c("tb_sims",
             "realised_bbs_sd",
             "bbs_sd_trends",
             "realised_bbs_sw_freq",
             "bbs_continental_trends",
             "realised_bbs_politic_freq",
             "bbs_politic_trends",
             "trends_abs_all",
             "freq_brks")

to_rem <- ls()[-which(ls() %in% to_save)]
rm(list = to_rem)




# Plotting selected priors ------------------------------------------------


trends_abs_all <- NULL

for(i in 1:nrow(tb_sims)){
  M = tb_sims[i,"model"]
  spat = tb_sims[i,"spatial"]
  hier = tb_sims[i,"hierarchical"]
  fl = tb_sims[i,"fl"]
  
  load(paste0("output/",fl))
  
  if(M == "GAMYE"){
    trends_out <- trends_out %>% 
      filter(prior_scale != 0.5)
  }
  modl <- paste("Non Hierarchical",M)
  if(spat & hier){
    modl <- paste("Spatial Hierarchical",M)
  }
  if(!spat & hier){
    modl <- paste("Non Spatial Hierarchical",M)
  }
  if(i %in% c(3,4)){
    if(spat){
      modl <- paste("Spatial Hierarchical",M,"with YE")
    }else{
      modl <- paste("Non Spatial Hierarchical",M,"with YE")
      
    }
    }
  
  
  trends_abs <- trends_out %>% 
    mutate(abs_trend = abs(trend),
           t_years = factor(paste0(nyears,"-year trends"),
                            levels = c("1-year trends",
                                       "5-year trends",
                                       "10-year trends",
                                       "20-year trends",
                                       "53-year trends"),
                            ordered = TRUE),
           Model = M,
           model = modl,
           spatial = spat,
           hierarchical = hier) %>% 
    filter(distribution %in% c("t3"))
  
 trends_abs_all <- bind_rows(trends_abs_all,trends_abs)
  
  print(i)
}


trends_abs_all <- trends_abs_all %>% 
  mutate(scale_factor = factor(prior_scale,ordered = TRUE))

save(list = "trends_abs_all",
     file = "output/allpriorsimulationtrends.RData")




# Select the best priors by model -----------------------------------------


pdf(file = "Figures/Selected_prior_distributions.pdf",
    height = 8.5,
    width = 11)
# 
# mod_pal <- scale_colour_viridis_d(begin = 0.2,end = 0.9,
#                                     guide_legend(title = "Model"))
# 
#     
# mod_pal <- scale_color_brewer(type = "div",
#                               palette = "PuOr",
#                               guide = "legend",
#                               guide_legend(title = "Model"))
  
 my_cols <- RColorBrewer::brewer.pal(7,"Paired")
 names(my_cols) <- c("Spatial Hierarchical GAMYE",
                     "Spatial Hierarchical GAMYE with YE",
                     "Non Spatial Hierarchical GAMYE",
                     "Non Spatial Hierarchical GAMYE with YE",
                     "Spatial Hierarchical Difference",
                     "Non Spatial Hierarchical Difference",
                     "Non Hierarchical Difference")
mod_pal <- scale_color_manual(values = my_cols,
                              guide = "legend",
                              guide_legend(title = "Model"))


  sel_prs <- tb_sims %>% 
  select(model,spatial,hierarchical,prior_sel) %>% 
  rename(Model = model,
         prior_scale = prior_sel) 
sel_prs[,"model"] <- c("Spatial Hierarchical GAMYE",
                       "Non Spatial Hierarchical GAMYE",
                       "Non Spatial Hierarchical GAMYE with YE",
                       "Spatial Hierarchical GAMYE with YE",
                       "Spatial Hierarchical Difference",
                       "Non Spatial Hierarchical Difference",
                       "Non Hierarchical Difference")


comp_trends <- trends_abs_all %>% 
  filter(param %in% c("NComp","NSmoothComp"))%>% 
  inner_join(.,sel_prs,by = c("Model",
                              "model",
                              "spatial",
                              "hierarchical",
                              "prior_scale")) 

plot_composite <- realised_bbs_sw_freq +
  geom_freqpoly(data = comp_trends,
                aes(abs_trend,after_stat(density),
                    colour = model),
                breaks = freq_brks,
                center = 0,
                alpha = 0.8)+
  mod_pal+
  facet_grid(cols = vars(t_years),
             rows = vars(param))+
  labs(title = paste("Composite Survey-wide trends"))

print(plot_composite)


Hyper_trends <- trends_abs_all %>% 
  filter(param %in% c("N","NSMOOTH")) %>% 
  inner_join(.,sel_prs,by = c("Model",
                              "model",
                              "spatial",
                              "hierarchical",
                              "prior_scale")) 
plot_hyper <- realised_bbs_sw_freq +
  geom_freqpoly(data = Hyper_trends,
                aes(abs_trend,after_stat(density),
                    colour = model),
                breaks = freq_brks,
                center = 0,
                alpha = 0.8)+
  mod_pal+
  facet_grid(cols = vars(t_years),
             rows = vars(param))+
  labs(title = paste("Hyper-parameter Survey-wide trends"))

print(plot_hyper)



# stratum level trends ----------------------------------------------------
strata_trendsn <- trends_abs_all %>% 
  filter(param %in% c("n"))  %>% 
  inner_join(.,sel_prs,by = c("Model",
                              "model",
                              "spatial",
                              "hierarchical",
                              "prior_scale"))
 
strata_trendsm <- trends_abs_all %>% 
  filter(param %in% c("nsmooth"))  %>% 
  inner_join(.,sel_prs,by = c("Model",
                              "model",
                              "spatial",
                              "hierarchical",
                              "prior_scale"))

strata_trends <- bind_rows(strata_trendsn,
                           strata_trendsm)

plot_strat_trends <- realised_bbs_politic_freq +
  geom_freqpoly(data = strata_trends,
                aes(abs_trend,after_stat(density),
                    colour = model),
                breaks = freq_brks,
                center = 0,
                alpha = 0.8)+
  mod_pal+
  facet_grid(cols = vars(t_years),
             rows = vars(param))+
  labs(title = paste("Stratum trends"))
  
print(plot_strat_trends)
# SD of trends ------------------------------------------------------------

trends_abs_sel <- trends_abs_all %>% 
  inner_join(.,sel_prs,by = c("Model",
                              "model",
                              "spatial",
                              "hierarchical",
                              "prior_scale"))

trends_sd <- trends_abs_sel %>% 
  filter(param %in% c("nsmooth","n")) %>% 
  group_by(.draw,t_years,distribution,
           Model,model,prior_scale,param) %>% 
  summarise(sd_trends = sd(trend),
            min_trend = min(trend),
            max_trend = max(trend),
            q5_trend = quantile(trend,0.05),
            q95_trend = quantile(trend,0.95),
            .groups = "keep") %>% 
  filter(is.finite(sd_trends))

comp_plot_sd <- realised_bbs_sd +
  geom_freqpoly(data = trends_sd,
                aes(sd_trends,after_stat(density),
                    colour = model),
                breaks = freq_brks,
                center = 0,
                alpha = 0.8)+
  mod_pal+
  facet_grid(cols = vars(t_years),
             rows = vars(param))+
  labs(title = paste("SD trends"))

print(comp_plot_sd)

dev.off()



# Compare selected prior graphs for most relevant situations









