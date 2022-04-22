## Summary of the trend estimates from prior simulations



load("output/Non_Hierarchical_Difference_prior_sim_summary.RData")



trends_abs <- trends_out %>% 
  filter(!is.na(Stratum_Factored),
         distribution != "t4",
         prior_scale != 0.5) %>% 
  mutate(abs_trend = abs(trend),
         scale_factor = factor(prior_scale,ordered = TRUE),
         nyears_factor = factor(nyears,levels = c(1,5,10,20,53),ordered = TRUE))





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
