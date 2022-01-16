### prior simulation of spatial GAM sd

library(bbsBayes)
library(tidyverse)
library(cmdstanr)

setwd("C:/Users/adam_/OneDrivedelete/Documents/GitHub/GAMYE_Elaboration")


# fit model with fixed data for all parameters except the local sm --------

species <- "Pine Warbler"
species_f <- gsub(species,pattern = " ",replacement = "_")

dd <- NULL
for(rr in c(0.5,1,2,4)){
 tmp = data.frame(p = rgamma(10000,2,rr),
                  i = 1:10000,
                  pr = paste0("gamma_",rr)) %>% 
   mutate()
 dd <- bind_rows(dd,tmp)
 
}
for(rr in rev(c(0.5,1,2,4))){
  tmp = data.frame(p = abs(rnorm(10000,0,rr)),
                   i = 1:10000,
                   pr = paste0("norm_",rr))
  dd <- bind_rows(dd,tmp)
  
}

hists = ggplot(data = dd)+
  geom_histogram(aes(x = p))+
  facet_wrap(~pr,nrow = 2,ncol = 4)
print(hists)

for(pp in c("gamma","norm")){
  for(prior_scale in c(0.5,1,2,4)){
    
  tp = paste0(pp,"_rate_",prior_scale)

  #STRATA_True <- log(2)
  output_dir <- "output/"
  out_base <- paste0(species_f,"_sim_",tp,"_BBS")
  csv_files <- paste0(out_base,"-",1:3,".csv")
  
  if(pp == "gamma"){
    pnorm <- 0
  }else{
    pnorm <- 1
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



