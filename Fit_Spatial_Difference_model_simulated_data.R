### Fitting model to simulated BBS data
library(bbsBayes)
library(tidyverse)
library(cmdstanr)

setwd("C:/GitHub/GAMYE_Elaboration")

#species <- "Pine Warbler"

species = "Yellow-headed Blackbird"  

  species_f <- gsub(species,pattern = " ",replacement = "_")
 
  tp = "breakpoint_cycle"
  
MAs <- round(log(c(0.1,0.5,1,5,10,50)),2)# true mean abundances for different simulations


for(ma in MAs[c(6,1,2,3,4,5)]){  
  
  
  output_dir <- "output/"
  out_base <- paste0("sim_Spatial_Difference",tp,"_",ma,"_BBS")
  csv_files <- paste0(out_base,"-",1:3,".csv")
  
  load(paste0("Data/Simulated_data_",ma,"_",tp,"_BBS.RData"))
  
  
        if(!file.exists(paste0(output_dir,csv_files[1]))){
          
        
tmp_data = realised
  
nsites = max(routes_df$Route_Factored)
nstrata = max(strata_df$Stratum_Factored)
ncounts = nrow(tmp_data)
nyears = max(tmp_data$Year_Index)

count = tmp_data$count
strat = tmp_data$Stratum_Factored
year = tmp_data$Year_Index
site = tmp_data$Route_Factored

N_edges = neighbours$N_edges
node1 = neighbours$node1
node2 = neighbours$node2

nobservers = max(tmp_data$Observer_Factored)
observer = tmp_data$Observer_Factored

midyear = floor(nyears/2)

Iy1 = c((midyear-1):1)
Iy2 = c((midyear+1):nyears)
nIy1 = length(Iy1)
nIy2 = length(Iy2)


nroutes_strata <- routes_df %>% arrange(Stratum_Factored,
                                     Route_Factored) %>% 
  group_by(Stratum_Factored) %>% 
  summarise(nroutes = n())



nsites_strata <- as.integer(nroutes_strata$nroutes)
maxnsites_strata <- max(nsites_strata)

ste_mat <- matrix(data = 0,
                  nrow = nstrata,
                  ncol = maxnsites_strata)
for(i in 1:nstrata){
  ste_mat[i,1:nsites_strata[i]] <- routes_df[which(routes_df$Stratum_Factored == i),"Route_Factored"]
}

nonzeroweight <- rep(1,nstrata)
#nu <- 100 #fixed df for t-distributed noise

stan_data = list(#scalar indicators
                 nsites = nsites,
                 nstrata = nstrata,
                 ncounts = ncounts,
                 nyears = nyears,
                 nyears_m1 = nyears-1,
                 
                 #basic data
                 count = count,
                 strat = strat,
                 year = year,
                 site = site,
                 
                 #spatial structure
                 N_edges = N_edges,
                 node1 = node1,
                 node2 = node2,
                 
                 #Temporal Center and indexing
                 midyear = midyear,
                 Iy1 = Iy1,
                 Iy2 = Iy2,
                 nIy1 = nIy1,
                 nIy2 = nIy2,
                 
                 
                 #Observer information
                 nobservers = nobservers,
                 observer = observer,
                 
                 #Ragged array information to link sites to strata
                 nsites_strata = nsites_strata,
                 maxnsites_strata = maxnsites_strata,
                 ste_mat = ste_mat,
                 #nu = nu,
                 
                 #weights
                 nonzeroweight = nonzeroweight,
                 
                 #vector of zeros to fill midyear beta values
                 zero_betas = rep(0,nstrata)
)




# Fit model ---------------------------------------------------------------

print(paste("beginning",out_base,Sys.time()))


  mod.file = "models/difference_iCAR_sim.stan"
  
  
 

  init_def <- function(){ list(noise_raw = rnorm(ncounts,0,0.1),
                               strata_raw = rnorm(nstrata,0,0.1),
                               STRATA = 0,
                               sdstrata = runif(1,0.01,0.1),
                               #eta = 0,
                               obs_raw = rnorm(nobservers,0,0.01),
                               ste_raw = rnorm(nsites,0,0.01),
                               sdnoise = runif(1,0.01,0.2),
                               sdobs = runif(1,0.01,0.1),
                               sdste = runif(1,0.01,0.2),
                               sdbeta = runif(nyears-1,0.01,0.1),
                               sdBETA = runif(1,0.01,0.1),
                               BETA_raw = rnorm(nyears-1,0,0.1),
                               beta_raw = matrix(rnorm((nyears-1)*nstrata,0,0.01),
                                                 nrow = nstrata,
                                                 ncol = nyears-1))}
  

## compile model
model <- cmdstan_model(mod.file)


stanfit <- model$sample(
  data=stan_data,
  refresh=200,
  chains=3, iter_sampling=1500,
  iter_warmup=1500,
  parallel_chains = 3,
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

}#end ma loop

