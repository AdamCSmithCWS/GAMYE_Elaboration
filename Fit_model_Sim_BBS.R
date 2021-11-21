### Fitting model to simulated BBS data
library(bbsBayes)
library(tidyverse)
library(cmdstanr)

load("Simulated_data_BBS.RData")



for(smpl in c("balanced","realised")){
# GEnerate data -----------------------------------------------------------
if(smpl == "balanced"){tmp_data = balanced}
  if(smpl == "realised"){tmp_data = realised}
  
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

nknots_year = GAM_year$nknots_Year
year_basis = GAM_year$Year_basis

nobservers = max(tmp_data$Observer_Factored)
observer = tmp_data$Observer_Factored


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

nonzeroweight <- rep(0,nstrata)

stan_data = list(#scalar indicators
                 nsites = nsites,
                 nstrata = nstrata,
                 ncounts = ncounts,
                 nyears = nyears,
                 
                 #basic data
                 count = count,
                 strat = strat,
                 year = year,
                 site = site,
                 
                 #spatial structure
                 N_edges = N_edges,
                 node1 = node1,
                 node2 = node2,
                 
                 #GAM structure
                 nknots_year = nknots_year,
                 year_basis = year_basis,
                 
                 #Observer information
                 nobservers = nobservers,
                 observer = observer,
                 
                 #Ragged array information to link sites to strata
                 nsites_strata = nsites_strata,
                 maxnsites_strata = maxnsites_strata,
                 ste_mat = ste_mat,
                 
                 #weights
                 nonzeroweight = nonzeroweight
)




# Fit model ---------------------------------------------------------------

print(paste("beginning",smpl,"with",nstrata,"strata",Sys.time()))

mod.file = "models/gamye_iCAR.stan"

## compile model
model <- cmdstan_model(mod.file)


# Initial Values ----------------------------------------------------------


init_def <- function(){ list(noise_raw = rnorm(ncounts,0,0.1),
                             strata_raw = rnorm(nstrata,0,0.1),
                             STRATA = 0,
                             #eta = 0,
                             yeareffect_raw = matrix(rnorm(nstrata*nyears,0,0.1),nrow = nstrata,ncol = nyears),
                             obs_raw = rnorm(nobservers,0,0.1),
                             ste_raw = rnorm(nsites,0,0.1),
                             sdnoise = 0.2,
                             sdobs = 0.1,
                             sdste = 0.2,
                             sdbeta = runif(nknots_year,0.01,0.1),
                             sdBETA = 0.1,
                             sdyear = runif(nstrata,0.01,0.1),
                             BETA_raw = rnorm(nknots_year,0,0.1),
                             beta_raw = matrix(rnorm(nknots_year*nstrata,0,0.01),nrow = nstrata,ncol = nknots_year))}

output_dir <- "output/"
out_base <- paste0(smpl,"_BBS")

stanfit <- model$sample(
  data=stan_data,
  refresh=200,
  chains=3, iter_sampling=1000,
  iter_warmup=1000,
  parallel_chains = 3,
  #pars = parms,
  adapt_delta = 0.8,
  max_treedepth = 14,
  seed = 123,
  init = init_def,
  output_dir = output_dir,
  output_basename = out_base)


csv_files <- paste0(out_base,"-",1:3,".csv")


save(list = c("stanfit","stan_data","csv_files"),
     file = paste0(output_dir,"/",out_base,"_gamye_iCAR.RData"))





}#end smpl
