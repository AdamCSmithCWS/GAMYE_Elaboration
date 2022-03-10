// This is a Stan implementation of the gamye model
// with iCAR component for the stratum-level intercepts and smooth parameters

// Consider moving annual index calculations outside of Stan to 
// facilitate the ragged array issues

// iCAR function, from Morris et al. 2019
// Morris, M., K. Wheeler-Martin, D. Simpson, S. J. Mooney, A. Gelman, and C. DiMaggio (2019). 
// Bayesian hierarchical spatial models: Implementing the Besag York Mollié model in stan. 
// Spatial and Spatio-temporal Epidemiology 31:100301.

// functions {
//   real icar_normal_lpdf(vector bb, int ns, int[] n1, int[] n2) {
//     return -0.5 * dot_self(bb[n1] - bb[n2])
//       + normal_lpdf(sum(bb) | 0, 0.001 * ns); //soft sum to zero constraint on bb
//  }
// }


data {
  int<lower=1> nsites;
  int<lower=1> nstrata;
  int<lower=1> ncounts;
  int<lower=1> nyears;

  int<lower=0> count[ncounts];              // count observations
  int<lower=1> strat[ncounts];               // strata indicators
  int<lower=1> year[ncounts]; // year index
  int<lower=1> site[ncounts]; // site index
  
  int<lower=1> observer[ncounts];              // observer indicators
  int<lower=1> nobservers;

 int<lower=0> nsites_strata[nstrata]; // number of sites in each stratum
 int<lower=0> maxnsites_strata; //largest value of nsites_strata

  int ste_mat[nstrata,maxnsites_strata]; //matrix identifying which sites are in each stratum
  // above is actually a ragged array, but filled with 0 values so that it works
  // but throws an error if an incorrect strata-site combination is called
 
  // // spatial neighbourhood information
  // int<lower=1> N_edges;
  // int<lower=1, upper=nstrata> node1[N_edges];  // node1[i] adjacent to node2[i]
  // int<lower=1, upper=nstrata> node2[N_edges];  // and node1[i] < node2[i]



  // data to center the abundance estimate
  int midyear; //middle year of the time-series scaled to ~(nyears/2)

}

parameters {
  vector[ncounts] noise_raw;             // over-dispersion
 
 vector[nstrata] strata_raw;
  real STRATA; 

  vector[nobservers] obs_raw;    // sd of year effects
  vector[nsites] ste_raw;   // 
  real<lower=0> sdnoise;    // sd of over-dispersion
  real<lower=0> sdobs;    // sd of observer effects
  real<lower=0> sdste;    // sd of site effects
  real<lower=0> sdbeta[nyears-1];    // sd of annual changes among strata 
  real<lower=0> sdstrata;    // sd of intercepts
  real<lower=0> sdBETA;    // sd of overall annual changes

  vector[nyears-1] BETA_raw;//_raw; 
  matrix[nstrata,nyears-1] beta_raw;         // strata level parameters

}

transformed parameters { 
  vector[ncounts] E;           // log_scale additive likelihood
  matrix[nstrata,nyears] beta;         // spatial effect slopes (0-centered deviation from continental mean slope B)

  matrix[nstrata,nyears] yeareffect;
  vector[nyears-1] BETA;
  vector[nstrata] zero_betas = 0;
  vector[nstrata] strata;
  
  strata = (sdstrata*strata_raw) + STRATA;
  beta[,midyear] = zero_betas;
  yeareffect[,midyear] = strata;

  for(t in (midyear-1):1){
    beta[,t] = (sdbeta[t] * beta_raw[,t]) + BETA[t];
    yeareffect[,t] = yeareffect[,t+1] + beta[,t];
  }
 
   for(t in (midyear+1):nyears){
    beta[,t] = (sdbeta[t-1] * beta_raw[,t-1]) + BETA[t-1];//t-1 indicators to match dimensionality
    yeareffect[,t] = yeareffect[,t-1] + beta[,t];
  }
 
   



// intercepts and slopes

  


  for(i in 1:ncounts){
    real noise = sdnoise*noise_raw[i];
    real obs = sdobs*obs_raw[observer[i]];
  real ste = sdste*ste_raw[site[i]]; // site intercepts

    E[i] =  yeareffect[strat[i],year[i]] + ste + obs + noise;
  }
  
  }
  
  
  
model {
  sdnoise ~ normal(0,0.5); //prior on scale of extra Poisson log-normal variance
  sdobs ~ normal(0,0.5); //prior on sd of observer effects
  sdste ~ std_normal(); //prior on sd of site effects
  sdyear ~ gamma(2,2); // prior on sd of yeareffects - stratum specific, and boundary-avoiding with a prior mode at 0.5 (1/2) - recommended by https://doi.org/10.1007/s11336-013-9328-2 
  sdBETA ~ student_t(3,0,1); // prior on sd of yearly changes
  sdbeta ~ student_t(3,0,1); // prior on sd of strata variation in yearly changes
  sdstrata ~ std_normal(); //prior on sd of intercept variation
  
  noise_raw ~ normal(0,1);
  
  obs_raw ~ std_normal();//observer effects
  sum(obs_raw) ~ normal(0,0.001*nobservers);

  ste_raw ~ std_normal();//site effects
  sum(ste_raw) ~ normal(0,0.001*nsites);
 

  BETA_raw ~ std_normal();// prior on fixed effect mean GAM parameters
 
  STRATA ~ std_normal();// prior on fixed effect mean intercept

for(t in 1:(years-1)){
    beta_raw[,t] ~ std_normal();
  sum(beta_raw[,t]) ~ normal(0,0.001*nstrata);
    
}

    strata_raw ~ std_normal();
  sum(strata_raw) ~ normal(0,0.001*nstrata);
  
// for(t in 1:(years-1)){
//     beta_raw[,t] ~ icar_normal(nstrata, node1, node2);
// }
//    strata_raw ~ icar_normal(nstrata, node1, node2);

  count ~ poisson_log(E); //vectorized count likelihood with log-transformation


}

 generated quantities {

   real<lower=0> n[nstrata,nyears];
   real<lower=0> retrans_noise;
   real<lower=0> retrans_obs;
  // vector[ncounts] log_lik; // alternative value to track the observervation level log-likelihood
  // potentially useful for estimating loo-diagnostics, such as looic
  
  // for(i in 1:ncounts){
  // log_lik[i] = poisson_log_lpmf(count[i] | E[i]);
  // }
  
  
retrans_noise = 0.5*(sdnoise^2);
retrans_obs = 0.5*(sdobs^2);

for(y in 1:nyears){

      for(s in 1:nstrata){

  real n_t[nsites_strata[s]];
  real retrans_yr = 0.5*(sdyear[s]^2);

        for(t in 1:nsites_strata[s]){

  real ste = sdste*ste_raw[ste_mat[s,t]]; // site intercepts


      n_t[t] = exp(ste + yeareffect[s,y] + retrans_noise + retrans_obs);

       }
        n[s,y] = mean(n_t); //mean of exponentiated predictions across sites in a stratum
        //using the mean of hte exponentiated values, instead of including the log-normal
        // retransformation factor (0.5*sdste^2), because this retransformation makes 2 questionable assumptions:
          // 1 - assumes that sites are exchangeable among strata - e.g., that sdste is equal among all strata
          // 2 - assumes that the distribution of site-effects is normal



    }
  }



 }

