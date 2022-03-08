// This is a Stan implementation of the gamye model
// with iCAR component for the stratum-level intercepts and smooth parameters

// Consider moving annual index calculations outside of Stan to 
// facilitate the ragged array issues


data {
  int<lower=1> nsites;
  int<lower=1> nstrata;
  int<lower=1> ncounts;
  int<lower=1> nyears;

  int<lower=0> count[ncounts];              // count observations
  int<lower=1> strat[ncounts];               // strata indicators
  int<lower=1> year[ncounts]; // year index
  int<lower=1> site[ncounts]; // site index
  
 int<lower=0> nsites_strata[nstrata]; // number of sites in each stratum
 int<lower=0> maxnsites_strata; //largest value of nsites_strata

  int ste_mat[nstrata,maxnsites_strata]; //matrix identifying which sites are in each stratum
  // above is actually a ragged array, but filled with 0 values so that it works
  // but throws an error if an incorrect strata-site combination is called
 

  // CBC effort values - party_hours scaled to the mean across all surveys (party_hours/mean(party_hours))
  vector[ncounts] hours;

  // data for spline s(year)
  int<lower=1> nknots_year;  // number of knots in the basis function for year
  matrix[nyears, nknots_year] year_basis; // basis function matrix
  
  // circle inclusion scaling value
  vector[nstrata] nonzeroweight;

}

parameters {
  vector[ncounts] noise_raw;             // over-dispersion
 
  vector[nstrata] strata_raw;
  real STRATA; 
  matrix[nstrata,nyears] yeareffect_raw;

  vector[nsites] ste_raw;   // 
  real<lower=0> sdnoise;    // sd of over-dispersion
  real<lower=0> sdste;    // sd of site effects
  real<lower=0> sdbeta[nstrata];    // sd of GAM coefficients among strata 
  real<lower=0> sdstrata;    // sd of intercepts
  real<lower=0> sdBETA;    // sd of GAM coefficients
  real<lower=0> sdyear[nstrata];    // sd of year effects
  
  real B; //mean of the effort slope
  real P; //mean of the effort exponent
  vector[nstrata] b_raw; //stratum effort slopes
  vector[nstrata] p_raw; //stratum effort exponents
  real<lower=0> sdb;    // sd of effort slopes
  real<lower=0> sdp;    // sd of effort exponents
  
  vector[nknots_year] BETA_raw;//_raw; 
  matrix[nstrata,nknots_year] beta_raw;         // GAM strata level parameters

}

transformed parameters { 
  vector[ncounts] E;           // log_scale additive likelihood
  matrix[nstrata,nknots_year] beta;         // spatial effect slopes (0-centered deviation from continental mean slope B)
  matrix[nyears,nstrata] smooth_pred;
  vector[nyears] SMOOTH_pred;  

  matrix[nstrata,nyears] yeareffect;
  vector[nknots_year] BETA;
  
  
  
  BETA = sdBETA*BETA_raw;
  
  
  for(s in 1:nstrata){
    beta[s,] = (sdbeta[s] * beta_raw[s,]) + transpose(BETA);
  } 
  
  for(s in 1:nstrata){
     smooth_pred[,s] = year_basis * transpose(beta[s,]);
}

for(s in 1:nstrata){
    yeareffect[s,] = sdyear[s]*yeareffect_raw[s,];

}

// intercepts and slopes

  


  for(i in 1:ncounts){
    real noise = sdnoise*noise_raw[i];
    real strata = (sdstrata*strata_raw[strat[i]]) + STRATA;
    real ste = sdste*ste_raw[site[i]]; // site intercepts
    real b = sdb * b_raw[strat[i]] + B;  // effort slope for stratum i
    real p = sdp * p_raw[strat[i]] + P;  //effort coefficient for stratum i
    real effort_effect = (b*((hours[i]^p)-1))/p; //effort correction for count i
    
    E[i] =  smooth_pred[year[i],strat[i]] + strata + yeareffect[strat[i],year[i]] + ste + effort_effect + noise;
  }
  
  }
  
  
  
model {
  sdnoise ~ normal(0,0.5); //prior on scale of extra Poisson log-normal variance
  noise_raw ~ normal(0,1);
  sdste ~ std_normal(); //prior on sd of site effects
  sdyear ~ gamma(2,2); // prior on sd of yeareffects - stratum specific, and boundary-avoiding with a prior mode at 0.5 (1/2) - recommended by https://doi.org/10.1007/s11336-013-9328-2 
  sdBETA ~ student_t(3,0,2); // prior on sd of GAM parameters
  sdbeta ~ student_t(3,0,2); // prior on sd of GAM parameters
  
  b_raw ~ std_normal();//effort slopes by stratum
  sum(b_raw) ~ normal(0,0.001*nstrata);
  sdb ~ normal(0,0.5); //prior on scale of effort slopes
  
  p_raw ~ std_normal();//effort exponents by stratum
  sum(p_raw) ~ normal(0,0.001*nstrata);
  sdp ~ normal(0,0.5); //prior on scale of effort exponents
  
  P ~ std_normal();//effort mean exponent
  B ~ std_normal();//effort mean slope
  

  ste_raw ~ std_normal();//site effects
  sum(ste_raw) ~ normal(0,0.001*nsites);
 
 for(s in 1:nstrata){

  yeareffect_raw[s,] ~ std_normal();
  //soft sum to zero constraint on year effects within a stratum
  sum(yeareffect_raw[s,]) ~ normal(0,0.001*nyears);
  
 }
  
  BETA_raw ~ std_normal();// prior on fixed effect mean GAM parameters
  //sum to zero constraint
  // not necessary because built into the basis function
  //sum(BETA_raw) ~ normal(0,0.001*nknots_year);
  
  STRATA ~ std_normal();// prior on fixed effect mean intercept

  //spatial iCAR intercepts and gam parameters by strata
  sdstrata ~ std_normal(); //prior on sd of intercept variation

for(s in 1:nstrata){
    beta_raw[s,] ~ normal(0,1);
}
   strata_raw ~ normal(0,1);
    sum(strata_raw) ~ normal(0,0.001*nstrata);
 
  count ~ poisson_log(E); //vectorized count likelihood with log-transformation


}

 generated quantities {

   real<lower=0> n[nstrata,nyears];
   real<lower=0> nsmooth[nstrata,nyears];
   real<lower=0> retrans_noise;
  // vector[ncounts] log_lik; // alternative value to track the observervation level log-likelihood
  // potentially useful for estimating loo-diagnostics, such as looic
  
  // for(i in 1:ncounts){
  // log_lik[i] = poisson_log_lpmf(count[i] | E[i]);
  // }
  
  
retrans_noise = 0.5*(sdnoise^2);

for(y in 1:nyears){

      for(s in 1:nstrata){

  real n_t[nsites_strata[s]];
  real nsmooth_t[nsites_strata[s]];
  real retrans_yr = 0.5*(sdyear[s]^2);
  real strata = (sdstrata*strata_raw[s]) + STRATA;
  
        for(t in 1:nsites_strata[s]){

  real ste = sdste*ste_raw[ste_mat[s,t]]; // site intercepts


      n_t[t] = exp(strata+ smooth_pred[y,s] + ste + yeareffect[s,y] + retrans_noise);
      nsmooth_t[t] = exp(strata + smooth_pred[y,s] + ste + retrans_yr + retrans_noise);
        }
        n[s,y] = mean(n_t) * nonzeroweight[s]; //mean of exponentiated predictions across sites in a stratum
        nsmooth[s,y] = mean(nsmooth_t) * nonzeroweight[s]; //mean of exponentiated predictions across sites in a stratum
        //using the mean of hte exponentiated values, instead of including the log-normal
        // retransformation factor (0.5*sdste^2), because this retransformation makes 2 questionable assumptions:
          // 1 - assumes that sites are exchangeable among strata - e.g., that sdste is equal among all strata
          // 2 - assumes that the distribution of site-effects is normal



    }
  }



 }

