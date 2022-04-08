// This is a full hierarchical GAM time-series, with spatial gam parameters 
// as well as a gam-based Seasonal adjustment and sruvey-wide random year-effects

// iCAR function, from Morris et al. 2019
// Morris, M., K. Wheeler-Martin, D. Simpson, S. J. Mooney, A. Gelman, and C. DiMaggio (2019). 
// Bayesian hierarchical spatial models: Implementing the Besag York Mollié model in stan. 
// Spatial and Spatio-temporal Epidemiology 31:100301.

functions {
  real icar_normal_lpdf(vector bb, int ns, int[] n1, int[] n2) {
    return -0.5 * dot_self(bb[n1] - bb[n2])
      + normal_lpdf(sum(bb) | 0, 0.001 * ns); //soft sum to zero constraint on bb
 }
}



data {
  int<lower=1> nsites;
  int<lower=1> nstrata;
  int<lower=1> ncounts;
  int<lower=1> nyears;

  int<lower=0> count[ncounts];              // count observations
  int<lower=1> strat[ncounts];               // strata indicators
  int<lower=1> year[ncounts]; // year index
  int<lower=1> site[ncounts]; // site index
 
 
  // spatial neighbourhood information
  int<lower=1> N_edges;
  int<lower=1, upper=nstrata> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=nstrata> node2[N_edges];  // and node1[i] < node2[i]

  
  //indexes for re-scaling predicted counts within strata based on site-level intercepts
 int<lower=0> nsites_strata[nstrata]; // number of sites in each stratum
 int<lower=0> maxnsites_strata; //largest value of nsites_strata
 int ste_mat[nstrata,maxnsites_strata]; //matrix identifying which sites are in each stratum
  // above is actually a ragged array, but filled with 0 values so that it works
  // but throws an error if an incorrect strata-site combination is called
 
   // data for spline s(year)
  int<lower=1> nknots_year;  // number of knots in the basis function for year
  matrix[nyears, nknots_year] year_basis; // basis function matrix
   
  // data for spline s(date)
  int<lower=1> ndays;  // number of days in the basis function for season
  int<lower=1> nknots_season;  // number of knots in the basis function for season
  matrix[ndays, nknots_season] season_basis; // basis function matrix
  int<lower=1> date[ncounts];  // day indicator in the season
  int<lower=1> seas_strat[ncounts];
 // indices for identifying number of seasons for each strata
  int<lower=1> seasons[nstrata,2]; //matrix of which strata have season pattern 1 or 2
 
 
}

parameters {
  vector[ncounts] noise_raw;             // over-dispersion

  vector[nsites] ste_raw;             // site intercepts
  real STRATA; // overall intercept
  
  vector[nyears] yeareffect_raw;             // continental year-effects - vector only
  
  real<lower=0> sdnoise;    // sd of over-dispersion
  real<lower=0> sdste;    // sd of site effects
  real<lower=0> sdbeta[nknots_year];    // sd of GAM coefficients among strata 
  real<lower=0> sdBETA;    // sd of GAM coefficients
  real<lower=0> sdyear;    // sd of year effects
  real<lower=0> sdseason[2];    // sd of season effects

  vector[nknots_year] BETA_raw;//_raw; 
  matrix[nstrata,nknots_year] beta_raw;         // GAM strata level parameters

  vector[nknots_season] beta_raw_season_1;         // GAM coefficients
  vector[nknots_season] beta_raw_season_2;         // GAM coefficients

}

transformed parameters { 
  vector[ncounts] E;           // log_scale additive likelihood
  matrix[nstrata,nknots_year] beta;         // spatial effect slopes (0-centered deviation from continental mean slope B)
  matrix[nyears,nstrata] smooth_pred;
  vector[nyears] SMOOTH_pred;  

  vector[nknots_year] BETA;

  vector[ndays] season_smooth1 = season_basis*(sdseason[1]*beta_raw_season_1);
  vector[ndays] season_smooth2 = season_basis*(sdseason[2]*beta_raw_season_2);
  matrix[ndays,2] season_smooth;
  vector[nyears] yeareffect;             // continental year-effects
 
  BETA = sdBETA*BETA_raw;
  
  for(k in 1:nknots_year){
    beta[,k] = (sdbeta[k] * beta_raw[,k]) + BETA[k];
  }
  SMOOTH_pred = year_basis * BETA; 
  
      for(s in 1:nstrata){
     smooth_pred[,s] = year_basis * transpose(beta[s,]);
  }
  yeareffect = sdyear*yeareffect_raw;
  

    season_smooth[,1] = season_smooth1;
    season_smooth[,2] = season_smooth2;
 
  for(i in 1:ncounts){
          real noise = sdnoise*noise_raw[i];             // over-dispersion
          real ste = sdste*ste_raw[site[i]]; // site intercepts

    // E[i] = ALPHA1 + year_pred[year_raw[i],strat[i]] + alpha[site[i]] + year_effect[year_raw[i]] + season_pred[date[i],seas_strat[i]] + noise;
    E[i] = STRATA + smooth_pred[year[i],strat[i]] + ste + yeareffect[year[i]] + season_smooth[date[i],seas_strat[i]] + noise;
  }
  
  }
  
  
model { 
  sdnoise ~ std_normal(); //prior on scale of extra Poisson log-normal variance
  ste_raw ~ student_t(3,0,1); // fixed site-effects
  sum(ste_raw) ~ normal(0,0.001*nsites);//sum to zero constraint on site-effects
  noise_raw ~ student_t(3,0,1);//std_normal(); // extra Poisson log-normal variance
  sum(noise_raw) ~ normal(0,0.001*ncounts);//sum to zero constraint on site-effects
  sdyear ~ normal(0,0.2); //prior on scale of annual fluctuations - 
  // above is informative so that 95% of the prior includes yearly fluctuations fall
  // between 33% decrease and a 50% increase
  sdste ~ normal(0,2); //prior on scale of site level variation
  
  //sdste ~ student_t(3,0,2); //prior on scale of site level variation
  sdBETA ~ student_t(10,0,1); // prior on sd of GAM parameters


 // sdbeta ~ normal(0,1); //prior on sd of gam hyperparameters
  sdbeta ~ gamma(2,4);//boundary avoiding prior 
  //sdbeta ~ student_t(3,0,1);// prior on spatial variation of spline parameters 
  
  sdseason ~ std_normal();//variance of GAM parameters
  beta_raw_season_1 ~ std_normal();//GAM parameters
  beta_raw_season_2 ~ std_normal();//GAM parameters
  STRATA ~ std_normal();// overall species intercept 
 
  BETA_raw ~ std_normal();// prior on GAM hyperparameters
  yeareffect_raw ~ std_normal(); //prior on ▲annual fluctuations
  sum(yeareffect_raw) ~ normal(0,0.001*nyears);//sum to zero constraint on year-effects

  
for(k in 1:nknots_year){
    beta_raw[,k] ~ icar_normal(nstrata, node1, node2);
}
 
   count ~ poisson_log(E); //vectorized count likelihood
  

}

generated quantities {
  real<lower=0> N[nyears];
  real<lower=0> NSmooth[nyears];
  real<lower=0> n[nstrata,nyears];
  real<lower=0> nsmooth[nstrata,nyears];
  real<lower=0> retrans_noise;
  real<lower=0> retrans_year;
  real<lower=0> retrans_ste;
  
 real seas_max1 = max(season_smooth1)/2;
 real seas_max2 = max(season_smooth2)/2;
 vector[2] seas_max;
 
 seas_max[1] = seas_max1;
 seas_max[2] = seas_max2;
 
 retrans_noise = 0.5*(sdnoise^2);
 retrans_year = 0.5*(sdyear^2);
 retrans_ste = 0.5*(sdste^2);
 
for(y in 1:nyears){

      for(s in 1:nstrata){

   real n_t[nsites_strata[s]];
   real nsmooth_t[nsites_strata[s]];
           
       for(t in 1:nsites_strata[s]){
          real ste = sdste*ste_raw[ste_mat[s,t]]; // site intercepts
      n_t[t] = exp(STRATA + smooth_pred[y,s] + ste + yeareffect[y] + seas_max[seasons[s,2]] + retrans_noise);
      nsmooth_t[t] = exp(STRATA + smooth_pred[y,s] + ste + seas_max[seasons[s,2]] + retrans_noise);

        //   atmp[j] = exp(ALPHA1 + year_pred[y,s] + year_effect[y] + seas_max[seasons[s,2]] + 0.5*(sdnoise^2) + alpha[sites[j,s]]);
        //   atmp_smo[j] = exp(ALPHA1 + year_pred[y,s] + seas_max[seasons[s,2]] + 0.5*(sdyear^2)  + 0.5*(sdnoise^2) + alpha[sites[j,s]]);
         }
        n[s,y] = mean(n_t); //mean of exponentiated predictions across sites in a stratum
        nsmooth[s,y] = mean(nsmooth_t); //mean of exponentiated predictions across sites in a stratum
    }
  }
  
    for(y in 1:nyears){

      N[y] = exp(STRATA + SMOOTH_pred[y] + yeareffect[y] + seas_max[1] + retrans_ste + retrans_noise );
      NSmooth[y] = exp(STRATA + SMOOTH_pred[y] + seas_max[1] + retrans_year + retrans_ste + retrans_noise );
      
    }
    
}

