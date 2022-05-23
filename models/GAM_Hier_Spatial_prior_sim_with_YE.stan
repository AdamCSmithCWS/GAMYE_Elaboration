// simple GAM prior simulation


functions {
  real icar_normal_lpdf(vector bb, int ns, int[] n1, int[] n2) {
    return -0.5 * dot_self(bb[n1] - bb[n2])
      + normal_lpdf(sum(bb) | 0, 0.001 * ns); //soft sum to zero constraint on bb
 }
}

data {
  int<lower=1> nstrata;
  int<lower=1> nyears;
  int<lower=3> df; // df of t-prior
  real<lower=0>  prior_scale_B; //scale of the prior distribution
  real<lower=0>  prior_scale_b; //scale of the prior distribution
  real<lower=0>  prior_scale_y; //scale of the prior distribution for year-effects
  int<lower=0,upper=1> pnorm; // indicator for the prior distribution 0 = t, 1 = normal
  
   // spatial neighbourhood information
  int<lower=1> N_edges;
  int<lower=1, upper=nstrata> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=nstrata> node2[N_edges];  // and node1[i] < node2[i]

  // data for spline s(year)
  int<lower=1> nknots_year;  // number of knots in the basis function for year
  matrix[nyears, nknots_year] year_basis; // basis function matrix
}

parameters {

  real<lower=0> sdyear[nstrata];    // sd of GAM coefficients among strata 
  matrix[nstrata,nyears] yeareffect_raw;         // GAM strata level parameters
  
  real<lower=0> sdbeta[nknots_year];    // sd of GAM coefficients among strata 
  matrix[nstrata,nknots_year] beta_raw;         // GAM strata level parameters

  real<lower=0> sdBETA;    // sd of spline coefficients  
  vector[nknots_year] BETA_raw;         // unscaled spline coefficients
}
 
transformed parameters { 
  vector[nknots_year] BETA;         // spatial effect slopes (0-centered deviation from continental mean slope B)
  vector[nyears] SMOOTH_pred;

  matrix[nstrata,nknots_year] beta;         // spatial effect slopes (0-centered deviation from continental mean slope B)
  matrix[nyears,nstrata] smooth_pred;

  matrix[nstrata,nyears] yeareffect;         // GAM strata level parameters
  
    BETA = sdBETA * BETA_raw; //scaling the spline parameters
     SMOOTH_pred = year_basis * BETA; //log-scale smooth trajectory
  
  for(k in 1:nknots_year){
    beta[,k] = (sdbeta[k] * beta_raw[,k]) + BETA[k];
  }
    
   for(s in 1:nstrata){
   
    yeareffect[s,] = (sdyear[s] * yeareffect_raw[s,]);
    
  } 
  
      for(s in 1:nstrata){
     smooth_pred[,s] = year_basis * transpose(beta[s,]);
}


  }
  
model {

//Conditional statements to select the prior distribution
if(pnorm == 1){
 sdBETA ~ normal(0,prior_scale_B); //prior on sd of GAM parameter variation
 sdbeta ~ normal(0,prior_scale_b); //prior on sd of GAM parameter variation
}
if(pnorm == 0){
  sdBETA ~ student_t(df,0,prior_scale_B); //prior on sd of GAM parameter variation
  sdbeta ~ student_t(df,0,prior_scale_b); //prior on sd of GAM parameter variation
}

   BETA_raw ~ normal(0,1); //non-centered parameterisation

 sdyear ~ gamma(2,prior_scale_y); //prior on sd of yeareffects

for(k in 1:nknots_year){
    beta_raw[,k] ~ icar_normal(nstrata, node1, node2);
}

for(s in 1:nstrata){
   yeareffect_raw[s,] ~ normal(0,1);
    sum(yeareffect_raw[s,]) ~ normal(0,0.001*nyears);
}


}

 generated quantities {
  //estimated smooth on a count-scale
   vector[nyears] NSMOOTH = exp(SMOOTH_pred);
    
    
       real<lower=0> nsmooth[nstrata,nyears];
        real<lower=0> n[nstrata,nyears];

for(y in 1:nyears){

      for(s in 1:nstrata){

 
      nsmooth[s,y] = exp(smooth_pred[y,s]);
      n[s,y] = exp(smooth_pred[y,s] + yeareffect[s,y]);
        }



    }
    
  }



 

