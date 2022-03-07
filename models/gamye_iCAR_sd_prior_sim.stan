// This is a Stan implementation of the gamye model
// with iCAR component for the stratum-level intercepts and smooth parameters

// Consider moving annual index calculations outside of Stan to 
// facilitate the ragged array issues

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
 
  int<lower=1> nstrata;
  int<lower=1> nyears;
  
  real<lower=0>  prior_scale;
  int<lower=0,upper=2> pnorm; // indicator for the prior distribution 2 = t, 1 = normal, 0 = gamma

 // spatial neighbourhood information
  int<lower=1> N_edges;
  int<lower=1, upper=nstrata> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=nstrata> node2[N_edges];  // and node1[i] < node2[i]



  // data for spline s(year)
  int<lower=1> nknots_year;  // number of knots in the basis function for year
  matrix[nyears, nknots_year] year_basis; // basis function matrix
 

}



parameters {
 //  vector[ncounts] noise_raw;             // over-dispersion
 // 

  real<lower=0> sdbeta[nknots_year];    // sd of GAM coefficients among strata 

  matrix[nstrata,nknots_year] beta_raw;         // GAM strata level parameters

}

 
transformed parameters { 
  matrix[nstrata,nknots_year] beta;         // spatial effect slopes (0-centered deviation from continental mean slope B)
  matrix[nyears,nstrata] smooth_pred;


  for(k in 1:nknots_year){
    beta[,k] = (sdbeta[k] * beta_raw[,k]);
  }
 
      for(s in 1:nstrata){
     smooth_pred[,s] = year_basis * transpose(beta[s,]);
}


// intercepts and slopes

  



  }
  
  
  
model {


if(pnorm == 1){
 sdbeta ~ normal(0,prior_scale); //prior on sd of GAM parameter variation
}
if(pnorm == 0){
 
  sdbeta ~ gamma(2,prior_scale); //prior on sd of GAM parameter variation
 
}
if(pnorm == 2){
 
  sdbeta ~ student_t(3,0,prior_scale); //prior on sd of GAM parameter variation
 
}

for(k in 1:nknots_year){
    beta_raw[,k] ~ icar_normal(nstrata, node1, node2);
}



}

 generated quantities {

   real<lower=0> nsmooth[nstrata,nyears];

for(y in 1:nyears){

      for(s in 1:nstrata){

 
      nsmooth[s,y] = exp(smooth_pred[y,s]);
        }



    }
  }



 

