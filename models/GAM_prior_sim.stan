// simple GAM prior simulation

data {
 
 int<lower=1> nyears;
  
  real<lower=0>  prior_scale;
  int<lower=0,upper=2> pnorm; // indicator for the prior distribution 2 = t, 1 = normal, 0 = gamma


  // data for spline s(year)
  int<lower=1> nknots_year;  // number of knots in the basis function for year
  matrix[nyears, nknots_year] year_basis; // basis function matrix
 

}



parameters {
 //  vector[ncounts] noise_raw;             // over-dispersion
 // 

  real<lower=0> sdbeta;    // sd of GAM coefficients among strata 

  vector[nknots_year] BETA_raw;         // GAM strata level parameters

}

 
transformed parameters { 
  vector[nknots_year] BETA;         // spatial effect slopes (0-centered deviation from continental mean slope B)
  vector[nyears] smooth_pred;


    BETA = sdbeta * BETA_raw;

 
     smooth_pred = year_basis * BETA;



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

   BETA_raw ~ normal(0,1);



}

 generated quantities {

   vector[nyears] nsmooth;


 
      nsmooth = exp(smooth_pred);
        



  }



 

