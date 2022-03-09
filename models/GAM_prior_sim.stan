// simple GAM prior simulation

data {
  int<lower=1> nyears;
  real<lower=0>  prior_scale; //scale of the prior distribution
  int<lower=0,upper=1> pnorm; // indicator for the prior distribution 0 = t, 1 = normal
  // data for spline s(year)
  int<lower=1> nknots_year;  // number of knots in the basis function for year
  matrix[nyears, nknots_year] year_basis; // basis function matrix
}

parameters {
  real<lower=0> sdbeta;    // sd of spline coefficients  
  vector[nknots_year] BETA_raw;         // unscaled spline coefficients
}
 
transformed parameters { 
  vector[nknots_year] BETA;         // spatial effect slopes (0-centered deviation from continental mean slope B)
  vector[nyears] smooth_pred;

    BETA = sdbeta * BETA_raw; //scaling the spline parameters
     smooth_pred = year_basis * BETA; //log-scale smooth trajectory
  }
  
model {

//Conditional statements to select the prior distribution
if(pnorm == 1){
 sdbeta ~ normal(0,prior_scale); //prior on sd of GAM parameter variation
}
if(pnorm == 0){
  sdbeta ~ student_t(3,0,prior_scale); //prior on sd of GAM parameter variation
}

   BETA_raw ~ normal(0,1); //non-centered parameterisation

}

 generated quantities {
  //estimated smooth on a count-scale
   vector[nyears] nsmooth = exp(smooth_pred);
    
  }



 

