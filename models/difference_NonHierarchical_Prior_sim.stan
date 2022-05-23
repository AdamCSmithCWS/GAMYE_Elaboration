// This is a Stan implementation of a prior simulation of
// the first difference time-series


data {
  int<lower=1> nstrata;
  int<lower=1> nyears;
  int<lower=1> nyears_m1;

  int<lower=3> df; // df of t-prior
  real<lower=0>  prior_scale_b; //scale of the prior distribution
  int<lower=0,upper=1> pnorm; // indicator for the prior distribution 0 = t, 1 = normal
 
 

  // data to center the abundance estimate
  int midyear; //middle year of the time-series scaled to ~(nyears/2)
  int nIy1; //indexing vector dimension - number of years before midyear
  int nIy2; //indexing vector dimension - number of years after midyear
  int Iy1[nIy1];//indexing vector - indices of years before midyear
  int Iy2[nIy2];//indexing vector - indices of years after midyear
  
  // just a vector of zeros to fill fixed beta values for midyear
  vector[nstrata] zero_betas;
}

parameters {
  
  real<lower=0> sdbeta[nstrata];    // sd of annual changes among strata 

  matrix[nstrata,nyears_m1] beta_raw;         // strata level parameters

}

transformed parameters { 
  matrix[nstrata,nyears] beta;         // spatial effect slopes (0-centered deviation from continental mean slope B)

  matrix[nstrata,nyears] yeareffect;


  beta[,midyear] = zero_betas;
  yeareffect[,midyear] = zero_betas;

for(s in 1:nstrata){
  for(t in Iy1){
    beta[s,t] = (sdbeta[s] * beta_raw[s,t]);
    yeareffect[s,t] = yeareffect[s,t+1] + beta[s,t];
  }
 
   for(t in Iy2){
    beta[s,t] = (sdbeta[s] * beta_raw[s,t-1]);//t-1 indicators to match dimensionality
    yeareffect[s,t] = yeareffect[s,t-1] + beta[s,t];
  }
 
}
  
  }
  
  
  
model {
 
//Conditional statements to select the prior distribution
if(pnorm == 1){
 sdbeta ~ normal(0,prior_scale_b); //prior on sd of GAM parameter variation
}
if(pnorm == 0){
  sdbeta ~ student_t(df,0,prior_scale_b); //prior on sd of GAM parameter variation
}



for(s in 1:(nstrata)){
    beta_raw[s,] ~ std_normal();
  //sum(beta_raw[s,]) ~ normal(0,0.001*nyears);

}

}

 generated quantities {

  //estimated trajectory on a count-scale
      real<lower=0> n[nstrata,nyears];

for(y in 1:nyears){

      for(s in 1:nstrata){
      n[s,y] = exp(yeareffect[s,y]);
        }

    }
 
 }

