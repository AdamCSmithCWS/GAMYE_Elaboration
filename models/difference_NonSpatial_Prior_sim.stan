// This is a Stan implementation of a prior simulation of
// the first difference time-series


data {
  int<lower=1> nstrata;
  int<lower=1> nyears;
  int<lower=1> nyears_m1;

  int<lower=3> df; // df of t-prior
  real<lower=0>  prior_scale_B; //scale of the prior distribution
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
  
  real<lower=0> sdbeta;    // sd of annual changes among strata 
  real<lower=0> sdBETA;    // sd of overall annual changes

  vector[nyears_m1] BETA_raw;//_raw; 
  matrix[nstrata,nyears_m1] beta_raw;         // strata level parameters

}

transformed parameters { 
  matrix[nstrata,nyears] beta;         // spatial effect slopes (0-centered deviation from continental mean slope B)

  matrix[nstrata,nyears] yeareffect;
  vector[nyears_m1] BETA;
  vector[nyears] YearEffect;

  BETA = sdBETA * BETA_raw;

  beta[,midyear] = zero_betas;
  yeareffect[,midyear] = zero_betas;
  YearEffect[midyear] = 0;

  for(t in Iy1){
    beta[,t] = (sdbeta * beta_raw[,t]) + BETA[t];
    yeareffect[,t] = yeareffect[,t+1] + beta[,t];
    YearEffect[t] = YearEffect[t+1] + BETA[t]; 
  }
 
   for(t in Iy2){
    beta[,t] = (sdbeta * beta_raw[,t-1]) + BETA[t-1];//t-1 indicators to match dimensionality
    yeareffect[,t] = yeareffect[,t-1] + beta[,t];
    YearEffect[t] = YearEffect[t-1] + BETA[t-1]; 
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

   BETA_raw ~ std_normal(); //non-centered parameterisation


for(t in 1:(nyears_m1)){
    beta_raw[,t] ~ std_normal();
  sum(beta_raw[,t]) ~ normal(0,0.001*nstrata);

}

}

 generated quantities {

  //estimated trajectory on a count-scale
   vector[nyears] N = exp(YearEffect);
       real<lower=0> n[nstrata,nyears];

for(y in 1:nyears){

      for(s in 1:nstrata){
      n[s,y] = exp(yeareffect[s,y]);
        }

    }
 
 }

