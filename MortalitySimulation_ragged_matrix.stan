data {
  int<lower=0> N; //Amount of sampled sites or trees
  int<lower=0> obs; //amount of non zero observations
  int<lower=0> N_p; //number of plots
  int<lower=0,upper=1> surv[N]; //survival of each trees
  vector[N] dbh;               //initial DBH
  int n_nb[N];    //number of non-zero values per row
  int pos[N];     //vector giving the position of non-zero values
  vector [obs] size_vector;        // vector of size
  vector [obs] dist_vector;        // vector of distances
  int  pl[N]; // random effect of plots
}
parameters{
  real b0; 
  real b1;
  real b2; 
  real bp0_m; 
  real<lower=0> bp0_sigma; 
  vector[N_p] zp; 
  real<lower=0> rho;
}
transformed parameters{
  real<lower=0,upper=1>  psi[N];    //probability survival
  vector[N] cri;
  vector[N] cri_scale;
  vector [obs] d_vector;
  real bp0 [N_p]; //random effect coefficients

    for (i in 1:N_p){
    bp0[i]=bp0_m+bp0_sigma*zp[i];
  }
  for (l in 1:obs){
  d_vector[l]=dist_vector[l]^2;
  }
    for(k in 1:N) {
    cri[k] = sum(segment(size_vector, pos[k], n_nb[k]) .* exp(segment(d_vector, pos[k], n_nb[k])));
  }
  cri_scale = (cri - mean(cri)) / (2*sd(cri));

  for (n in 1:N){
    psi[n]=inv_logit(b0+bp0[pl[n]] + b1 * dbh[n] + b2 * cri_scale[n]);
    }
}
model{
  b0~normal(0,1);
  b1~normal(0,1);
  b2~normal(0,1);
  rho~normal(0,1);
  bp0_m~normal(0,1);
  bp0_sigma~normal(0,1); //Unit verctor, always follows the normal(0,1) distribution
  zp~normal(0,1);

surv~bernoulli_lpmf(psi);
}
