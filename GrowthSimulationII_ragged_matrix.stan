data {
  int<lower=0> N;   //amount of sampled trees
  int<lower=0> N_p;                 // number of plots
  int<lower=0> obs; //total amount of non zero observations in local neighborhoods
  vector[N] gr;    //survival of each trees
  vector[N] dbh;  //initial DBH
  int n_nb[N];    //number of non-zero values per row
  int pos[N];     //vector giving the position of non-zero values
  vector [obs] dist_vector;        // vector of distance to neighbors
  vector [obs] size_vector;        // vector of neighbor sizes
  int pl[N];                      // random effect of plots
}
transformed data{
  // square distances
  vector [obs] d_vector;
  for(i in 1:obs){
    d_vector[i] = dist_vector[i]^2;
  }
}
parameters{
  real b0; 
  real b1; 
  real b2; 
  real kb0_m; 
  real<lower=0> k_sigma; 
  real<lower=0> sigma; 
  real rho;
  real zp[N_p];
}
transformed parameters{
  vector[N] mu;
  vector[N] cri;
  vector[N] cri_scale;
  real bp0 [N_p];

  for (i in 1:N_p){
    bp0[i]=kb0_m+k_sigma*zp[i];
  }

  for(n in 1:N) {
    cri[n] = sum( segment(size_vector, pos[n], n_nb[n]) .* exp(-1/rho^2 * segment(d_vector, pos[n], n_nb[n])));
  }
  cri_scale = (cri - mean(cri)) / (2*sd(cri));
  for (n in 1:N){
    mu[n] = bp0[pl[n]] + b1 * dbh[n] + b2 * cri_scale[n];
  }
}
model{
  b0 ~ normal(0, 1);
  b1 ~ normal(0, 1);
  b2 ~ normal(0, 1);
  sigma ~ exponential(2);
  k_sigma ~ exponential(2);
  rho ~ normal(2, .5);
  kb0_m~normal(0,1);
  zp~normal(0,1);
  gr ~ normal(mu, sigma);
}
