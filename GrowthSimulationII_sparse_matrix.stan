data {
  int<lower=0> N;       //number of sampled trees
  int<lower=0> N_p;                 // number of plots
  vector[N] gr;         //survival of each trees
  vector[N] dbh;        //initial tree DBH
  vector [N] dist[N];   //matrix of neighbor distances
  vector [N] siz[N];    //matrix of neighbor size
  int pl[N];            //random effect of plots
}
parameters{
  real b0; 
  real b1; 
  real b2; 
  real kb0_m; 
  real<lower=0> k_sigma; 
  real<lower=0> sigma; 
  real rho;
}
transformed parameters{
  vector[N] mu;
  vector[N] cri;
  vector[N] cri_scale;
  vector [N] d_vector[N];
  real bp0 [N_p]]; 

  for (i in 1:N_p){
    bp0[i]=bp0_m+bp0_sigma*zp[i];
  }
  for (o in 1:N){
    for (l in 1:N){
    d_vector[o,l] = (-1/rho^2) *dist[o,l]^2;
  }
 }
  for(n in 1:N) {
    cri[n] = sum(siz[n] .* exp( d_vector[n]));
  }
  cri_scale = (cri - mean(cri)) / (2*sd(cri));
  for (n in 1:N){
    mu[n] = kb0[pl[n]] + b1 * dbh[n] + b2 * cri_scale[n];
  }
}
model{
  b0 ~ normal(0, 1);
  kb0 ~ normal(b0, k_sigma);
  b1 ~ normal(0, 1);
  b2 ~ normal(0, 1);
  sigma ~ exponential(2);
  k_sigma ~ exponential(2);
  rho ~ normal(2, .5);

  gr ~ normal(mu, sigma);
}
