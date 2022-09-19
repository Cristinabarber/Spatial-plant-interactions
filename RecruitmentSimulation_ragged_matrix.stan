data {
  int<lower=0> N; //Amount of sampled sites or trees
  int<lower=0> obs; //Amount of non zero observations
  int<lower=0> N_p;
  int<lower=0> recruit[N]; //Survival of each trees
  int n_nb[N];
  int pos[N];
  vector [obs] size_vector;        // long vector of pairwise
  vector [obs] dist_vector;        // long vector of pairwise distances
  int  pl[N];
}
parameters{
  real b2; 
  real bp0_m; 
  real<lower=0> bp0_sigma;
  vector[N_p] zp; //Scale parameter
  real<lower=0> rho;
  
}

transformed parameters{
  real<lower=0>  mu[N];   
  vector[N] cri;
  vector[N] cri_scale;
  real bp0 [N_p];
  vector [obs] d_vector;

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
    mu[n]=exp(bp0[pl[n]] + b2*cri_scale[n]);
    }
}
model{
  b2~normal(0,1);
  rho~gamma(3,2);
  bp0_sigma~normal(0,1); 
  bp0_m~normal(0,1);
  zp~normal(0,1);
  

  recruit~poisson(mu);
}
