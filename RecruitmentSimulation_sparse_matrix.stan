data {
  int<lower=0> N; //number of sampled plots
  int<lower=0> M; //number of parent trees
  int<lower=0> N_p; //number of plots
  int<lower=0> recruit[N]; //number of recruitments on each plot
  vector [N] siz[M];        // vector of size
  vector [N] dist[M];        // vector of distances
  int  pl[N];  // random effect of plots
}
parameters{
  real b2; 
  real b0;
  real bp0_m; 
  real<lower=0> bp0_sigma; 
  real<lower=0> rho;
  vector[N_p] zp;
}
transformed parameters{
  real<lower=0>  mu[N];   
  vector[N] cri;
  vector[N] cri_scale;
  vector [N] d_vector[N];
  real bp0 [N_p]; 

    for (i in 1:N_p){
    bp0[i]=bp0_m+bp0_sigma*zp[i];
  }
  for (o in 1:N){
    for (l in 1:N){
    d_vector[o,l] = (-1/rho^2) *dist[o,l]^2;
  }
 }
  for(k in 1:N) {
    cri[k] = sum(siz[k] .* exp(d_vector[k]));
  }
  cri_scale = (cri - mean(cri)) / (2*sd(cri));
  for (n in 1:N){
    mu[n]=exp(bp0[pl[n]] + b2*cri_scale[n]);
    }
}
model{
  b2~normal(0,1);
  rho~gamma(3,2);
  b0~normal(0,1);
  bp0_sigma~exponential(1); 
  bp0~normal(b0,bp0_sigma);

  recruit~poisson(mu);
}
