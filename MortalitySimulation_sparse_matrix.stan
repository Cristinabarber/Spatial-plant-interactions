data {
  int<lower=0> N; //number of sampled sites or trees
  int<lower=0> N_p; //number of plots
  int<lower=0,upper=1> surv[N]; //survival of each trees
  vector[N] dbh;    //sontinuous Variable we think that we believe that affects mortality
  vector [N] siz[N];        // matrix of size
  vector [N] dist[N];        // matrix of distances
  int  pl[N];          // random effect of plots
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
  real<lower=0,upper=1>  psi[N];    
  vector[N] cri;
  vector[N] cri_scale;
  vector [N] d_vector[N];
  real bp0 [N_p]; 

    for (i in 1:N_p){
    bp0[i]=bp0_m+bp0_sigma*zp[i];
  }
  for (o in 1:N){
    for (l in 1:N){
      d_vector[o,l]=(-1/rho^2)*dist[o,l]^2;
    }
  }
    for(k in 1:N) {
    cri[k] = sum(siz[k] .* exp(d_vector[k]));
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
