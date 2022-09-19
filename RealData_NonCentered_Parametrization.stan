data{
 int N;      //number of plots
 int K;      //number of non-zero parent trees
 int M;      //number of random levels
 vector [K] sizeN;    //matrix of neighbor size
 vector [K] distN;    //matrix of neighbor distances
 int x[N];      //number of seedlings
 int seeds[N];     //number of seeds
 int am[N];     //vector giving the number of non-zero values
 int pos [N];    //vector giving the position of non-zero values
 int Cseedlings [N];   //number of conspecific seedlings
 int plots[N];     //random effect of plots
}
parameters {
 real a;
 real b;
 real<lower=0> ger;
 real mu;
 real mu_omega;
 real slope_omega[M];
 real scale_omega;

 
}
transformed parameters{
 real<lower=0, upper=1> s[N];
 real g[N];
 real omega [M];
 
 for (n in 1:M){
 omega[n]= mu_omega +slope_omega[n]*scale_omega;
 }
 for (n in 1:N){
 if (am[n]==0){g[n]=mu+ b*Cseedlings[n]+omega[plots[n]];}
 else{
  g[n] = mu + b*Cseedlings[n] + a* sum(segment(sizeN,pos[n],am[n]) ./
    exp(ger*log(segment(distN,pos[n],am[n])))) + omega[plots[n]];
  }
 }
}
model{
 a~normal(0,1);
 b~normal(0,1);
 ger~normal(0,1);
 mu~normal(0,1);
 slope_omega~normal(0,1);
 mu_omega~normal(0,1);
 scale_omega~normal(0,1);
 
 for (n in 1:N){
 x~binomial_logit(seeds,g);
 }
}
