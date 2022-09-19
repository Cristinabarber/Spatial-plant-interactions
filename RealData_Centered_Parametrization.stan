data {
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
 real e [M];
 real sigma_plot;
 
}
transformed parameters{
 real<lower=0, upper=1> s[N];
 real g[N];

 for (n in 1:N){
 if (am[n]==0){g[n]=mu+ b*Cseedlings[n]+e[plots[n]];}
 else{
  g[n] = mu + b*Cseedlings[n] + a* sum(segment(sizeN,pos[n],am[n]) ./
    exp(ger*log(segment(distN,pos[n],am[n])))) + e[plots[n]];
  }
 }
 
}
model{
 a~normal(0,1);
 b~normal(0,1);
 ger~normal(0,1);
 mu~nomal(0,1);
 e~normal(0,sigma_plot);
 sigma_plot~normal(0,1);

 x~binomial_logit(seeds,g);
 }

