data {
 int N;      //number of available places for recruitment
 int K;      //number of parent trees 
 vector [K] dist[N];   //array containing N vectors with K distances 
 int x[N];      //count of seedlings
 int CP[N];     //count of cabbage palms
 vector [N] ones;    //vector of ones for the division
}
parameters {
 real<lower=0> a;
 real<lower=0> b;
 real<lower=0> c;
 real<lower=0> phy;
}
transformed parameters{
 real mu[N];
 real aa[N];
 
 for(i in 1:N) {
 aa[i] = sum(ones[i]./(c+dist[i]));
 }
 for (n in 1:N){
  mu[n]=(a+b*aa[n])*CP[n];
 }
}
model{
 a~normal(0,100);
 b~normal(0,100);
 c~normal(0,100);
 phy~exponential(0.5);
 
 x~neg_binomial_2(mu,phy);
 }

