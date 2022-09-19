data {
 int N;      //number of available places for recruitment
 int K;      //number of non-zero parent trees 
 vector [K] distrag;   //vector containing all the non-zero distances
 int x[N];      //count of seedlings
 int CP[N];     //count of cabbage palms
 int n_nb[N];     //vector giving the amount of non-zero values
 int pos [N];    //vector giving the position of non-zero values
 vector [N] one;    //vector of ones for the division
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
 aa[i] = sum(one[i] ./(c +segment(distrag,pos[i],n_nb[i])));
 }
 
 for (n in 1:N){
 if (n_nb[n]==0){mu[n]=a;}
 else{
  mu[n]=(a+b*aa[n])*CP[n];
 }
 }
}
model{
 a~normal(0,100);
 b~normal(0,100);
 c~normal(0,100);
 phy~exponential(0.5);
 
 x~neg_binomial_2(mu,phy);
}