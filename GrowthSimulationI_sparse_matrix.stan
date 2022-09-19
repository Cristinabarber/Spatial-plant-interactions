data{ 
 int N;      // number of individuals 
 vector [N] size_t0;   // size of focal plants 
 vector [N] growth;   // growth of focal plants, response 
 vector [N] sizemat[N];  // full size matrix 
 vector [N] distmat[N];  // full distance matrix 
}
parameters{
 real alpha;      
 real beta;      
 real sigma;      
 real<lower=0> a1;     
 real a3;       
 real<lower=0> a2;     
 }
transformed parameters{ 
 vector[N] kernel; 
 vector[N] mu;
 vector[N] smat[N];
 vector[N] dmat[N];
 for(i in 1:N){
 for(j in 1:N){
  smat[i,j]=sizemat[i,j]^a1; 
  dmat[i,j]=distmat[i,j]^2*a2;
 }}
 for(n in 1:N)
  kernel[n]=sum(smat[n]./exp(dmat[n]));
 for(n in 1:N)
  mu[n]=alpha+size_t0[n]*beta+a3*kernel[n];
}
model{
 alpha~normal(0,5); 
 beta~normal(0,5);
 a1~normal(0,5);
 a2~normal(0,5);
 a3~normal(0,5);
 sigma~exponential(1);
 
 growth ~ normal(mu,sigma);
}
