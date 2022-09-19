data{ 
 int N;      // number of individuals 
 int obs;    // number of non-zero values in the sparse matrices
 vector [N] size_t0;   // size of focal plants 
 vector [N] growth;   // growth of focal plants, response 
 vector [obs] size_vector;  //vector of non-zero size observations 
 vector [obs] dist_vector;  //vector of non-zero distance observations  
 int pos[N];     // order of the first non-zero values
 int n_nb[N];     //number of non-zero values per row
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
 vector[obs] size_vec;
 vector[obs] dist_vec;
 for (i in 1:obs){
  dist_vec[i]=dist_observations[i]^2; 
  size_vec[i]=size_observations[i]^a1;
 }
 for(n in 1:N)
  kernel[n]=sum(segment(size_vec, pos[n], n_nb[n]))/
       exp(segment(dist_vec, pos[n], n_nb[n])*a2));
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
