// Stan object to produce mixture parameters for a new data set Mixture Model (Soft KMeans)

data {
  int<lower=0> N_tilde; // number of data points for predictions
  int<lower=1> M; // number of posterior draws
  int<lower=1> D; // number of dimensions
  int<lower=1> K; // number of clusters
  vector[D] x_tilde[N_tilde]; // inputs for data points for predictions
  
  ordered[K] sigma[M];
  matrix[K, D] beta[M];
}
model{
  
}
generated quantities {
  vector[K] mixture_vec[M, N_tilde];
  vector[K] mixture_vecmeans[N_tilde]; 
  
  for(m in 1:M) {
    for(n in 1:N_tilde) {
      vector[K] mix;
      mix = softmax(beta[m]*x_tilde[n]);
      mixture_vec[m, n] = mix;
    }
  }
  for(n in 1:N_tilde){
    for(k in 1:K) {
      mixture_vecmeans[n, k] = mean(mixture_vec[, n, k]);
    }
  }
}  
  

