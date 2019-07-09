// Stan object to produce Leave-One-Out predictions with nugget option
data {
  int<lower=1> N1; // number of ensemble members in design set.
  int<lower=1> p; // number of inputs.
  int<lower=1> M; // number of posterior draws.
  int<lower=1> Np; // number of regression functions.
  int<lower=1> L; // number of centres of mass.
  
  row_vector[p] X1[N1]; // design matrix.
  vector[N1] y1; // vector of simulator evaluations at design matrix.
  matrix[N1, Np] H1; //regression matrix for design matrix.
  vector[N1] A[L]; // array of weights vectors.

  vector[Np] beta[M];
  row_vector[p] delta[M, L];
  real sigma[M, L];
  real nugget[M, L];
}
transformed data {
  real length_scale;
  int LOO[N1, N1-1];
  
  length_scale = pow(sqrt(2.0), -1 );
  for(i in 1:N1) {
    for(j in 1:(i-1)) LOO[i, j] = j;
    for(j in (i+1):N1) LOO[i, j-1] = j;
  }
}
model {
}
generated quantities {
  vector[N1] tmeans;
  vector[N1] tsds;
  matrix[M, N1] predict_y;
  for(m in 1:M) {
    matrix[N1, N1] K;
    vector[N1] Mu;
    vector[N1] e;
    row_vector[p] XS1[N1];
    // compute mean for design matrix.
    Mu = H1*beta[m];
    // compute error vector
    e = y1 - Mu;
    for(n in 1:N1) XS1[n] = X1[n] ./ delta[m, 1];
    K = quad_form_diag(cov_exp_quad(XS1, sigma[m, 1], length_scale), A[1]);
    for(l in 2:L) {
      row_vector[p] XS[N1];
      for(n in 1:N1) XS[n] = X1[n] ./ delta[m, l];
      K = K + quad_form_diag(cov_exp_quad(XS, sigma[m, l], length_scale), A[l]);
    }
    for(i in 1:N1) K[i, i] = K[i, i] + nugget[m, sort_indices_desc(A[, i])[1]];
    for(n in 1:N1) {
      real mu;
      real sigma_squared;
      matrix[N1-1, N1-1] Sigma;// new covariance matrix. 
      row_vector[N1-1] kt;
      row_vector[N1-1] k;
      
      Sigma = K[LOO[n], LOO[n]];
      kt = K[n, LOO[n]];
      k = mdivide_right_spd(kt, Sigma);
 
      mu = Mu[n] + k*e[LOO[n]];
      sigma_squared = K[n, n] - k*kt';
      predict_y[m, n] = normal_rng(mu, sqrt(sigma_squared));
   }
  }
    for(k in 1:N1){
      tmeans[k] = mean(predict_y[,k]);
      tsds[k] = sd(predict_y[,k]);
  }
}
    
