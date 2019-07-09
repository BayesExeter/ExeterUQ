data {
  int<lower=1> N; // number of data points.
  vector[N] If; // vector of field implausibilities.
  vector[N] Ic; // vector of coefficient implausibilities.
  real<lower=0> beta0prior; // prior variance for intercept
  real<lower=0> alpha0prior; // prior variance for intercept (log model)
  real<lower=0> beta1shape; // prior shape parameter for If slope
  real<lower=0> beta1scale; // prior scale parameter for If slope
  real<lower=0> alpha1shape; // prior shape parameter for Ic slope
  real<lower=0> alpha1scale; // prior scale parameter for Ic slope
  real<lower=0> T_f; // bound on the field
}
parameters{
  real<lower=0> sigma_sq;
  real beta0;
  real beta1;
  real alpha0;
  real alpha1;
}
transformed parameters{
  vector[N] Mu;
  vector[N] Sigma;
  Mu = beta0 + beta1 * If;
  for(i in 1:N) Sigma[i] = sigma_sq / exp(alpha0 + alpha1 * log(If[i]));
}
model{
  //sigma_sq ~ (1 / sigma_sq);
  sigma_sq ~ gamma(1, 1);
  beta0 ~ normal(0, beta0prior);
  alpha0 ~ normal(0, alpha0prior);
  beta1 ~ gamma(beta1shape, beta1scale);
  -alpha1 ~ gamma(alpha1shape, alpha1scale);
  Ic ~ normal(Mu, Sigma);
}
generated quantities {
  real T_c;
  T_c = beta0 + beta1 * T_f;
}