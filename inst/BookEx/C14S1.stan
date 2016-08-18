functions{
// downside function
real f1(real x, real Lambda, real eta){
real valf1;
valf1 <- (Lambda - 1) * exp(-x * eta^2) / sqrt(2 * pi());
return(valf1);
}
// upside function
real f2(real x, real Lambda, real eta){
real valf2;
valf2 <- (1 + Lambda) / 2 - (Lambda - 1) / 2 * erfc(sqrt(x) * eta);
return(valf2);
}
// utility function (log probabilistic specification)
real pu_log(vector w, vector mu, matrix Sigma, real Lambda, real L, int iperiod, real nu){
real ll;
real M;
real S;
real eta;
real dt;
vector[iperiod] val;
M <- w' * mu - L;
S <- sqrt(w' * Sigma * w);
eta <- M / (2 * S);
dt <- 1 / (1.0 * iperiod);
for(i in 1:iperiod){
val[i] <- i * dt * M * f2(i * dt, Lambda, eta) - sqrt(i * dt) * S * f1(i * dt, Lambda, eta);
}
return(nu * mean(val));
}
}

data{
int<lower=0> N; // number of assets
int<lower=0> iperiod; // investment period
real<lower=0> nu; // convergence rate
vector[N] mu; // vector of expected returns
matrix[N, N] Sigma; // variance-covariance matrix
real<lower=0> Lambda; // risk aversion
real<lower=0> L; // target portfolio return
}

parameters {
simplex[N] w; // weight vector
}

model{
increment_log_prob(pu_log(w, mu, Sigma, Lambda, L, iperiod, nu));
}