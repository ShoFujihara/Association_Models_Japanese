# Test corrected Unidiff Bayesian estimation
library(tidyverse)
library(gnm)
library(cmdstanr)
library(bayesplot)

# Data
Freq_yamaguchi <- c(
  1275, 364, 274, 272, 17, 1055, 597, 394, 443, 31,
  1043, 587, 1045, 951, 47, 1159, 791, 1323, 2046, 52,
  666, 496, 1031, 1632, 646, 474, 129, 87, 124, 11,
  300, 218, 171, 220, 8, 438, 254, 669, 703, 16,
  601, 388, 932, 1789, 37, 76, 56, 125, 295, 191,
  127, 101, 24, 30, 12, 86, 207, 64, 61, 13,
  43, 73, 122, 60, 13, 35, 51, 62, 66, 11,
  109, 206, 184, 253, 325)

d_xie <- tibble(
  O = gl(5, k = 5, length = 75),
  D = gl(5, k = 1, length = 75),
  C = gl(3, k = 25, length = 75),
  Freq = Freq_yamaguchi
) |>
  mutate(
    Diag = case_when(O == D ~ O, .default = "0") |> factor(),
    i_O = as.integer(O), i_D = as.integer(D), i_C = as.integer(C),
    cell_OD = (i_O - 1) * 5 + i_D,
    diag_cell_idx = if_else(O == D, (i_C - 1) * 5 + as.integer(O), 0L)
  )

# gnm reference
cat("\n========== gnm (reference) ==========\n")
set.seed(12345)
fit_gnm <- gnm(Freq ~ O*C + D*C + Diag*C + Mult(Exp(C), O*D),
               family = poisson, data = d_xie, tolerance = 1e-12)
coef_gnm <- coef(fit_gnm)
exp_c_coefs <- coef_gnm[grep("^Mult.*Exp.*\\.C", names(coef_gnm))]
gnm_beta <- exp(exp_c_coefs - exp_c_coefs[1])
names(gnm_beta) <- c("US", "Britain", "Japan")
cat("gnm beta:", gnm_beta, "\n")
cat("gnm L²:", deviance(fit_gnm), "df:", df.residual(fit_gnm), "\n")

# Stan model (corrected)
stan_code <- '
data {
  int<lower=1> N; int<lower=1> I; int<lower=1> J; int<lower=1> K;
  int<lower=1> IJ; int<lower=1> N_diag;
  array[N] int<lower=1> row_idx; array[N] int<lower=1> col_idx;
  array[N] int<lower=1> layer_idx; array[N] int<lower=1> cell_idx;
  array[N] int<lower=0, upper=N_diag> diag_cell_idx;
  array[N] int<lower=0> y;
}
parameters {
  real alpha;
  vector[I-1] alpha_row_raw; vector[J-1] alpha_col_raw; vector[K-1] alpha_layer_raw;
  matrix[I-1, K-1] alpha_row_layer_raw; matrix[J-1, K-1] alpha_col_layer_raw;
  vector[N_diag-1] alpha_diag_raw;
  vector[IJ-1] psi_raw;
  vector[K-1] beta_raw;
}
transformed parameters {
  vector[I] alpha_row; vector[J] alpha_col; vector[K] alpha_layer;
  matrix[I, K] alpha_row_layer; matrix[J, K] alpha_col_layer;
  vector[N_diag] alpha_diag; vector[IJ] psi; vector[K] beta; vector[N] log_mu;

  alpha_row[1]=0; alpha_col[1]=0; alpha_layer[1]=0;
  for(i in 2:I) alpha_row[i]=alpha_row_raw[i-1];
  for(j in 2:J) alpha_col[j]=alpha_col_raw[j-1];
  for(k in 2:K) alpha_layer[k]=alpha_layer_raw[k-1];
  for(k in 1:K){alpha_row_layer[1,k]=0; alpha_col_layer[1,k]=0;}
  for(i in 1:I) alpha_row_layer[i,1]=0;
  for(j in 1:J) alpha_col_layer[j,1]=0;
  for(i in 2:I) for(k in 2:K) alpha_row_layer[i,k]=alpha_row_layer_raw[i-1,k-1];
  for(j in 2:J) for(k in 2:K) alpha_col_layer[j,k]=alpha_col_layer_raw[j-1,k-1];
  alpha_diag[1]=0; for(d in 2:N_diag) alpha_diag[d]=alpha_diag_raw[d-1];
  psi[1]=0; for(ij in 2:IJ) psi[ij]=psi_raw[ij-1];
  beta[1]=1; for(k in 2:K) beta[k]=exp(beta_raw[k-1]);

  for(n in 1:N){
    log_mu[n]=alpha+alpha_row[row_idx[n]]+alpha_col[col_idx[n]]+alpha_layer[layer_idx[n]]
              +alpha_row_layer[row_idx[n],layer_idx[n]]+alpha_col_layer[col_idx[n],layer_idx[n]];
    if(diag_cell_idx[n]>0) log_mu[n]+=alpha_diag[diag_cell_idx[n]];
    else log_mu[n]+=beta[layer_idx[n]]*psi[cell_idx[n]];
  }
}
model {
  alpha~normal(0,10); alpha_row_raw~normal(0,5); alpha_col_raw~normal(0,5);
  alpha_layer_raw~normal(0,5); alpha_diag_raw~normal(0,5);
  to_vector(alpha_row_layer_raw)~normal(0,2); to_vector(alpha_col_layer_raw)~normal(0,2);
  psi_raw~normal(0,2); beta_raw~normal(0,0.5);
  y~poisson_log(log_mu);
}
generated quantities {
  vector[N] expected; for(n in 1:N) expected[n]=exp(log_mu[n]);
}
'

stan_data <- list(N=75,I=5,J=5,K=3,IJ=25,N_diag=15,
                  row_idx=d_xie$i_O,col_idx=d_xie$i_D,layer_idx=d_xie$i_C,
                  cell_idx=d_xie$cell_OD,diag_cell_idx=d_xie$diag_cell_idx,y=d_xie$Freq)

cat("\n========== Stan MCMC ==========\n")
model <- cmdstan_model(write_stan_file(stan_code))
fit <- model$sample(data=stan_data, seed=123, chains=4, parallel_chains=4,
                    iter_warmup=1000, iter_sampling=2000, refresh=500)

# Results
cat("\n========== Results ==========\n")
fit$summary("beta")

stan_beta <- fit$summary("beta")$mean
stan_beta_sd <- fit$summary("beta")$sd

cat("\nComparison:\n")
comparison <- tibble(
  Country = c("US", "Britain", "Japan"),
  gnm = gnm_beta,
  Stan_mean = stan_beta,
  Stan_sd = stan_beta_sd,
  LEM = c(1.0000, 1.0398, 0.7989)
)
print(comparison)

# Expected frequencies
expected_stan <- fit$summary("expected")$mean
expected_gnm <- fitted(fit_gnm)
cat("\nExpected frequencies correlation:", cor(expected_gnm, expected_stan), "\n")

# Diagonal check
cat("\nDiagonal cells (first 5):\n")
d_xie |>
  mutate(gnm = expected_gnm, stan = expected_stan) |>
  filter(O == D) |>
  select(O, D, C, Freq, gnm, stan) |>
  mutate(gnm_diff = Freq - gnm, stan_diff = Freq - stan) |>
  head(5) |>
  print()
