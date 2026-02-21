# Model 7 vs Model 8 職業スコア比較
library(tidyverse)
library(cmdstanr)

# Data (Table 5.4)
Freq <- c(1096,  1847,  1255,  925,    3321,  6123,  6830,  5524,
          1541,  3134,  3145, 3300,    1915,  4429,  7035,  9421,
          4183,  5139,  1857, 1272,    8080,  8586,  4788,  4294,
          6033,  9211,  5046, 1058,   28130, 44589, 20074,  3408,
          4354, 13430, 18670, 9821,    2250,  9075, 18286, 14358,
         14587, 31470, 16390, 3751,    8242, 17532, 12825,  3956,
          1517,  5820,  6197, 2372,     721,  2909,  4141,  2070,
          3581,  9268,  5463, 1007,    1341,  3808,  3163,   815,
          1454,  3109,  1055,  888,     563,  1909,  1018,  1051,
          3237,  3851,   377,  102,     731,   858,   247,    84,
         14882, 22182,  5363, 1136,   11650, 15818,  5524,  2122,
          6033,  3475,    63,   18,    1603,  1005,    30,    16,

          5968,  8783,  7701, 6483,    8733, 14329, 19386, 28143,
          1011,  2162,  3536, 6649,     697,  1299,  2362, 10796,
          3214,  3621,  2485, 3177,     793,  1134,  1292,  3597,
         11532, 16837,  6975, 1839,    2563,  2995,  2060,  1600,
          1009,  2719,  3521, 3409,     296,   503,   626,  1273,
          1586,  3025,  1726,  668,     245,   415,   238,   218,
           387,   941,   564,  316,      86,   138,    79,    48,
           994,  1988,   542,  145,     158,   259,   101,    56,
           171,   409,   223,  245,      65,   172,    99,   174,
           293,   290,    67,   31,      32,    62,    18,    30,
          4288,  4916,  1452,  766,     616,   794,   347,   300,
           370,   186,     3,    4,      67,    37,     5,     2)

income <- gl(4, 1, 192)
occ <- gl(12, 8, 192)
edu1 <- rep(c(1,2), each = 4, length.out = 96)
edu2 <- rep(c(3,4), each = 4, length.out = 96)
educ <- c(edu1, edu2) |> factor()

# Create 3D array
Y_array <- array(0, dim = c(12, 4, 4))
for (n in 1:192) {
  Y_array[as.integer(occ[n]), as.integer(educ[n]), as.integer(income[n])] <- Freq[n]
}
stan_data <- list(N_occ = 12, N_edu = 4, N_inc = 4, Y = Y_array)

# ========== Model 7 (2φ) ==========
cat("\n========== Model 7 (2φ: gnm Model 7相当) ==========\n")
stan_code_model7 <- '
data {
  int<lower=1> N_occ;
  int<lower=1> N_edu;
  int<lower=1> N_inc;
  array[N_occ, N_edu, N_inc] int<lower=0> Y;
}

parameters {
  real alpha0;
  vector[N_edu - 1] alpha_edu;
  vector[N_occ - 1] alpha_occ;
  vector[N_inc - 1] alpha_inc;

  real phi1;
  vector[N_occ] nu;
  vector[N_edu] mu1;
  vector[N_inc - 1] eta1;

  real phi2;
  vector[N_edu] mu2;
  vector[N_inc] eta2;
}

model {
  alpha0 ~ normal(0, 10);
  alpha_edu ~ normal(0, 5);
  alpha_occ ~ normal(0, 5);
  alpha_inc ~ normal(0, 5);
  phi1 ~ normal(0, 5);
  phi2 ~ normal(0, 2);
  nu ~ normal(0, 2);
  mu1 ~ normal(0, 2);
  eta1 ~ normal(0, 2);
  mu2 ~ normal(0, 5);
  eta2 ~ normal(0, 2);

  for (j in 1:N_occ) {
    for (i in 1:N_edu) {
      for (k in 1:N_inc) {
        real main_eff = alpha0;
        if (i > 1) main_eff += alpha_edu[i - 1];
        if (j > 1) main_eff += alpha_occ[j - 1];
        if (k > 1) main_eff += alpha_inc[k - 1];
        real eta1_k = (k == 1) ? 0.0 : eta1[k - 1];
        real mult1 = phi1 * nu[j] * (mu1[i] + eta1_k);
        real mult2 = phi2 * mu2[i] * eta2[k];
        Y[j, i, k] ~ poisson_log(main_eff + mult1 + mult2);
      }
    }
  }
}

generated quantities {
  vector[N_occ] nu_std;
  nu_std = (nu - mean(nu)) / sd(nu);
  vector[N_occ] nu_final;
  if (nu_std[1] > 0) {
    nu_final = nu_std;
  } else {
    nu_final = -nu_std;
  }
}
'

model7 <- cmdstan_model(write_stan_file(stan_code_model7))
fit7 <- model7$optimize(data = stan_data, seed = 123, algorithm = "lbfgs",
                        iter = 10000, tol_rel_grad = 1e-12)
nu_model7 <- fit7$summary() |> filter(grepl("^nu_final\\[", variable)) |> pull(estimate)
cat("Model 7 occupation scores (nu_final):\n")
print(round(nu_model7, 3))

# ========== Model 8 (3φ) ==========
cat("\n========== Model 8 (3φ: スコア全共通) ==========\n")
stan_code_model8 <- '
data {
  int<lower=1> N_occ;
  int<lower=1> N_edu;
  int<lower=1> N_inc;
  array[N_occ, N_edu, N_inc] int<lower=0> Y;
}

parameters {
  real alpha0;
  vector[N_edu - 1] alpha_edu;
  vector[N_occ - 1] alpha_occ;
  vector[N_inc - 1] alpha_inc;
  
  vector[N_edu] mu;
  vector[N_occ] nu;
  vector[N_inc] eta;
  
  real phi_EO;
  real phi_OI;
  real phi_EI;
}

model {
  alpha0 ~ normal(0, 10);
  alpha_edu ~ normal(0, 5);
  alpha_occ ~ normal(0, 5);
  alpha_inc ~ normal(0, 5);
  mu ~ normal(0, 2);
  nu ~ normal(0, 2);
  eta ~ normal(0, 2);
  phi_EO ~ normal(0, 2);
  phi_OI ~ normal(0, 2);
  phi_EI ~ normal(0, 2);

  for (j in 1:N_occ) {
    for (i in 1:N_edu) {
      for (k in 1:N_inc) {
        real main_eff = alpha0;
        if (i > 1) main_eff += alpha_edu[i - 1];
        if (j > 1) main_eff += alpha_occ[j - 1];
        if (k > 1) main_eff += alpha_inc[k - 1];
        real log_lambda = main_eff
                        + phi_EO * mu[i] * nu[j]
                        + phi_OI * nu[j] * eta[k]
                        + phi_EI * mu[i] * eta[k];
        Y[j, i, k] ~ poisson_log(log_lambda);
      }
    }
  }
}

generated quantities {
  vector[N_occ] nu_std;
  nu_std = (nu - mean(nu)) / sd(nu);
  vector[N_occ] nu_final;
  if (nu_std[1] > 0) {
    nu_final = nu_std;
  } else {
    nu_final = -nu_std;
  }
}
'

model8 <- cmdstan_model(write_stan_file(stan_code_model8))
fit8 <- model8$optimize(data = stan_data, seed = 123, algorithm = "lbfgs",
                        iter = 10000, tol_rel_grad = 1e-12)
nu_model8 <- fit8$summary() |> filter(grepl("^nu_final\\[", variable)) |> pull(estimate)
cat("Model 8 occupation scores (nu_final):\n")
print(round(nu_model8, 3))

# ========== 比較 ==========
cat("\n========== Model 7 vs Model 8 比較 ==========\n")
occ_labels <- c("専門職", "管理職", "販売", "事務", "熟練工", 
                "半熟練", "サービス", "労働者", "農業", "農業労働",
                "軍人", "無職")

comparison <- tibble(
  category = 1:12,
  label = occ_labels,
  model7 = nu_model7,
  model8 = nu_model8,
  diff = model7 - model8
)
print(comparison, n = 12)

cat("\n相関係数:", cor(nu_model7, nu_model8), "\n")
