# Test 3-way RC(1) model - Version 2: Match gnm parameterization exactly
library(tidyverse)
library(gnm)
library(cmdstanr)

# Data
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

tab5.4 <- data.frame(Freq = Freq, occ = occ, educ = educ, income = income)

# gnm Model 7
cat("\n========== gnm Model 7 ==========\n")
set.seed(123)
fit_gnm <- gnm(Freq ~ educ + occ + income +
                 Mult(1, occ, educ + income) +
                 Mult(1, educ, income),
               data = tab5.4,
               family = poisson,
               tolerance = 1e-12)

cat("gnm L²:", deviance(fit_gnm), "df:", df.residual(fit_gnm), "\n")

# Extract gnm scores
mu_gnm_occ <- getContrasts(model = fit_gnm,
                           set = pickCoef(fit_gnm, "[.]occ"),
                           ref = "mean",
                           scaleRef = "mean",
                           scaleWeights = "unit")
gnm_occ_scores <- mu_gnm_occ$qvframe[, 1]

# ========== Stan Model (matching gnm exactly) ==========
# Key changes:
# 1. phi can be negative (no lower=0)
# 2. No standardization during estimation
# 3. Corner constraints matching gnm

stan_code_v2 <- '
data {
  int<lower=1> N_occ;
  int<lower=1> N_edu;
  int<lower=1> N_inc;
  array[N_occ, N_edu, N_inc] int<lower=0> Y;
}

parameters {
  // Main effects (corner constraints: first = 0)
  real alpha0;                    // intercept
  vector[N_edu - 1] alpha_edu;    // educ2-4
  vector[N_occ - 1] alpha_occ;    // occ2-12
  vector[N_inc - 1] alpha_inc;    // income2-4

  // Mult term 1: phi1 * nu[occ] * (mu1[educ] + eta1[income])
  real phi1;                       // can be negative!
  vector[N_occ] nu;                // all 12 occ scores
  vector[N_edu] mu1;               // all 4 educ scores for term 1
  vector[N_inc - 1] eta1;          // income2-4 for term 1 (income1 = 0)

  // Mult term 2: phi2 * mu2[educ] * eta2[income]
  real phi2;                       // can be negative
  vector[N_edu] mu2;               // all 4 educ scores for term 2
  vector[N_inc] eta2;              // all 4 income scores for term 2
}

model {
  // Priors
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

  // Likelihood
  for (j in 1:N_occ) {
    for (i in 1:N_edu) {
      for (k in 1:N_inc) {
        // Main effects with corner constraints
        real main_eff = alpha0;
        if (i > 1) main_eff += alpha_edu[i - 1];
        if (j > 1) main_eff += alpha_occ[j - 1];
        if (k > 1) main_eff += alpha_inc[k - 1];

        // Mult term 1: income1 score = 0 in this term
        real eta1_k = (k == 1) ? 0.0 : eta1[k - 1];
        real mult1 = phi1 * nu[j] * (mu1[i] + eta1_k);

        // Mult term 2: all income scores estimated
        real mult2 = phi2 * mu2[i] * eta2[k];

        real log_lambda = main_eff + mult1 + mult2;
        Y[j, i, k] ~ poisson_log(log_lambda);
      }
    }
  }
}

generated quantities {
  // Standardized occupation scores (for comparison with gnm)
  vector[N_occ] nu_std;
  real nu_mean = mean(nu);
  real nu_sd_val = sd(nu);
  nu_std = (nu - nu_mean) / nu_sd_val;
}
'

# Create 3D array
Y_array <- array(0, dim = c(12, 4, 4))
for (n in 1:192) {
  Y_array[as.integer(occ[n]), as.integer(educ[n]), as.integer(income[n])] <- Freq[n]
}

stan_data <- list(N_occ = 12, N_edu = 4, N_inc = 4, Y = Y_array)

model <- cmdstan_model(write_stan_file(stan_code_v2))

# Stan MLE
cat("\n--- Stan MLE ---\n")
fit_mle <- model$optimize(data = stan_data, seed = 123, algorithm = "lbfgs",
                          iter = 10000, tol_rel_grad = 1e-12)

mle_summary <- fit_mle$summary()

# Raw occ scores
nu_raw <- mle_summary |> filter(grepl("^nu\\[", variable)) |> pull(estimate)
cat("Raw occ scores (nu):\n")
print(nu_raw)

# Standardized scores
nu_std <- mle_summary |> filter(grepl("^nu_std\\[", variable)) |> pull(estimate)
cat("\nStandardized occ scores (nu_std):\n")
print(nu_std)

# 符号反転: 専門職（カテゴリ1）が高いスコアになるよう調整
# RC(1)モデルでは符号の不定性があるため、解釈しやすい向きに揃える
if (nu_std[1] > 0) {
  # 専門職が正 = 高スコア → そのまま
  nu_final <- nu_std
  cat("\n符号調整: 不要（専門職が既に高スコア）\n")
} else {
  # 専門職が負 → 符号反転
  nu_final <- -nu_std
  cat("\n符号調整: 反転（専門職を高スコアに）\n")
}

cat("\n最終職業スコア（専門職=高）:\n")
occ_labels <- c("1:専門職", "2:管理職", "3:販売", "4:事務", "5:熟練",
                "6:半熟練", "7:サービス", "8:労働者", "9:農業", "10:農業労働",
                "11:軍", "12:無職")
names(nu_final) <- occ_labels
print(round(nu_final, 3))

# Compare with gnm
cat("\n========== Comparison ==========\n")
cat("gnm occupation scores:\n")
print(gnm_occ_scores)

cat("\nCorrelation (gnm vs Stan nu_std):", cor(gnm_occ_scores, nu_std), "\n")

# Check phi values
phi1_est <- mle_summary |> filter(variable == "phi1") |> pull(estimate)
phi2_est <- mle_summary |> filter(variable == "phi2") |> pull(estimate)
cat("\nStan phi1:", phi1_est, "(gnm: -3.852)\n")
cat("Stan phi2:", phi2_est, "(gnm: 0.470)\n")

# Calculate expected frequencies and L²
cat("\n========== Expected Frequencies & L² ==========\n")
expected_gnm <- fitted(fit_gnm)

# Manually calculate Stan expected
alpha0_est <- mle_summary |> filter(variable == "alpha0") |> pull(estimate)
alpha_edu_est <- mle_summary |> filter(grepl("^alpha_edu\\[", variable)) |> pull(estimate)
alpha_occ_est <- mle_summary |> filter(grepl("^alpha_occ\\[", variable)) |> pull(estimate)
alpha_inc_est <- mle_summary |> filter(grepl("^alpha_inc\\[", variable)) |> pull(estimate)
mu1_est <- mle_summary |> filter(grepl("^mu1\\[", variable)) |> pull(estimate)
eta1_est <- mle_summary |> filter(grepl("^eta1\\[", variable)) |> pull(estimate)
mu2_est <- mle_summary |> filter(grepl("^mu2\\[", variable)) |> pull(estimate)
eta2_est <- mle_summary |> filter(grepl("^eta2\\[", variable)) |> pull(estimate)

expected_stan <- numeric(192)
for (n in 1:192) {
  j <- as.integer(occ[n])
  i <- as.integer(educ[n])
  k <- as.integer(income[n])

  # Main effects
  main_eff <- alpha0_est
  if (i > 1) main_eff <- main_eff + alpha_edu_est[i - 1]
  if (j > 1) main_eff <- main_eff + alpha_occ_est[j - 1]
  if (k > 1) main_eff <- main_eff + alpha_inc_est[k - 1]

  # Mult term 1
  eta1_k <- if (k == 1) 0 else eta1_est[k - 1]
  mult1 <- phi1_est * nu_raw[j] * (mu1_est[i] + eta1_k)

  # Mult term 2
  mult2 <- phi2_est * mu2_est[i] * eta2_est[k]

  expected_stan[n] <- exp(main_eff + mult1 + mult2)
}

cat("Expected frequencies correlation:", cor(expected_gnm, expected_stan), "\n")
cat("Max absolute difference:", max(abs(expected_gnm - expected_stan)), "\n")

# L²
G2_stan <- 2 * sum(Freq * log(Freq / expected_stan))
cat("\ngnm L²:", deviance(fit_gnm), "\n")
cat("Stan L²:", G2_stan, "\n")
