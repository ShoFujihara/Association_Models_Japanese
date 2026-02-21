# Test 3-way RC(1) model - Fixed to match gnm Model 7 structure
# Key fix: Separate education/income scores for each Mult term
library(tidyverse)
library(gnm)
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

tab5.4 <- data.frame(Freq = Freq, occ = occ, educ = educ, income = income)

# ========== gnm Model 7 ==========
cat("\n========== gnm Model 7 ==========\n")
set.seed(123)
fit_gnm <- gnm(Freq ~ educ + occ + income +
                 Mult(1, occ, educ + income) +
                 Mult(1, educ, income),
               data = tab5.4,
               family = poisson,
               tolerance = 1e-12)

cat("gnm L²:", deviance(fit_gnm), "df:", df.residual(fit_gnm), "\n")

# Extract gnm occupation scores
mu_gnm_occ <- getContrasts(model = fit_gnm,
                           set = pickCoef(fit_gnm, "[.]occ"),
                           ref = "mean",
                           scaleRef = "mean",
                           scaleWeights = "unit")
cat("\ngnm occupation scores:\n")
print(mu_gnm_occ$qvframe[, 1])

# ========== Stan Model (Fixed) ==========
cat("\n========== Stan Model 7 (Fixed) ==========\n")

# Model structure matching gnm:
# log μ = main effects + φ₁ × ν_occ × (μ1_educ + η1_income) + φ₂ × μ2_educ × η2_income
# Where: mu1, eta1 are for Mult term 1; mu2, eta2 are for Mult term 2 (SEPARATE scores!)

stan_code_fixed <- '
data {
  int<lower=1> N_occ;
  int<lower=1> N_edu;
  int<lower=1> N_inc;
  array[N_occ, N_edu, N_inc] int<lower=0> Y;
}

parameters {
  // Main effects
  vector[N_edu] alpha;    // education main effect
  vector[N_occ] beta;     // occupation main effect
  vector[N_inc] gamma;    // income main effect

  // Mult term 1: Mult(1, occ, educ + income)
  // Occupation scores (shared, this is the "consistent" part)
  vector[N_occ - 1] nu_free;
  // Education scores for term 1
  vector[N_edu - 1] mu1_free;
  // Income scores for term 1
  vector[N_inc - 1] eta1_free;

  // Mult term 2: Mult(1, educ, income)
  // Education scores for term 2 (SEPARATE from term 1)
  vector[N_edu - 1] mu2_free;
  // Income scores for term 2 (SEPARATE from term 1)
  vector[N_inc - 1] eta2_free;

  // Association strength parameters
  real<lower=0> phi1;   // Mult term 1 strength
  real<lower=0> phi2;   // Mult term 2 strength
}

transformed parameters {
  // Occupation scores (standardized)
  vector[N_occ] nu_raw;
  vector[N_occ] nu;
  nu_raw[1] = 0;
  nu_raw[2:N_occ] = nu_free;
  nu = (nu_raw - mean(nu_raw)) / sd(nu_raw);

  // Education scores for Mult term 1
  vector[N_edu] mu1_raw;
  vector[N_edu] mu1;
  mu1_raw[1] = 0;
  mu1_raw[2:N_edu] = mu1_free;
  mu1 = (mu1_raw - mean(mu1_raw)) / sd(mu1_raw);

  // Income scores for Mult term 1
  vector[N_inc] eta1_raw;
  vector[N_inc] eta1;
  eta1_raw[1] = 0;
  eta1_raw[2:N_inc] = eta1_free;
  eta1 = (eta1_raw - mean(eta1_raw)) / sd(eta1_raw);

  // Education scores for Mult term 2
  vector[N_edu] mu2_raw;
  vector[N_edu] mu2;
  mu2_raw[1] = 0;
  mu2_raw[2:N_edu] = mu2_free;
  mu2 = (mu2_raw - mean(mu2_raw)) / sd(mu2_raw);

  // Income scores for Mult term 2
  vector[N_inc] eta2_raw;
  vector[N_inc] eta2;
  eta2_raw[1] = 0;
  eta2_raw[2:N_inc] = eta2_free;
  eta2 = (eta2_raw - mean(eta2_raw)) / sd(eta2_raw);
}

model {
  // Priors
  alpha ~ normal(0, 3);
  beta ~ normal(0, 3);
  gamma ~ normal(0, 3);

  nu_free ~ normal(0, 2);
  mu1_free ~ normal(0, 2);
  eta1_free ~ normal(0, 2);
  mu2_free ~ normal(0, 2);
  eta2_free ~ normal(0, 2);

  phi1 ~ lognormal(0, 0.5);
  phi2 ~ lognormal(0, 0.5);

  // Likelihood
  // Model: log μ = main + φ₁ × ν[occ] × (μ1[edu] + η1[inc]) + φ₂ × μ2[edu] × η2[inc]
  for (j in 1:N_occ) {
    for (i in 1:N_edu) {
      for (k in 1:N_inc) {
        real log_lambda = alpha[i] + beta[j] + gamma[k]
                        + phi1 * nu[j] * (mu1[i] + eta1[k])
                        + phi2 * mu2[i] * eta2[k];
        Y[j, i, k] ~ poisson_log(log_lambda);
      }
    }
  }
}

generated quantities {
  array[N_occ, N_edu, N_inc] real expected;
  real G2 = 0;

  for (j in 1:N_occ) {
    for (i in 1:N_edu) {
      for (k in 1:N_inc) {
        real log_lambda = alpha[i] + beta[j] + gamma[k]
                        + phi1 * nu[j] * (mu1[i] + eta1[k])
                        + phi2 * mu2[i] * eta2[k];
        expected[j, i, k] = exp(log_lambda);
        if (Y[j, i, k] > 0) {
          G2 += 2 * Y[j, i, k] * (log(Y[j, i, k]) - log_lambda);
        }
      }
    }
  }
}
'

# Create 3D array: Y[occ, edu, inc]
Y_array <- array(0, dim = c(12, 4, 4))
for (n in 1:192) {
  Y_array[as.integer(occ[n]), as.integer(educ[n]), as.integer(income[n])] <- Freq[n]
}

stan_data <- list(N_occ = 12, N_edu = 4, N_inc = 4, Y = Y_array)

model <- cmdstan_model(write_stan_file(stan_code_fixed))

# Stan MLE
cat("\n--- Stan MLE ---\n")
fit_mle <- model$optimize(data = stan_data, seed = 123, algorithm = "lbfgs",
                          iter = 10000, tol_rel_grad = 1e-12)

mle_summary <- fit_mle$summary()
nu_stan_mle <- mle_summary |> filter(grepl("^nu\\[", variable)) |> pull(estimate)
G2_stan <- mle_summary |> filter(variable == "G2") |> pull(estimate)

cat("Stan MLE occupation scores:\n")
print(nu_stan_mle)
cat("\nStan MLE L²:", G2_stan, "\n")

cat("\nCorrelation (gnm vs Stan MLE):", cor(mu_gnm_occ$qvframe[, 1], nu_stan_mle), "\n")

# Detailed comparison
cat("\n========== Score Comparison ==========\n")
occ_labels <- c("1:Prof", "2:Mgr", "3:Sales", "4:Clerical", "5:Craft",
                "6:Oper", "7:Serv", "8:Labor", "9:Farm", "10:FarmLab",
                "11:Milit", "12:NoOcc")

comparison <- tibble(
  category = 1:12,
  label = occ_labels,
  gnm = mu_gnm_occ$qvframe[, 1],
  stan_mle = nu_stan_mle
)
print(comparison)

# Expected frequencies comparison
cat("\n========== Expected Frequencies ==========\n")
expected_gnm <- fitted(fit_gnm)
expected_stan_array <- mle_summary |>
  filter(grepl("^expected\\[", variable)) |>
  pull(estimate)

# Flatten Stan expected to match gnm order
expected_stan <- numeric(192)
idx <- 1
for (n in 1:192) {
  j <- as.integer(occ[n])
  i <- as.integer(educ[n])
  k <- as.integer(income[n])
  # Stan array index: (j-1)*4*4 + (i-1)*4 + k
  stan_idx <- (j-1)*16 + (i-1)*4 + k
  expected_stan[n] <- expected_stan_array[stan_idx]
}

cat("Expected frequencies correlation:", cor(expected_gnm, expected_stan), "\n")
cat("Max absolute difference:", max(abs(expected_gnm - expected_stan)), "\n")
