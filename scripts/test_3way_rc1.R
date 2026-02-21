# Test 3-way RC(1) model: gnm vs Stan comparison
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

# Index variables
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
               family = poisson,
               data = tab5.4,
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

# ========== Stan Model ==========
cat("\n========== Stan Model 7 ==========\n")

stan_code <- '
data {
  int<lower=1> N_occ;
  int<lower=1> N_edu;
  int<lower=1> N_inc;
  array[N_occ, N_edu, N_inc] int<lower=0> Y;
}

parameters {
  vector[N_edu] alpha;
  vector[N_occ] beta;
  vector[N_inc] gamma;

  vector[N_edu - 1] mu_free;
  vector[N_occ - 1] nu_free;
  vector[N_inc - 1] eta_free;

  real<lower=0> phi_O;
  real<lower=0> phi_EI;
}

transformed parameters {
  vector[N_edu] mu_raw;
  vector[N_occ] nu_raw;
  vector[N_inc] eta_raw;
  vector[N_edu] mu;
  vector[N_occ] nu;
  vector[N_inc] eta;

  mu_raw[1] = 0;
  mu_raw[2:N_edu] = mu_free;
  nu_raw[1] = 0;
  nu_raw[2:N_occ] = nu_free;
  eta_raw[1] = 0;
  eta_raw[2:N_inc] = eta_free;

  mu = (mu_raw - mean(mu_raw)) / sd(mu_raw);
  nu = (nu_raw - mean(nu_raw)) / sd(nu_raw);
  eta = (eta_raw - mean(eta_raw)) / sd(eta_raw);
}

model {
  alpha ~ normal(0, 3);
  beta ~ normal(0, 3);
  gamma ~ normal(0, 3);
  mu_free ~ normal(0, 2);
  nu_free ~ normal(0, 2);
  eta_free ~ normal(0, 2);
  phi_O ~ lognormal(0, 0.5);
  phi_EI ~ lognormal(0, 0.5);

  for (j in 1:N_occ) {
    for (i in 1:N_edu) {
      for (k in 1:N_inc) {
        real log_lambda = alpha[i] + beta[j] + gamma[k]
                        + phi_O * nu[j] * (mu[i] + eta[k])
                        + phi_EI * mu[i] * eta[k];
        Y[j, i, k] ~ poisson_log(log_lambda);
      }
    }
  }
}

generated quantities {
  array[N_occ, N_edu, N_inc] real expected;
  for (j in 1:N_occ) {
    for (i in 1:N_edu) {
      for (k in 1:N_inc) {
        real log_lambda = alpha[i] + beta[j] + gamma[k]
                        + phi_O * nu[j] * (mu[i] + eta[k])
                        + phi_EI * mu[i] * eta[k];
        expected[j, i, k] = exp(log_lambda);
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

model <- cmdstan_model(write_stan_file(stan_code))

# First try MLE
cat("\n--- Stan MLE ---\n")
fit_mle <- model$optimize(data = stan_data, seed = 123, algorithm = "lbfgs",
                          iter = 10000, tol_rel_grad = 1e-12)

mle_summary <- fit_mle$summary()
nu_stan_mle <- mle_summary |> filter(grepl("^nu\\[", variable)) |> pull(estimate)
cat("Stan MLE occupation scores:\n")
print(nu_stan_mle)

cat("\nCorrelation (gnm vs Stan MLE):", cor(mu_gnm_occ$qvframe[, 1], nu_stan_mle), "\n")

# Now try MCMC
cat("\n--- Stan MCMC ---\n")
fit_mcmc <- model$sample(data = stan_data, seed = 123, chains = 4,
                         parallel_chains = 4, iter_warmup = 1000,
                         iter_sampling = 2000, refresh = 500)

nu_stan_mcmc <- fit_mcmc$summary("nu")$mean
cat("Stan MCMC occupation scores (posterior mean):\n")
print(nu_stan_mcmc)

cat("\nCorrelation (gnm vs Stan MCMC):", cor(mu_gnm_occ$qvframe[, 1], nu_stan_mcmc), "\n")

# Detailed comparison
cat("\n========== Score Comparison ==========\n")
occ_labels <- c("1:Prof", "2:Mgr", "3:Sales", "4:Clerical", "5:Craft",
                "6:Oper", "7:Serv", "8:Labor", "9:Farm", "10:FarmLab",
                "11:Milit", "12:NoOcc")

comparison <- tibble(
  category = 1:12,
  label = occ_labels,
  gnm = mu_gnm_occ$qvframe[, 1],
  stan_mle = nu_stan_mle,
  stan_mcmc = nu_stan_mcmc
)
print(comparison)

# Check phi parameters
cat("\n========== Phi Parameters ==========\n")
phi_mle <- mle_summary |> filter(grepl("^phi", variable))
print(phi_mle)

phi_mcmc <- fit_mcmc$summary(c("phi_O", "phi_EI"))
print(phi_mcmc)
