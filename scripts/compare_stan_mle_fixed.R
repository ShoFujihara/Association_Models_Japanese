# Stan MLE with properly saturated diagonal cells
library(tidyverse)
library(gnm)
library(cmdstanr)

# Data preparation
Freq_yamaguchi <- c(
  1275, 364, 274, 272, 17,
  1055, 597, 394, 443, 31,
  1043, 587, 1045, 951, 47,
  1159, 791, 1323, 2046, 52,
  666, 496, 1031, 1632, 646,
  474, 129, 87, 124, 11,
  300, 218, 171, 220, 8,
  438, 254, 669, 703, 16,
  601, 388, 932, 1789, 37,
  76, 56, 125, 295, 191,
  127, 101, 24, 30, 12,
  86, 207, 64, 61, 13,
  43, 73, 122, 60, 13,
  35, 51, 62, 66, 11,
  109, 206, 184, 253, 325
)

d_xie <- tibble(
  O = gl(5, k = 5, length = 3 * 5 * 5),
  D = gl(5, k = 1, length = 3 * 5 * 5),
  C = gl(3, k = 5 * 5, length = 3 * 5 * 5),
  Freq = Freq_yamaguchi
) |>
  mutate(
    Diag = case_when(O == D ~ O, .default = "0") |> factor(),
    is_diag = O == D,
    i_O = as.integer(O),
    i_D = as.integer(D),
    i_C = as.integer(C),
    cell_OD = (i_O - 1) * 5 + i_D,
    # Diagonal cell index: 1-15 for diagonal, 0 for off-diagonal
    diag_cell_idx = if_else(O == D, (as.integer(C) - 1) * 5 + as.integer(O), 0L)
  )

# gnm
cat("\n========== gnm ==========\n")
set.seed(12345)
fit_gnm <- gnm(Freq ~ O*C + D*C + Diag*C + Mult(Exp(C), O*D),
               family = poisson, data = d_xie, tolerance = 1e-12)
expected_gnm <- fitted(fit_gnm)

# Extract gnm scaling parameters
coef_gnm <- coef(fit_gnm)
exp_c_coefs <- coef_gnm[grep("^Mult.*Exp.*\\.C", names(coef_gnm))]
psi_gnm <- exp(exp_c_coefs - exp_c_coefs[1])
names(psi_gnm) <- c("US", "Britain", "Japan")
cat("gnm scaling (psi):\n")
print(psi_gnm)
cat("gnm L²:", round(deviance(fit_gnm), 4), "df:", df.residual(fit_gnm), "\n")

# Stan with FULLY saturated diagonal (15 independent parameters)
cat("\n========== Stan MLE (fixed diagonal) ==========\n")

stan_code_fixed <- "
data {
  int<lower=1> N;
  int<lower=1> I;
  int<lower=1> J;
  int<lower=1> K;
  int<lower=1> IJ;
  int<lower=1> N_diag;  // Number of diagonal cells (I * K = 15)
  array[N] int<lower=1> row_idx;
  array[N] int<lower=1> col_idx;
  array[N] int<lower=1> layer_idx;
  array[N] int<lower=1> cell_idx;
  array[N] int<lower=0, upper=N_diag> diag_cell_idx;  // 0 for off-diagonal, 1-15 for diagonal
  array[N] int<lower=0> y;
}

parameters {
  real alpha;
  vector[I-1] alpha_row_raw;
  vector[J-1] alpha_col_raw;
  vector[K-1] alpha_layer_raw;
  matrix[I-1, K-1] alpha_row_layer_raw;
  matrix[J-1, K-1] alpha_col_layer_raw;
  vector[N_diag-1] alpha_diag_raw;  // Fully saturated diagonal (14 free params)
  vector[IJ-1] psi_raw;
  vector[K-1] beta_raw;
}

transformed parameters {
  vector[I] alpha_row;
  vector[J] alpha_col;
  vector[K] alpha_layer;
  matrix[I, K] alpha_row_layer;
  matrix[J, K] alpha_col_layer;
  vector[N_diag] alpha_diag;  // Full diagonal saturation
  vector[IJ] psi;
  vector[K] beta;
  vector[N] log_mu;

  // Main effects with corner constraints
  alpha_row[1] = 0;
  alpha_col[1] = 0;
  alpha_layer[1] = 0;
  for (i in 2:I) alpha_row[i] = alpha_row_raw[i-1];
  for (j in 2:J) alpha_col[j] = alpha_col_raw[j-1];
  for (k in 2:K) alpha_layer[k] = alpha_layer_raw[k-1];

  // Interactions
  for (k in 1:K) { alpha_row_layer[1,k] = 0; alpha_col_layer[1,k] = 0; }
  for (i in 1:I) alpha_row_layer[i,1] = 0;
  for (j in 1:J) alpha_col_layer[j,1] = 0;
  for (i in 2:I) for (k in 2:K) alpha_row_layer[i,k] = alpha_row_layer_raw[i-1,k-1];
  for (j in 2:J) for (k in 2:K) alpha_col_layer[j,k] = alpha_col_layer_raw[j-1,k-1];

  // Fully saturated diagonal (15 cells, 14 free parameters)
  alpha_diag[1] = 0;
  for (d in 2:N_diag) alpha_diag[d] = alpha_diag_raw[d-1];

  // OD association pattern
  psi[1] = 0;
  for (ij in 2:IJ) psi[ij] = psi_raw[ij-1];

  // Layer scaling (Unidiff)
  beta[1] = 1;
  for (k in 2:K) beta[k] = exp(beta_raw[k-1]);

  // Expected frequencies
  for (n in 1:N) {
    log_mu[n] = alpha + alpha_row[row_idx[n]] + alpha_col[col_idx[n]]
                + alpha_layer[layer_idx[n]]
                + alpha_row_layer[row_idx[n], layer_idx[n]]
                + alpha_col_layer[col_idx[n], layer_idx[n]];

    if (diag_cell_idx[n] > 0) {
      // Diagonal cell: use saturated diagonal parameter
      log_mu[n] += alpha_diag[diag_cell_idx[n]];
    } else {
      // Off-diagonal cell: apply Unidiff
      log_mu[n] += beta[layer_idx[n]] * psi[cell_idx[n]];
    }
  }
}

model {
  y ~ poisson_log(log_mu);
}

generated quantities {
  vector[N] expected;
  for (n in 1:N) expected[n] = exp(log_mu[n]);
}
"

stan_data <- list(
  N = nrow(d_xie),
  I = 5, J = 5, K = 3, IJ = 25,
  N_diag = 15,  # 5 diagonal positions × 3 countries
  row_idx = d_xie$i_O,
  col_idx = d_xie$i_D,
  layer_idx = d_xie$i_C,
  cell_idx = d_xie$cell_OD,
  diag_cell_idx = d_xie$diag_cell_idx,
  y = d_xie$Freq
)

model <- cmdstan_model(write_stan_file(stan_code_fixed))
fit_mle <- model$optimize(data = stan_data, seed = 123, algorithm = "lbfgs",
                          iter = 10000, tol_rel_grad = 1e-12)

mle_summary <- fit_mle$summary()
beta_stan <- mle_summary |> filter(grepl("^beta\\[", variable)) |> pull(estimate)
names(beta_stan) <- c("US", "Britain", "Japan")
cat("Stan MLE scaling (beta):\n")
print(beta_stan)

expected_stan <- mle_summary |> filter(grepl("^expected\\[", variable)) |> pull(estimate)
G2_stan <- 2 * sum(d_xie$Freq * log(d_xie$Freq / expected_stan), na.rm = TRUE)
cat("Stan MLE L²:", round(G2_stan, 4), "\n")

# Comparison
cat("\n========== Comparison ==========\n")
comparison <- tibble(
  Country = c("US", "Britain", "Japan"),
  gnm = psi_gnm,
  Stan_MLE = beta_stan,
  LEM = c(1.0000, 1.0398, 0.7989)
)
print(comparison)

cat("\nExpected frequencies correlation:", cor(expected_gnm, expected_stan), "\n")

# Check diagonal cells
cat("\n========== Diagonal cells check ==========\n")
d_xie |>
  filter(is_diag) |>
  mutate(
    gnm = expected_gnm[row_number()],
    stan = expected_stan[row_number()]
  ) |>
  select(O, D, C, Freq, gnm, stan) |>
  mutate(
    gnm_match = abs(Freq - gnm) < 0.01,
    stan_match = abs(Freq - stan) < 0.01
  ) |>
  print(n = 15)
