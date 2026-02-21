# Compare expected frequencies between gnm and Stan MLE
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
    diag_idx = if_else(O == D, as.integer(O), 0L)
  )

# gnm
set.seed(12345)
fit_gnm <- gnm(Freq ~ O*C + D*C + Diag*C + Mult(Exp(C), O*D),
               family = poisson, data = d_xie, tolerance = 1e-12)
expected_gnm <- fitted(fit_gnm)

# Stan (read from previous run or re-run)
stan_code_mle <- "
data {
  int<lower=1> N;
  int<lower=1> I;
  int<lower=1> J;
  int<lower=1> K;
  int<lower=1> IJ;
  int<lower=1> D_max;
  array[N] int<lower=1> row_idx;
  array[N] int<lower=1> col_idx;
  array[N] int<lower=1> layer_idx;
  array[N] int<lower=1> cell_idx;
  array[N] int<lower=0, upper=D_max> diag_idx;
  array[N] int<lower=0> y;
}

parameters {
  real alpha;
  vector[I-1] alpha_row_raw;
  vector[J-1] alpha_col_raw;
  vector[K-1] alpha_layer_raw;
  matrix[I-1, K-1] alpha_row_layer_raw;
  matrix[J-1, K-1] alpha_col_layer_raw;
  matrix[D_max-1, K-1] alpha_diag_layer_raw;
  vector[D_max-1] alpha_diag_raw;
  vector[IJ-1] psi_raw;
  vector[K-1] beta_raw;
}

transformed parameters {
  vector[I] alpha_row;
  vector[J] alpha_col;
  vector[K] alpha_layer;
  matrix[I, K] alpha_row_layer;
  matrix[J, K] alpha_col_layer;
  vector[D_max] alpha_diag;
  matrix[D_max, K] alpha_diag_layer;
  vector[IJ] psi;
  vector[K] beta;
  vector[N] log_mu;

  alpha_row[1] = 0; alpha_col[1] = 0; alpha_layer[1] = 0; alpha_diag[1] = 0;
  for (i in 2:I) alpha_row[i] = alpha_row_raw[i-1];
  for (j in 2:J) alpha_col[j] = alpha_col_raw[j-1];
  for (k in 2:K) alpha_layer[k] = alpha_layer_raw[k-1];
  for (d in 2:D_max) alpha_diag[d] = alpha_diag_raw[d-1];

  for (k in 1:K) { alpha_row_layer[1,k]=0; alpha_col_layer[1,k]=0; alpha_diag_layer[1,k]=0; }
  for (i in 1:I) alpha_row_layer[i,1] = 0;
  for (j in 1:J) alpha_col_layer[j,1] = 0;
  for (d in 1:D_max) alpha_diag_layer[d,1] = 0;
  for (i in 2:I) for (k in 2:K) alpha_row_layer[i,k] = alpha_row_layer_raw[i-1,k-1];
  for (j in 2:J) for (k in 2:K) alpha_col_layer[j,k] = alpha_col_layer_raw[j-1,k-1];
  for (d in 2:D_max) for (k in 2:K) alpha_diag_layer[d,k] = alpha_diag_layer_raw[d-1,k-1];

  psi[1] = 0;
  for (ij in 2:IJ) psi[ij] = psi_raw[ij-1];

  beta[1] = 1;
  for (k in 2:K) beta[k] = exp(beta_raw[k-1]);

  for (n in 1:N) {
    real diag_effect = 0;
    if (diag_idx[n] > 0) {
      diag_effect = alpha_diag[diag_idx[n]] + alpha_diag_layer[diag_idx[n], layer_idx[n]];
    }
    log_mu[n] = alpha + alpha_row[row_idx[n]] + alpha_col[col_idx[n]]
                + alpha_layer[layer_idx[n]]
                + alpha_row_layer[row_idx[n], layer_idx[n]]
                + alpha_col_layer[col_idx[n], layer_idx[n]]
                + diag_effect
                + beta[layer_idx[n]] * psi[cell_idx[n]];
  }
}

model { y ~ poisson_log(log_mu); }

generated quantities {
  vector[N] expected;
  for (n in 1:N) expected[n] = exp(log_mu[n]);
}
"

stan_data <- list(
  N = nrow(d_xie), I = 5, J = 5, K = 3, IJ = 25, D_max = 5,
  row_idx = d_xie$i_O, col_idx = d_xie$i_D, layer_idx = d_xie$i_C,
  cell_idx = d_xie$cell_OD, diag_idx = d_xie$diag_idx, y = d_xie$Freq
)

model <- cmdstan_model(write_stan_file(stan_code_mle))
fit_mle <- model$optimize(data = stan_data, seed = 123, algorithm = "lbfgs",
                          iter = 10000, tol_rel_grad = 1e-12)
mle_summary <- fit_mle$summary()
expected_stan <- mle_summary |> filter(grepl("^expected\\[", variable)) |> pull(estimate)

# Comparison
comparison <- d_xie |>
  mutate(
    gnm = expected_gnm,
    stan = expected_stan,
    diff = stan - gnm,
    pct_diff = (stan - gnm) / gnm * 100,
    abs_diff = abs(diff)
  ) |>
  select(O, D, C, is_diag, Freq, gnm, stan, diff, pct_diff)

# Show all cells sorted by difference
cat("\n========== All cells sorted by absolute difference ==========\n")
comparison |>
  arrange(desc(abs(diff))) |>
  print(n = 75)

# Summary by diagonal vs off-diagonal
cat("\n========== Summary by diagonal status ==========\n")
comparison |>
  group_by(is_diag) |>
  summarise(
    n = n(),
    mean_diff = mean(diff),
    sd_diff = sd(diff),
    max_abs_diff = max(abs(diff)),
    mean_pct_diff = mean(pct_diff)
  ) |>
  print()

# Summary by country
cat("\n========== Summary by country ==========\n")
comparison |>
  group_by(C) |>
  summarise(
    n = n(),
    mean_diff = mean(diff),
    sd_diff = sd(diff),
    max_abs_diff = max(abs(diff)),
    mean_pct_diff = mean(pct_diff)
  ) |>
  print()

# Summary by country and diagonal status
cat("\n========== Summary by country and diagonal ==========\n")
comparison |>
  group_by(C, is_diag) |>
  summarise(
    n = n(),
    mean_diff = mean(diff),
    max_abs_diff = max(abs(diff)),
    .groups = "drop"
  ) |>
  print()

# Largest differences
cat("\n========== Top 10 largest differences ==========\n")
comparison |>
  arrange(desc(abs(diff))) |>
  head(10) |>
  print()
