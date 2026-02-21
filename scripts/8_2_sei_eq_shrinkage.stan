// rc_3way.stan
// Notation: edu = alpha/mu (i), occ = beta/nu (j), inc = gamma/eta (k)
// Data array: Y[occ, edu, inc] = Y[j, i, k]

data {
  int<lower=1> N_occ;    // Number of occupations (J)
  int<lower=1> N_edu;    // Number of education categories (I)
  int<lower=1> N_inc;    // Number of income categories (K)
  array[N_occ, N_edu, N_inc] int<lower=0> Y;
}

parameters {
  // Main effects
  vector[N_edu] alpha;   // Education main effect (α_i)
  vector[N_occ] beta;    // Occupation main effect (β_j)
  vector[N_inc] gamma;   // Income main effect (γ_k)
  
  // Scores (raw parameters)
  ordered[N_edu] mu_ordered;      // Education score (μ_i)
  vector[N_occ] nu_raw;           // Occupation score (ν_j)
  ordered[N_inc] eta_ordered;     // Income score (η_k)

  // Association parameters
  real<lower=0> phi_EO;  // Education-Occupation (φ^EO)
  real<lower=0> phi_OI;  // Occupation-Income (φ^OI)
  real<lower=0> phi_EI;  // Education-Income (φ^EI)
}

transformed parameters {
  vector[N_edu] mu;
  vector[N_occ] nu;
  vector[N_inc] eta;
  
  // Standardization (mean=0, sd=1)
  mu = (mu_ordered - mean(mu_ordered)) / sd(mu_ordered);
  nu = (nu_raw - mean(nu_raw)) / sd(nu_raw);
  eta = (eta_ordered - mean(eta_ordered)) / sd(eta_ordered);
}

model {
  // Priors
  alpha ~ normal(0, 3);
  beta ~ normal(0, 3);
  gamma ~ normal(0, 3);
  
  mu_ordered ~ normal(0, 1);
  nu_raw ~ student_t(4, 0, 1);  // 階層的 + 頑健
  eta_ordered ~ normal(0, 1);
  
  phi_EO ~ lognormal(0, 0.5);
  phi_OI ~ lognormal(0, 0.5);
  phi_EI ~ lognormal(0, 0.5);
  
  // Likelihood: Y[occ, edu, inc] = Y[j, i, k]
  for (j in 1:N_occ) {
    for (i in 1:N_edu) {
      for (k in 1:N_inc) {
        real log_lambda = alpha[i] + beta[j] + gamma[k]
                        + phi_EO * mu[i] * nu[j]
                        + phi_OI * nu[j] * eta[k]
                        + phi_EI * mu[i] * eta[k];
        Y[j, i, k] ~ poisson_log(log_lambda);
      }
    }
  }
}

generated quantities {
  // 偏差値スコア（平均50、SD10）
  vector[N_occ] SEI;
  real nu_mean = mean(nu);
  real nu_sd = sd(nu);
  
  // ゼロ除算を防ぐ
  if (nu_sd > 0.0001) {
    SEI = 50 + 10 * (nu - nu_mean) / nu_sd;
  } else {
    SEI = rep_vector(50, N_occ);
  }
}
