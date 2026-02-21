# Model 7 自由度チェック
library(tidyverse)
library(gnm)

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

cat("========== gnm Model 7 ==========\n")
set.seed(123)
m7 <- gnm(Freq ~ educ + occ + income +
            Mult(1, occ, educ + income) +
            Mult(1, educ, income),
          data = tab5.4,
          family = poisson,
          tolerance = 1e-12)

cat("L²:", deviance(m7), "\n")
cat("df.residual:", df.residual(m7), "\n")
cat("Total cells:", nrow(tab5.4), "\n")
cat("Estimated parameters (192 - df):", nrow(tab5.4) - df.residual(m7), "\n")

# Count coefficients
all_coef <- coef(m7)
cat("\nTotal coefficients:", length(all_coef), "\n")
cat("NA coefficients:", sum(is.na(all_coef)), "\n")
cat("Non-NA coefficients:", sum(!is.na(all_coef)), "\n")

cat("\n========== Coefficient breakdown ==========\n")
# Main effects
cat("Main effects (intercept, educ, occ, income):\n")
main_names <- grep("^(\\(Intercept\\)|educ|occ|income)[0-9]*$", names(all_coef), value = TRUE)
cat("  Count:", length(main_names), "\n")

# Mult term 1
cat("\nMult term 1 (occ, educ + income):\n")
mult1_phi <- grep("Mult\\(\\., occ", names(all_coef), value = TRUE)
mult1_occ <- grep("Mult\\(1, \\., educ \\+ income\\)", names(all_coef), value = TRUE)
mult1_educ <- grep("Mult\\(1, occ, \\. \\+ income\\)", names(all_coef), value = TRUE)
mult1_inc <- grep("Mult\\(1, occ, educ \\+ \\.\\)", names(all_coef), value = TRUE)
cat("  phi1:", length(mult1_phi), "\n")
cat("  occ scores:", length(mult1_occ), "\n")
cat("  educ scores:", length(mult1_educ), "\n")
cat("  income scores:", length(mult1_inc), "\n")
cat("  Subtotal term 1:", length(mult1_phi) + length(mult1_occ) + length(mult1_educ) + length(mult1_inc), "\n")

# Mult term 2
cat("\nMult term 2 (educ, income):\n")
mult2_phi <- grep("Mult\\(\\., educ, income\\)", names(all_coef), value = TRUE)
mult2_educ <- grep("Mult\\(1, \\., income\\)", names(all_coef), value = TRUE)
mult2_inc <- grep("Mult\\(1, educ, \\.\\)", names(all_coef), value = TRUE)
cat("  phi2:", length(mult2_phi), "\n")
cat("  educ scores:", length(mult2_educ), "\n")
cat("  income scores:", length(mult2_inc), "\n")
cat("  Subtotal term 2:", length(mult2_phi) + length(mult2_educ) + length(mult2_inc), "\n")

cat("\n========== Stan Model 7 parameters ==========\n")
# Stan Model 7 as currently implemented
stan_params <- c(
  "alpha0 (intercept)" = 1,
  "alpha_edu" = 3,
  "alpha_occ" = 11,
  "alpha_inc" = 3,
  "phi1" = 1,
  "nu (occ scores)" = 12,
  "mu1 (educ scores term1)" = 4,
  "eta1 (income scores term1)" = 3,
  "phi2" = 1,
  "mu2 (educ scores term2)" = 4,
  "eta2 (income scores term2)" = 4
)

for (nm in names(stan_params)) {
  cat(sprintf("  %s: %d\n", nm, stan_params[nm]))
}
cat("  Total Stan params:", sum(stan_params), "\n")

cat("\n========== 識別制約の問題 ==========\n")
cat("Mult項には位置・尺度の識別制約が必要:\n")
cat("  φ × a × b は φc × (a/c) × b や φ × a × (b*d) / d と同等\n")
cat("\ngnmの制約方式:\n")
cat("  - 各スコアベクトルに1つの値を固定（通常 first = 0）\n")
cat("  - φに符号・尺度を吸収\n")
cat("\nStan現行モデルの問題点:\n")
cat("  - スコアに明示的な制約なし → 識別不定\n")
cat("  - 事前分布で弱く制約しているが、dfが一致しない可能性\n")

cat("\n========== 必要な修正 ==========\n")
cat("gnm Model 7と完全に一致させるには:\n")
cat("  1. gnmと同じ制約方式を採用\n")
cat("  2. または、dfが一致することを確認してモデルを調整\n")
